"""后端 / 相对反射率 / 可迁移性的纯数值单元测试。

不依赖 23 GB 成果：用 32x32 合成辐亮度立方体验证 IARR 的数值逻辑、
产品标注、防"配置漂移"保护，以及后端与 UTM 带判定。
"""
import sys
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import rasterio
from rasterio.transform import from_origin

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ahsi import atmos, preflight, relative  # noqa: E402
from ahsi.config import load_config  # noqa: E402

CFG_PATH = ROOT / "config" / "pipeline.yaml"
SENSORS_PATH = ROOT / "config" / "sensors.yaml"


def _synthetic_cube(tmp: Path) -> dict:
    h = w = 32
    prof = {"driver": "GTiff", "width": w, "height": h, "count": 3, "dtype": "float32",
            "crs": "EPSG:32650", "transform": from_origin(0, 320, 10, 10), "nodata": -9999.0}
    cube, mask = tmp / "radiance.tif", tmp / "coverage.tif"
    with rasterio.open(cube, "w", **prof) as d:
        d.write(np.full((h, w), 10.0, "float32"), 1)
        d.write(np.full((h, w), 20.0, "float32"), 2)
        d.write(np.zeros((h, w), "float32"), 3)   # 全 0 -> bad band
    with rasterio.open(mask, "w", **{**prof, "count": 1, "dtype": "uint8", "nodata": 0}) as d:
        d.write(np.ones((h, w), "uint8"), 1)
    rows = [{"analysis_candidate": True, "smoke_output_band": i, "output_band": i,
             "wavelength_nm": float(wl), "fwhm_nm": 5.0, "segment": "VNIR",
             "source_band": i, "source_file": "synthetic", "gain": 1.0, "offset": 0.0}
            for i, wl in enumerate((500, 600, 700), 1)]
    return {"radiance_file": str(cube), "coverage_mask": str(mask),
            "band_records": rows, "width": w, "height": h, "scene_id": "SYNTH"}


class TestBackendRegistry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config(CFG_PATH)

    def test_default_is_dependency_free(self):
        self.assertEqual(self.cfg["atmospheric"]["method"], "iarr")
        b = atmos.detect(self.cfg)
        self.assertTrue(b["available"])
        self.assertFalse(b["requires_external_software"])
        self.assertEqual(b["external_dependencies"], [])
        self.assertTrue(b["produces_analysis_product"])

    def test_quac_without_engine_is_unavailable(self):
        c = deepcopy(self.cfg)
        c["atmospheric"]["method"] = "envi_quac"
        c["atmospheric"]["envi_taskengine"] = ""
        b = atmos.detect(c)
        self.assertFalse(b["available"])
        self.assertTrue(b["requires_external_software"])
        self.assertEqual(len(b["external_dependencies"]), 2)

    def test_none_backend_delivers_radiance_only(self):
        c = deepcopy(self.cfg)
        c["atmospheric"]["method"] = "none"
        b = atmos.detect(c)
        self.assertTrue(b["available"])
        self.assertFalse(b["produces_analysis_product"])

    def test_invalid_method_rejected_by_config(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "pipeline.yaml"
            p.write_text(CFG_PATH.read_text(encoding="utf-8").replace("method: iarr", "method: bogus"),
                         encoding="utf-8")
            with self.assertRaises(ValueError):
                load_config(p, SENSORS_PATH)

    def test_no_hardcoded_machine_paths_in_code_or_config(self):
        markers = ("C:/Users", "C:\\Users", "F:/", "F:\\", "/home/", "/Users/")
        for p in list((ROOT / "ahsi").glob("*.py")) + list((ROOT / "config").glob("*.yaml")):
            text = p.read_text(encoding="utf-8")
            for marker in markers:
                self.assertNotIn(marker, text, f"{p.name} 含本机绝对路径标记 {marker!r}")


class TestIARRNumerics(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config(CFG_PATH)
        cls.cfg["runtime"] = dict(cls.cfg["runtime"])
        cls.cfg["runtime"]["tiled"] = True
        cls.cfg["runtime"]["block_size"] = 16

    def test_factors_scaling_and_bad_band(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            cube = _synthetic_cube(tmp)
            scene = SimpleNamespace(scene_id="SYNTH")
            rec = relative.compute_factors(scene, self.cfg, cube, tmp / "IARR")
            self.assertAlmostEqual(rec["factors"][0], 10.0, places=5)
            self.assertAlmostEqual(rec["factors"][1], 20.0, places=5)
            self.assertIsNone(rec["factors"][2])
            self.assertEqual(rec["good_bands"], 2)
            self.assertEqual(rec["external_dependencies"], [])
            out = relative.finalize_relative_reflectance(scene, self.cfg, tmp, cube, rec)
            self.assertEqual(out["good_bands"], 2)
            self.assertEqual(out["bad_band_wavelengths_nm"], [700.0])
            self.assertFalse(out["surface_reflectance_available"])
            self.assertTrue(out["product_grade"].startswith("relative"))
            with rasterio.open(out["reflectance_file"]) as d:
                self.assertIn("NOT surface reflectance", d.tags()["product"])
                self.assertIn("relative_reflectance", Path(d.name).name)
                self.assertTrue(np.allclose(d.read(1), 1.0))
                self.assertTrue(np.allclose(d.read(2), 1.0))
                self.assertTrue(np.all(d.read(3) == -9999.0))
                self.assertEqual(d.tags(ns="ENVI")["bbl"], "{1,1,0}")

    def test_deterministic_recompute(self):
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            cube = _synthetic_cube(tmp)
            scene = SimpleNamespace(scene_id="SYNTH")
            a = relative.compute_factors(scene, self.cfg, cube, tmp / "a")
            b = relative.compute_factors(scene, self.cfg, cube, tmp / "b")
            self.assertEqual(a["factors"], b["factors"])
            self.assertEqual(a["band_key"], b["band_key"])

    def test_stale_factors_rejected(self):
        """换了波段表却复用旧因子必须报错，防止结论漂移。"""
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            cube = _synthetic_cube(tmp)
            scene = SimpleNamespace(scene_id="SYNTH")
            rec = relative.compute_factors(scene, self.cfg, cube, tmp / "IARR")
            changed = dict(cube)
            changed["band_records"] = [dict(r) for r in cube["band_records"]]
            changed["band_records"][0]["wavelength_nm"] = 501.0
            with self.assertRaisesRegex(ValueError, "不一致"):
                relative.finalize_relative_reflectance(scene, self.cfg, tmp, changed, rec)


class TestMigrationChecks(unittest.TestCase):
    def test_utm_zone_for_other_study_area(self):
        self.assertEqual(preflight._utm_crs(105.0, 34.0), "EPSG:32648")
        self.assertEqual(preflight._utm_crs(116.6, 41.2), "EPSG:32650")
        self.assertEqual(preflight._utm_crs(118.3, 41.0), "EPSG:32650")

    def test_checklist_is_study_area_agnostic(self):
        cfg = load_config(CFG_PATH)
        items = preflight.migration_checklist(cfg)
        # 清单会如实包含运行时解析出的绝对路径（例如 Linux CI 的 runner 工作目录），
        # 这是预期行为，不是硬编码。硬编码路径由
        # test_no_hardcoded_machine_paths_in_code_or_config 在源码/配置文本层面检查。
        self.assertTrue(all(c.get("item") and c.get("detail") for c in items))
        levels = {c["level"] for c in items}
        self.assertIn("必须人工准备", levels)
        self.assertIn("完全不用改", levels)


if __name__ == "__main__":
    unittest.main(verbosity=2)
