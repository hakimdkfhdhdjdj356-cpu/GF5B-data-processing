"""6SV LUT 后端的纯数值 / 契约回归单测（不调用 6S 可执行）。

覆盖：2.5 nm SRF 网格约束、反演公式与 6S 自报值对账、后端注册与降级、
配置校验、LUT 漂移保护。真实 6S 运行由正式 run 的 self_check 覆盖。
"""
import sys
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ahsi import atmos, sixs_lut  # noqa: E402
from ahsi.config import load_config  # noqa: E402

CFG_PATH = ROOT / "config" / "pipeline.yaml"
SENSORS_PATH = ROOT / "config" / "sensors.yaml"

# 本机 6S 实跑（2026-11-11）得到的一组系数与自报反射率，用作数值锚点
REAL_XA, REAL_XB, REAL_XC = 0.00287, 0.02232, 0.05935
REAL_L, REAL_SIXS_REPORTED = 100.0, 0.260880


class TestInversion(unittest.TestCase):
    def test_formula_matches_sixs_reported_within_print_precision(self):
        r = float(sixs_lut.surface_reflectance_from_radiance(
            REAL_L, REAL_XA, REAL_XB, REAL_XC))
        self.assertAlmostEqual(r, 0.260587, places=6)
        self.assertLess(abs(r - REAL_SIXS_REPORTED), 1e-3)

    def test_coefficients_do_not_depend_on_radiance(self):
        # 6S 实跑 L=50/100/200 得到同一组系数；反演式对 L 单调且系数固定
        vals = [float(sixs_lut.surface_reflectance_from_radiance(
            x, REAL_XA, REAL_XB, REAL_XC)) for x in (50.0, 100.0, 200.0)]
        self.assertTrue(all(a < b for a, b in zip(vals, vals[1:])))
        self.assertAlmostEqual(vals[1], 0.260587, places=6)

    def test_nonfinite_denominator_returns_nan(self):
        r = sixs_lut.surface_reflectance_from_radiance(
            np.array([100.0, np.nan, -9999.0]), REAL_XA, REAL_XB, REAL_XC)
        self.assertTrue(np.isfinite(r[0]))
        self.assertTrue(np.isnan(r[1]))

    def test_array_input(self):
        r = sixs_lut.surface_reflectance_from_radiance(
            np.array([50.0, 100.0]), REAL_XA, REAL_XB, REAL_XC)
        self.assertEqual(r.shape, (2,))

    def test_nonfinite_6s_coefficients_must_be_rejected(self):
        # 正式 run 实测：2018.996 nm + 2xFWHM 时 6S 返回全 NaN 系数
        self.assertTrue(sixs_lut._coefficients_finite(REAL_XA, REAL_XB, REAL_XC))
        self.assertFalse(sixs_lut._coefficients_finite(float("nan"), REAL_XB, REAL_XC))
        self.assertFalse(sixs_lut._coefficients_finite(REAL_XA, float("inf"), REAL_XC))
        self.assertFalse(sixs_lut._coefficients_finite(REAL_XA, REAL_XB, None))


class TestSRF(unittest.TestCase):
    def test_2p5nm_grid_and_peak(self):
        start, end, srf = sixs_lut.gaussian_srf(0.65, 0.010)
        self.assertEqual(len(srf), round((end - start) / 0.0025) + 1)
        self.assertAlmostEqual((end - start) / (len(srf) - 1), 0.0025, places=9)
        self.assertAlmostEqual(float(srf.max()), 1.0, places=6)
        self.assertGreater(start, 0.25)   # 6S 可接受的波长范围
        self.assertLess(end, 4.0)

    def test_off_grid_center_is_sampled_on_fixed_6sv_axis(self):
        start, stop, response = sixs_lut.gaussian_srf(.4212525, .00433)
        grid = start + .0025 * np.arange(len(response))
        ticks = (grid - .25) / .0025
        self.assertTrue(np.allclose(ticks, np.rint(ticks), rtol=0, atol=1e-10))
        self.assertAlmostEqual(grid[-1], stop)
        self.assertAlmostEqual(grid[np.argmax(response)], .4225)

    def test_widening_fwhm_increases_points(self):
        a = sixs_lut.gaussian_srf(2.2, 0.015, 1.0)
        b = sixs_lut.gaussian_srf(2.2, 0.015, 2.0)
        self.assertGreater(len(b[2]), len(a[2]))
        self.assertEqual(len(b[2]), round((b[1] - b[0]) / 0.0025) + 1)

    def test_invalid_fwhm_rejected(self):
        with self.assertRaises(ValueError):
            sixs_lut.gaussian_srf(0.65, 0.0)

    def test_sensitivity_indices_cover_ends(self):
        idx = sixs_lut._sensitivity_band_indices(254, 6)
        self.assertIn(1, idx)
        self.assertIn(254, idx)
        self.assertEqual(len(idx), 6)


class TestBackendContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config(CFG_PATH)

    def test_geometry_uses_scene_xml_in_a_multi_scene_root(self):
        fields = {"SolarZenith": 40, "SolarAzimuth": 150, "SatelliteZenith": 1,
                  "SatelliteAzimuth": 190, "SunEarthDistance": 1, "CenterLatitude": 34,
                  "CenterLongitude": 105, "StartTime": "2024-03-21 11:57:25"}
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            (root / "a.xml").write_text("<Product><SolarZenith>7</SolarZenith></Product>")
            folder = root / "scene_b"
            folder.mkdir()
            selected = folder / "b.xml"
            selected.write_text("<Product>" + "".join(f"<{k}>{v}</{k}>" for k, v in fields.items()) +
                                "</Product>")
            geom = sixs_lut.read_geometry(SimpleNamespace(files={"metadata": selected}))
            self.assertEqual(geom["solar_zenith_deg"], 40)
            self.assertEqual(geom["source_xml"], str(selected))

    def test_method_registered(self):
        self.assertIn("sixs_lut", atmos.METHOD_INFO)
        info = atmos.METHOD_INFO["sixs_lut"]
        self.assertTrue(info["produces_analysis_product"])
        self.assertFalse(info["requires_external_software"])
        self.assertTrue(info["requires_external_binary"])
        self.assertTrue(info["product_grade"].startswith("sixs_"))

    def test_detect_unavailable_when_binary_missing(self):
        c = deepcopy(self.cfg)
        c["atmospheric"]["method"] = "sixs_lut"
        c["atmospheric"]["sixs"] = {
            "executable": "C:/definitely/not/here/sixs.exe",
            "scenarios": [{"id": "S", "aot550": 0.1,
                           "water_g_cm2": 1.0, "ozone_atm_cm": 0.3, "nominal": True}]}
        b = atmos.detect(c)
        self.assertFalse(b["available"])
        self.assertIn("6SV1.1", b["external_dependencies"][0])

    def test_backend_degrades_without_fabricating_reflectance(self):
        c = deepcopy(self.cfg)
        c["atmospheric"]["method"] = "sixs_lut"
        c["atmospheric"]["sixs"] = {"executable": "C:/definitely/not/here/sixs.exe"}
        b = atmos.detect(c)
        rec = atmos.run_backend(c, SimpleNamespace(scene_id="SYNTH"), {}, Path("."), b)
        self.assertFalse(rec["successful"])
        self.assertEqual(rec["reason"], "sixs_executable_or_py6s_unavailable")

    def test_config_rejects_missing_nominal_scenario(self):
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "pipeline.yaml"
            txt = CFG_PATH.read_text(encoding="utf-8").replace(
                "method: iarr", "method: sixs_lut")
            txt += ("\natmospheric_sixs_placeholder: unused\n")
            p.write_text(txt, encoding="utf-8")
            with self.assertRaises(ValueError):
                load_config(p, SENSORS_PATH)

    def test_lut_drift_rejected_before_any_write(self):
        with tempfile.TemporaryDirectory() as td:
            lut = Path(td) / "sixs_lut.csv"
            lut.write_text("output_band,wavelength_nm\n1,650\n", encoding="utf-8")
            rec = {"successful": True, "backend_method": "sixs_lut",
                   "lut_file": str(lut), "lut_sha256": "0" * 64,
                   "band_key": [], "nominal_scenario": "S"}
            with self.assertRaisesRegex(ValueError, "内容已变化"):
                sixs_lut.finalize_sixs_reflectance(
                    SimpleNamespace(scene_id="SYNTH"), self.cfg,
                    Path(td), {"band_records": []}, rec)


if __name__ == "__main__":
    unittest.main(verbosity=2)
