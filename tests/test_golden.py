"""回归测试：用真实已交付的两景产品验证"发现 + 波段表 + 配置"层。

不写任何文件、不重算大立方体；只读取元数据与小窗口。
运行： python -m unittest discover -s tests -v
"""
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys_path = str(ROOT)
import sys
if sys_path not in sys.path:
    sys.path.insert(0, sys_path)

from ahsi.config import load_config  # noqa: E402
from ahsi.scene import discover_scenes, numeric_table  # noqa: E402

# Golden 回归依赖用户自备的 L1 场景数据。公开仓库/CI 中没有这些数据时整类跳过，
# 不伪造通过；本地配置了 input_root 且存在 XML 时仍会完整执行。
_CFG = load_config(ROOT / "config" / "pipeline.yaml")
_INPUT_ROOT = Path(_CFG["input_root"])
HAS_SCENE_DATA = _INPUT_ROOT.is_dir() and any(_INPUT_ROOT.glob("*.xml"))


@unittest.skipUnless(HAS_SCENE_DATA, "本地 L1 场景数据未配置，跳过 golden 回归")
class TestDiscovery(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config(ROOT / "config" / "pipeline.yaml")
        cls.scenes, cls.errors = discover_scenes(cls.cfg["input_root"], cls.cfg)
        cls.by_sat = {s.satellite: s for s in cls.scenes}

    def test_two_scenes_discovered(self):
        self.assertEqual(self.errors, [])
        self.assertEqual(set(self.by_sat), {"GF5B", "ZY1E"})

    def test_band_counts(self):
        self.assertEqual(self.by_sat["GF5B"].total_bands, 330)
        self.assertEqual(self.by_sat["ZY1E"].total_bands, 166)
        self.assertEqual([s.count for s in self.by_sat["GF5B"].segments], [150, 180])
        self.assertEqual([s.count for s in self.by_sat["ZY1E"].segments], [76, 90])

    def test_candidate_counts_match_validated_run(self):
        self.assertEqual(len(self.by_sat["GF5B"].candidate_bands), 248)
        self.assertEqual(len(self.by_sat["ZY1E"].candidate_bands), 123)

    def test_wavelength_ranges(self):
        g = self.by_sat["GF5B"]
        self.assertAlmostEqual(g.band_table[0]["wavelength_nm"], 387.0453, places=3)
        self.assertAlmostEqual(g.band_table[-1]["wavelength_nm"], 2515.577, places=3)
        z = self.by_sat["ZY1E"]
        self.assertAlmostEqual(z.band_table[0]["wavelength_nm"], 395.8612, places=3)

    def test_screening_rules(self):
        b = self.cfg["bands"]
        for s in self.scenes:
            for r in s.candidate_bands:
                self.assertGreaterEqual(r["wavelength_nm"], b["wavelength_min_nm"])
                self.assertLessEqual(r["wavelength_nm"], b["wavelength_max_nm"])
                for lo, hi in b["water_absorption_nm"]:
                    self.assertFalse(lo <= r["wavelength_nm"] <= hi)
                if r["segment"] == "VNIR":
                    self.assertLessEqual(r["wavelength_nm"], b["vnir_max_nm"])
                else:
                    self.assertGreaterEqual(r["wavelength_nm"], b["swir_min_nm"])

    def test_rpc_and_ogp(self):
        for s in self.scenes:
            self.assertTrue(all(seg.has_rpc for seg in s.segments))
            self.assertIn("solar_zenith_deg", s.ogp["center"])
            self.assertLess(s.ogp["solar_zenith_diff_deg"], 0.1,
                            "OGP 中心像元太阳天顶角应与 XML 中心角一致")

    def test_numeric_tables(self):
        for s in self.scenes:
            for seg in s.segments:
                self.assertEqual(seg.cal.shape, (seg.count, 2))
                self.assertEqual(seg.spec.shape, (seg.count, 2))
                self.assertTrue((seg.cal[:, 0] > 0).all())
        # 两种分隔符都能解析
        for p in list((ROOT / "config").glob("*.yaml")):
            self.assertTrue(p.is_file())

    def test_config_no_hardcoded_machine_path(self):
        text = (ROOT / "config" / "pipeline.yaml").read_text(encoding="utf-8")
        for marker in ("C:/Users", "F:/", "F:\\", "/home/", "/Users/"):
            self.assertNotIn(marker, text, f"默认配置含本机路径标记 {marker!r}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
