"""地形光照校正模块的几何与产品契约单测（纯合成数据，不读项目成果）。"""
from __future__ import annotations

import math
import unittest

import numpy as np

from ahsi import terrain


class TestSlopeAspect(unittest.TestCase):

    def test_east_facing_uphill_gives_west_aspect(self):
        # z 随 col 增大（东高西低）→ 下坡朝西 = 270°
        z = np.tile(np.arange(10, dtype=float) * 2.0, (10, 1))
        slope, aspect = terrain.slope_aspect(z, cell_m=30.0)
        self.assertAlmostEqual(math.degrees(slope[5, 5]), math.degrees(math.atan(2 / 30)), places=5)
        self.assertAlmostEqual(math.degrees(aspect[5, 5]), 270.0, places=4)

    def test_north_high_south_low_gives_south_aspect(self):
        # z 随 row 减小而增大（北高南低）→ dz_drow<0 → 下坡朝南 = 180°
        z = np.tile(np.arange(10, 0, -1, dtype=float)[:, None] * 2.0, (1, 10))
        _, aspect = terrain.slope_aspect(z, cell_m=30.0)
        self.assertAlmostEqual(math.degrees(aspect[5, 5]), 180.0, places=4)

    def test_flat_terrain_zero_slope(self):
        slope, _ = terrain.slope_aspect(np.full((8, 8), 1500.0), cell_m=30.0)
        self.assertTrue(np.allclose(slope, 0.0))


class TestIllumination(unittest.TestCase):

    def setUp(self):
        self.sun = terrain.SunGeometry(37.623436, 150.452459, source="unit_test")

    def test_flat_equals_cos_zenith(self):
        slope = np.zeros((4, 4))
        aspect = np.zeros((4, 4))
        ci = terrain.illumination_cosine(slope, aspect, self.sun)
        self.assertTrue(np.allclose(ci, math.cos(math.radians(self.sun.zenith_deg))))

    def test_sun_facing_slope_increases_illumination(self):
        # 坡向与太阳方位相同、坡度等于太阳天顶角时，太阳垂直入射
        slope = np.full((3, 3), math.radians(self.sun.zenith_deg))
        aspect = np.full((3, 3), math.radians(self.sun.azimuth_deg))
        ci = terrain.illumination_cosine(slope, aspect, self.sun)
        self.assertTrue(np.allclose(ci, 1.0, atol=1e-12))

    def test_sun_angle_validation(self):
        with self.assertRaises(ValueError):
            terrain.SunGeometry(95.0, 150.0).validate()
        with self.assertRaises(ValueError):
            terrain.SunGeometry(37.0, -5.0).validate()


class TestApplyCorrection(unittest.TestCase):

    def setUp(self):
        self.sun = terrain.SunGeometry(40.0, 150.0, source="unit_test")
        self.h, self.w, self.nb = 12, 12, 3
        self.dem = np.full((self.h, self.w), 1000.0, dtype="float32")
        self.slope, self.aspect = terrain.slope_aspect(self.dem, 30.0)
        self.cos_i = terrain.illumination_cosine(self.slope, self.aspect, self.sun)
        self.wl = np.array([600.0, 1600.0, 2200.0])
        self.rho = np.full((self.nb, self.h, self.w), 0.25, dtype="float32")

    def test_shape_and_flat_terrain_nearly_identity(self):
        r = terrain.apply_correction(self.rho, self.cos_i, self.slope, self.dem,
                                     self.sun, self.wl, min_regression_pixels=10)
        self.assertEqual(r.corrected.shape, self.rho.shape)
        self.assertEqual(r.qa.shape, (self.h, self.w))
        # 平坦地形 cos(i)=cos(theta_s)，SCS+C 分子也是 cos(theta_s) → factor=1
        self.assertTrue(np.allclose(r.corrected, self.rho, atol=1e-6))

    def test_low_illumination_pixels_not_amplified(self):
        cos_i = self.cos_i.copy()
        cos_i[0, 0] = -0.2
        self.rho[:, 0, 0] = 0.3
        r = terrain.apply_correction(self.rho, cos_i, self.slope, self.dem,
                                     self.sun, self.wl, min_regression_pixels=10)
        self.assertTrue(bool(r.qa[0, 0] & terrain.QA_BITS["low_illumination"]))
        self.assertTrue(np.allclose(r.corrected[:, 0, 0], 0.3, atol=1e-6))

    def test_nodata_preserved_and_qa_set(self):
        rho = self.rho.copy()
        rho[:, 1, 1] = -9999.0
        r = terrain.apply_correction(rho, self.cos_i, self.slope, self.dem,
                                     self.sun, self.wl, min_regression_pixels=10)
        self.assertTrue(np.all(r.corrected[:, 1, 1] == -9999.0))
        self.assertTrue(bool(r.qa[1, 1] & terrain.QA_BITS["input_invalid"]))

    def test_band_records_and_caveats_present(self):
        r = terrain.apply_correction(self.rho, self.cos_i, self.slope, self.dem,
                                     self.sun, self.wl, min_regression_pixels=10)
        self.assertEqual(len(r.band_records), self.nb)
        self.assertIn("caveats", r.meta)
        self.assertIn("c_fit_valid_bands", r.meta)
        self.assertFalse(r.meta["sun_geometry"]["pixelwise_ogp_angles_applied"])

    def test_invalid_dem_gradient_keeps_reflectance(self):
        dem = self.dem.copy()
        dem[4, 4] = np.nan
        slope, aspect = terrain.slope_aspect(dem, 30)
        cos_i = terrain.illumination_cosine(slope, aspect, self.sun)
        result = terrain.apply_correction(self.rho, cos_i, slope, dem, self.sun, self.wl)
        self.assertTrue(result.qa[4, 5] & terrain.QA_BITS["dem_invalid"])
        np.testing.assert_array_equal(result.corrected[:, 4, 5], self.rho[:, 4, 5])

    def test_c_interpolation_uses_wavelength_not_band_index(self):
        records = [
            {"wavelength_nm": 600, "c_fit_valid": True, "c_parameter": 1.0},
            {"wavelength_nm": 700, "c_fit_valid": False, "c_parameter": None},
            {"wavelength_nm": 1000, "c_fit_valid": True, "c_parameter": 3.0},
        ]
        terrain._interpolate_c(records)
        self.assertEqual(records[1]["c_parameter"], 1.5)

    def test_degenerate_scs_c_falls_back_to_cosine_and_sets_qa(self):
        slope = np.full_like(self.slope, .4)
        ci = np.full_like(self.cos_i, .5)
        result = terrain.apply_correction(self.rho, ci, slope, self.dem,
                                          self.sun, self.wl, min_regression_pixels=10)
        np.testing.assert_allclose(result.corrected, self.rho * math.cos(math.radians(40)) / .5)
        self.assertTrue(np.all(result.qa & terrain.QA_BITS["c_parameter_fallback"]))
        self.assertEqual(result.meta["c_fit_valid_bands"], 0)

    def test_unknown_method_rejected(self):
        with self.assertRaises(ValueError):
            terrain.apply_correction(self.rho, self.cos_i, self.slope, self.dem,
                                     self.sun, self.wl, method="magic")

    def test_cosine_method_forces_zero_c(self):
        r = terrain.apply_correction(self.rho, self.cos_i, self.slope, self.dem,
                                     self.sun, self.wl, method="cosine",
                                     min_regression_pixels=10)
        self.assertTrue(all(x["c_parameter"] == 0.0 for x in r.band_records))


if __name__ == "__main__":
    unittest.main()
