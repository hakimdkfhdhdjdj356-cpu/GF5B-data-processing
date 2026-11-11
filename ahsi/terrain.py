"""地形光照辐射校正：纯 numpy / rasterio 实现，零商业软件依赖。

科学定位
--------
本模块只做**一阶太阳直射光照的地形归一化**，不是完整地形辐射校正：
- 已做：DEM 重投影到产品格网 → 坡度/坡向 → 逐像元太阳入射角 cos(i)
  → 逐波段回归求 C 参数 → SCS+C / C-correction / cosine 校正 → QA 位。
- 未做：天空漫射/各向异性、邻近地形反射、投射阴影光线追踪、BRDF、
  大气参数随高程/坡面的再参数化、逐像元 OGP 观测角。
产品等级因此是「模型估计反射率 + 地形光照一阶校正」，**不是真值、未独立验证**。

方法
----
Soenen et al. (2005) SCS+C：
    rho_c = rho_obs * (cos(alpha) * cos(theta_s) + C_lambda) / (cos(i) + C_lambda)
其中 alpha 为坡度、theta_s 为太阳天顶角、C_lambda = b/m 来自逐波段
    rho_obs = m * cos(i) + b
的线性回归。C-correction 把分子换成 (cos(theta_s) + C)，cosine 则是 C=0。

重要边界：C_lambda 逐波段不同，因此校正因子逐波段不同，**会改变光谱向量
方向与谱形**。这是该方法的已知代价，必须随成果一起说明；用于 SAM / 光谱库
匹配前需重新评估，不得默认"校正后谱形不变"。
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject

QA_BITS = {
    "low_illumination": 1,       # cos(i) <= 阈值，自阴影/近切入射，不校正
    "corrected_above_1": 2,
    "corrected_below_0": 4,
    "c_parameter_fallback": 8,   # 像元至少一波段使用插值 C 或余弦回退
    "dem_invalid": 16,           # 高程或梯度无效
    "input_invalid": 32,
    "outside_coverage": 64,
}


@dataclass(frozen=True)
class SunGeometry:
    zenith_deg: float
    azimuth_deg: float            # 从北顺时针
    source: str = "xml_center"

    def validate(self) -> None:
        if not (0.0 <= float(self.zenith_deg) < 90.0):
            raise ValueError(f"太阳天顶角超出 [0,90): {self.zenith_deg}")
        if not (0.0 <= float(self.azimuth_deg) <= 360.0):
            raise ValueError(f"太阳方位角超出 [0,360]: {self.azimuth_deg}")


@dataclass
class TerrainResult:
    corrected: np.ndarray         # (bands, rows, cols) float32
    qa: np.ndarray                # (rows, cols) uint8
    slope_deg: np.ndarray
    aspect_deg: np.ndarray
    cos_i: np.ndarray
    dem_on_grid: np.ndarray
    band_records: list[dict]
    meta: dict


def slope_aspect(dem: np.ndarray, cell_m: float) -> tuple[np.ndarray, np.ndarray]:
    """由投影格网高程求坡度(rad)与坡向(rad, 从北顺时针, [0,2pi))。

    数组约定：row 向南增大、col 向东增大。推导：
        gx = dz/dx(东) = dz_dcol / cell ;  gy = dz/dy(北) = -dz_drow / cell
        下坡水平向量 = (-gx, -gy)，方位 = atan2(东分量, 北分量)
                      = atan2(-dz_dcol, dz_drow)
    纯东向升高(dz_dcol>0) → 下坡朝西 → 270°，由单测固定符号。
    """
    z = np.asarray(dem, dtype="float64")
    dz_drow, dz_dcol = np.gradient(z, float(cell_m), float(cell_m))
    slope = np.arctan(np.hypot(dz_dcol, dz_drow))
    aspect = np.mod(np.arctan2(-dz_dcol, dz_drow), 2.0 * np.pi)
    return slope, aspect


def illumination_cosine(slope: np.ndarray, aspect: np.ndarray,
                        sun: SunGeometry) -> np.ndarray:
    """逐像元 cos(i)：太阳直射相对坡面法线的入射余弦。"""
    ts = math.radians(float(sun.zenith_deg))
    ps = math.radians(float(sun.azimuth_deg))
    return (np.cos(slope) * math.cos(ts)
            + np.sin(slope) * math.sin(ts) * np.cos(ps - aspect))


def reproject_dem(dem_path: str | Path, crs, transform,
                  width: int, height: int,
                  resampling: Resampling = Resampling.bilinear) -> np.ndarray:
    """把 RPC 正射用的 DEM 重采样到反射率产品格网；返回 float32, 无效为 nan。"""
    out = np.full((height, width), np.nan, dtype="float32")
    with rasterio.open(dem_path) as src:
        if src.count < 1:
            raise ValueError("DEM 无波段")
        reproject(rasterio.band(src, 1), out,
                  src_crs=src.crs, src_transform=src.transform,
                  src_nodata=src.nodata,
                  dst_crs=crs, dst_transform=transform,
                  dst_nodata=np.nan, resampling=resampling,
                  num_threads=2, warp_mem_limit=256)
    return out


def _fit_c(rho: np.ndarray, cos_i: np.ndarray, mask: np.ndarray,
           min_samples: int) -> dict:
    """逐波段回归 rho = m*cos(i) + b；返回 C=b/m 与拟合质量。"""
    ok = mask & np.isfinite(rho) & np.isfinite(cos_i) & (rho > 0) & (rho < 1.5)
    n = int(ok.sum())
    rec = {"n_regression_pixels": n, "c_fit_valid": False,
           "slope_m": None, "intercept_b": None, "c_parameter": None,
           "r_squared": None}
    if n < min_samples:
        rec["reason"] = f"regression pixels {n} < {min_samples}"
        return rec
    x = cos_i[ok].astype("float64")
    y = rho[ok].astype("float64")
    x_mean, y_mean = float(x.mean()), float(y.mean())
    dx = x - x_mean
    ss_x = float(np.dot(dx, dx))
    if ss_x <= np.finfo("float64").eps * n:
        rec["reason"] = "insufficient incidence-angle variation"
        return rec
    m = float(np.dot(dx, y - y_mean) / ss_x)
    b = y_mean - m * x_mean
    yhat = m * x + b
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else None
    c = float(b / m) if m > 0 else None
    rec.update(slope_m=float(m), intercept_b=float(b),
               c_parameter=c, r_squared=r2)
    if c is None or not math.isfinite(c) or c <= 0:
        rec["reason"] = "non-positive or non-finite C; refuse to fabricate"
        return rec
    rec["c_fit_valid"] = True
    return rec


def _interpolate_c(records: list[dict]) -> None:
    """对回归失败的波段，用有效 C 的波长线性插值；无有效值则回退 C=0(cosine)。"""
    good = [i for i, r in enumerate(records) if r["c_fit_valid"]]
    if not good:
        for r in records:
            r["c_parameter"] = 0.0
            r["c_source"] = "fallback_cosine_zero_no_valid_regression"
        return
    wave = np.asarray([r["wavelength_nm"] for r in records], dtype=float)
    vals = [records[i]["c_parameter"] for i in good]
    interp = np.interp(wave, wave[good], vals)
    for i, r in enumerate(records):
        if r["c_fit_valid"]:
            r["c_source"] = "per_band_regression"
        else:
            r["c_parameter"] = float(interp[i])
            r["c_source"] = "interpolated_from_valid_bands"


def apply_correction(reflectance: np.ndarray, cos_i: np.ndarray,
                     slope: np.ndarray, dem: np.ndarray,
                     sun: SunGeometry, wavelength_nm: np.ndarray,
                     method: str = "scs_c",
                     min_cos_i: float = 0.10,
                     min_regression_pixels: int = 1000,
                     nodata: float = -9999.0,
                     coverage: np.ndarray | None = None,
                     aspect: np.ndarray | None = None) -> TerrainResult:
    """对 (bands, rows, cols) 反射率立方体做地形光照校正。

    low-illumination、DEM 无效、输入无效像元一律**不校正**（factor=1）并打 QA 位，
    绝不通过除以接近 0 的 cos(i) 放大噪声；输出不裁剪，越界只统计并标记。
    """
    sun.validate()
    if method not in ("scs_c", "c_correction", "cosine"):
        raise ValueError(f"未知地形校正方法 {method}")
    rho = np.asarray(reflectance, dtype="float32")
    if rho.ndim != 3:
        raise ValueError("reflectance 必须是 (bands, rows, cols)")
    nb, h, w = rho.shape
    if cos_i.shape != (h, w) or slope.shape != (h, w) or dem.shape != (h, w):
        raise ValueError("cos_i / slope / dem 与反射率格网不一致")
    wave = np.asarray(wavelength_nm, dtype=float)
    if wave.shape != (nb,) or not np.isfinite(wave).all() or not np.all(np.diff(wave) > 0):
        raise ValueError("波长必须有限、严格递增且与波段数一致")
    if not 0 <= min_cos_i < 1 or min_regression_pixels < 2:
        raise ValueError("需要 0 <= min_cos_i < 1 且至少两个回归像元")

    ts = math.radians(float(sun.zenith_deg))
    cos_ts = math.cos(ts)
    elevation_ok = np.isfinite(dem)
    dem_ok = elevation_ok & np.isfinite(slope) & np.isfinite(cos_i)
    input_ok = np.isfinite(rho) & (rho != nodata)
    illum_ok = dem_ok & (cos_i > min_cos_i)
    cov = np.ones((h, w), bool) if coverage is None else (np.asarray(coverage) > 0)
    if cov.shape != (h, w):
        raise ValueError("coverage 与反射率格网尺寸不一致")
    regression_mask = illum_ok & cov

    qa = np.zeros((h, w), dtype="uint8")
    qa[~dem_ok] |= QA_BITS["dem_invalid"]
    qa[~illum_ok & dem_ok] |= QA_BITS["low_illumination"]
    qa[~np.all(input_ok, axis=0)] |= QA_BITS["input_invalid"]
    qa[~cov] |= QA_BITS["outside_coverage"]

    records: list[dict] = []
    for b in range(nb):
        if method == "cosine":
            rec = {"c_fit_valid": False, "n_regression_pixels": 0,
                   "c_parameter": 0.0, "c_source": "cosine_method_forced_zero"}
        else:
            rec = _fit_c(rho[b], cos_i, regression_mask & input_ok[b], min_regression_pixels)
        rec.update(band=b + 1, wavelength_nm=float(wave[b]))
        records.append(rec)
    if method != "cosine":
        _interpolate_c(records)

    cos_alpha = np.cos(slope)
    numerator_base = (cos_alpha * cos_ts if method == "scs_c"
                      else np.full((h, w), cos_ts, dtype="float64"))
    out = rho.copy()
    for b, rec in enumerate(records):
        c = float(rec["c_parameter"])
        den = cos_i + c
        factor = np.ones((h, w), dtype="float32")
        usable = illum_ok & cov & input_ok[b] & np.isfinite(den) & (den > 1e-6)
        fallback_cosine = rec["c_source"] == "fallback_cosine_zero_no_valid_regression"
        numerator = cos_ts + c if fallback_cosine else numerator_base[usable] + c
        factor[usable] = (numerator / den[usable]).astype("float32")
        if method != "cosine" and not rec["c_fit_valid"]:
            qa[usable] |= QA_BITS["c_parameter_fallback"]
        out[b] = np.where(input_ok[b], rho[b] * factor, nodata).astype("float32")

    finite_out = np.isfinite(out) & (out != nodata)
    bad_hi = np.any(finite_out & (out > 1.0), axis=0)
    bad_lo = np.any(finite_out & (out < 0.0), axis=0)
    qa[bad_hi] |= QA_BITS["corrected_above_1"]
    qa[bad_lo] |= QA_BITS["corrected_below_0"]

    valid = np.all(input_ok, axis=0) & illum_ok & cov
    before = rho[:, valid]
    after = out[:, valid]
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(before > 0, (after - before) / before, np.nan)
    finite_rel = rel[np.isfinite(rel)]

    meta = {
        "method": method,
        "sun_geometry": {"zenith_deg": sun.zenith_deg,
                         "azimuth_deg": sun.azimuth_deg,
                         "source": sun.source,
                         "pixelwise_ogp_angles_applied": False},
        "min_cos_i": min_cos_i,
        "min_regression_pixels": min_regression_pixels,
        "bands": nb, "grid": [h, w],
        "dem_valid_fraction": float(elevation_ok.mean()),
        "dem_gradient_valid_fraction": float(dem_ok.mean()),
        "low_illumination_fraction": float((~illum_ok & dem_ok).mean()),
        "regression_pixel_count": int(regression_mask.sum()),
        "c_fit_valid_bands": int(sum(r["c_fit_valid"] for r in records)),
        "c_interpolated_bands": [r["band"] for r in records
                                 if r.get("c_source") == "interpolated_from_valid_bands"],
        "c_fallback_reasons": {str(r["band"]): r["reason"] for r in records
                               if not r["c_fit_valid"] and "reason" in r},
        "slope_deg": {"min": float(np.degrees(slope[dem_ok]).min()) if dem_ok.any() else None,
                      "median": float(np.median(np.degrees(slope[dem_ok]))) if dem_ok.any() else None,
                      "max": float(np.degrees(slope[dem_ok]).max()) if dem_ok.any() else None},
        "cos_i": {"min": float(cos_i[dem_ok].min()) if dem_ok.any() else None,
                  "median": float(np.median(cos_i[dem_ok])) if dem_ok.any() else None,
                  "max": float(cos_i[dem_ok].max()) if dem_ok.any() else None},
        "qa_fraction": {k: float((qa & bit > 0).mean()) for k, bit in QA_BITS.items()},
        "relative_change_median": float(np.median(finite_rel)) if finite_rel.size else None,
        "relative_change_p05": float(np.percentile(finite_rel, 5)) if finite_rel.size else None,
        "relative_change_p95": float(np.percentile(finite_rel, 95)) if finite_rel.size else None,
        "caveats": [
            "Only first-order direct-solar illumination normalization; diffuse sky, terrain "
            "reflected radiance, cast shadows, BRDF and adjacency are NOT modeled.",
            "C_lambda differs per band, so the correction CHANGES spectral vector direction; "
            "re-evaluate before SAM / library matching.",
            "Slope and aspect depend on DEM source, resolution and height reference.",
            "Single scene-center sun geometry from XML; pixelwise OGP angles are NOT applied.",
            "Low-illumination and invalid pixels are left uncorrected with QA bits; output is NOT clipped.",
        ],
    }
    aspect_deg = (np.degrees(aspect) if aspect is not None
                  else np.full((h, w), np.nan, dtype="float32"))
    return TerrainResult(corrected=out, qa=qa, slope_deg=np.degrees(slope),
                         aspect_deg=aspect_deg, cos_i=cos_i,
                         dem_on_grid=dem, band_records=records, meta=meta)
