"""6SV (6S V1.1) LUT 物理大气校正后端。

用本机开源 6S 可执行 + Py6S 生成 (xa, xb, xc) 查找表，把 at-sensor 辐亮度
反演为模型估计的地表反射率：

    y   = xa * L - xb                      L: W m-2 sr-1 um-1
    rho = y / (1 + xc * y)                 rho: 无量纲地表反射率

严格边界（产物 tag、QC、交付文档必须同时声明）：
- 产物是 **模型估计的地表反射率**，不是实测真值；未经独立精度验证。
- 本景缺少同步水汽/气溶胶/臭氧观测，故采用显式假设情景（scenarios）并
  输出敏感性带；**不得**把任一情景写成"实测大气参数"。
- 中心波长 + FWHM 只是近似 SRF；默认按高斯构造 2.5 nm 等间隔 SRF，
  并额外跑 2xFWHM 敏感性。任何回退到单色近似都必须落盘记录。
- 6S 输出的 xa/xb/xc 只有 5 位小数，反演相对误差约 1e-3 量级；
  每个 (波段, 情景) 都与 6S 自报反射率做公式对账并落盘。
- 未做地形光照校正、邻近效应校正、BRDF 校正、云/雪掩膜。
- 不读取、不修改任何许可文件；6S 为 conda-forge 开源本机二进制，无商业许可依赖。

分两段执行：build_lut 生成 LUT 与元数据（可被断点续跑复用），
finalize_sixs_reflectance 做第二遍写出，并用 LUT sha256 防止"换了配置复用旧表"。
"""
from __future__ import annotations

import hashlib
import io
import json
import logging
import math
import os
import xml.etree.ElementTree as ET
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import rasterio
from PIL import Image

from .core import save_json, write_csv

log = logging.getLogger(__name__)

SRF_STEP_UM = 0.0025          # 6S 强制要求 2.5 nm 等间隔
SRF_HALF_WIDTH_FWHM = 4.0    # 高斯截断于中心 ±4 FWHM
REFERENCE_RADIANCE = 100.0   # 仅用于触发 6S atmos-corr 输出，不影响 xa/xb/xc
SIXS_SELFTEST_EXPECTED = 619.158
PRODUCT_GRADE = "sixs_model_estimated_surface_reflectance_NOT_independently_validated"
LUT_SCHEMA_VERSION = "1.1"


# --------------------------------------------------------------------------
# 纯函数（可单测，不依赖 6S）
# --------------------------------------------------------------------------
def surface_reflectance_from_radiance(radiance, xa: float, xb: float, xc: float):
    """6S 大气校正反演式。radiance 单位必须是 W m-2 sr-1 um-1。"""
    L = np.asarray(radiance, dtype="float64")
    y = float(xa) * L - float(xb)
    denom = 1.0 + float(xc) * y
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = np.where(np.isfinite(denom) & (denom > 1e-12), y / denom, np.nan)
    return rho


def gaussian_srf(center_um: float, fwhm_um: float, scale: float = 1.0):
    """按 6S 要求构造 2.5 nm 等间隔高斯 SRF，返回 (start, end, values)。

    6SV 将输入滤波系数放入从 0.25 μm 起算的固定网格。
    """
    if not (np.isfinite(center_um) and .25 <= center_um <= 4.0
            and np.isfinite(fwhm_um) and fwhm_um > 0 and np.isfinite(scale) and scale > 0):
        raise ValueError("中心波长、FWHM 或响应宽度比例超出有效范围")
    sigma_fwhm = fwhm_um * float(scale)
    half = SRF_HALF_WIDTH_FWHM * sigma_fwhm
    lo = math.ceil((max(.25, center_um - half) - .25) / SRF_STEP_UM - 1e-9)
    hi = math.floor((min(4.0, center_um + half) - .25) / SRF_STEP_UM + 1e-9)
    if not 3 <= hi - lo + 1 <= 1501:
        raise ValueError("SRF 点数超出有效范围")
    grid = .25 + SRF_STEP_UM * np.arange(lo, hi + 1)
    values = np.exp(-4.0 * math.log(2.0) * ((grid - center_um) / sigma_fwhm) ** 2)
    return float(grid[0]), float(grid[-1]), values


def _coefficients_finite(xa, xb, xc) -> bool:
    """6S 在部分宽 SRF / 强吸收组合下输出 NaN；任何非有限/非数值系数都必须被拦截。"""
    try:
        return bool(all(np.isfinite(float(v)) for v in (xa, xb, xc)))
    except (TypeError, ValueError):
        return False


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------
# 输入参数：全部来自真实文件，配置里不写死数值
# --------------------------------------------------------------------------
def read_geometry(scene) -> dict:
    """读取当前场景 XML 的观测几何。"""
    path = Path(scene.files["metadata"])
    tree = ET.parse(path)
    fields = {element.tag.rsplit("}", 1)[-1]: element.text for element in tree.iter()}
    get = fields.get
    def num(tag):
        v = get(tag)
        if v is None:
            raise ValueError(f"XML 缺少字段 {tag}: {path}")
        return float(v)
    start = (get("StartTime") or "").strip()
    month = day = None
    if len(start) >= 10:
        try:
            month, day = int(start[5:7]), int(start[8:10])
        except ValueError:
            pass
    if month is None:
        raise ValueError(f"无法从 StartTime 解析月/日: {start!r}")
    return {
        "source_xml": str(path),
        "solar_zenith_deg": num("SolarZenith"),
        "solar_azimuth_deg": num("SolarAzimuth"),
        "view_zenith_deg": num("SatelliteZenith"),
        "view_azimuth_deg": num("SatelliteAzimuth"),
        "sun_earth_distance_au": num("SunEarthDistance"),
        "start_time": start,
        "month": month, "day": day,
        "center_lat_deg": num("CenterLatitude"),
        "center_lon_deg": num("CenterLongitude"),
        "note": ("SunEarthDistance 由 XML 记录但 6S/Py6S 用 month/day 自行计算日地距离；"
                 "本字段仅作交叉核对，不重复施加修正"),
    }


def read_target_altitude_km(cfg: dict) -> dict:
    """从 RPC DEM 取景内中位高程作为 6S 目标海拔（标量，非逐像元地形校正）。"""
    from .config import aux_path
    p = aux_path(cfg, "rpc_dem")
    with rasterio.open(p) as d:
        a = d.read(1).astype("float64")
        nd = d.nodata
        v = a[np.isfinite(a)] if nd is None else a[np.isfinite(a) & (a != nd)]
    if v.size == 0:
        raise ValueError(f"DEM 无有效高程: {p}")
    return {"source_dem": str(p), "median_m": float(np.median(v)),
            "min_m": float(v.min()), "max_m": float(v.max()),
            "target_altitude_km": float(np.median(v)) / 1000.0,
            "note": "整幅 DEM 中位数；未做逐像元高程/地形光照/邻近效应校正"}


def _scenarios(cfg: dict) -> list[dict]:
    raw = cfg["atmospheric"]["sixs"]["scenarios"]
    nominal = [s for s in raw if s.get("nominal")]
    if len(nominal) != 1:
        raise ValueError("sixs.scenarios 必须且只能有一个 nominal=true 的情景")
    return raw


def self_check(exe: Path) -> dict:
    """6S 官方自检；不单信返回码，另做数值严格比对。"""
    from Py6S import SixS
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = SixS.test(str(exe))
    text = buf.getvalue()
    actual = None
    for line in text.splitlines():
        if line.startswith("Actual result:"):
            try:
                actual = float(line.split(":")[1])
            except ValueError:
                pass
    version = None
    for line in text.splitlines():
        if line.startswith("6sV version:"):
            version = line.split(":", 1)[1].strip()
    ok = (rc == 0 and actual is not None
          and abs(actual - SIXS_SELFTEST_EXPECTED) < 1e-3)
    return {"passed": bool(ok), "return_code": int(rc),
            "expected": SIXS_SELFTEST_EXPECTED, "actual": actual,
            "sixs_version": version, "stdout": text.strip()}


# --------------------------------------------------------------------------
# 6S 单次运行
# --------------------------------------------------------------------------
def _run_one(exe: Path, geom: dict, altitude_km: float, scenario: dict,
             center_nm: float, fwhm_nm: float, srf_scale: float = 1.0) -> dict:
    from Py6S import (SixS, Geometry, Wavelength, AtmosProfile, AeroProfile,
                      AtmosCorr, Altitudes)

    s = SixS(str(exe))
    g = Geometry.User()
    g.solar_z = geom["solar_zenith_deg"]
    g.solar_a = geom["solar_azimuth_deg"]
    g.view_z = geom["view_zenith_deg"]
    g.view_a = geom["view_azimuth_deg"]
    g.month = geom["month"]
    g.day = geom["day"]
    s.geometry = g

    s.atmos_profile = AtmosProfile.UserWaterAndOzone(
        float(scenario["water_g_cm2"]), float(scenario["ozone_atm_cm"]))
    s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Continental)
    s.aot550 = float(scenario["aot550"])

    center_um, fwhm_um = center_nm / 1000.0, fwhm_nm / 1000.0
    srf_mode, n_points, start_um, end_um = "monochromatic_fallback", 0, None, None
    try:
        start_um, end_um, srf = gaussian_srf(center_um, fwhm_um, srf_scale)
        s.wavelength = Wavelength(start_um, end_um, srf.tolist())
        srf_mode, n_points = "gaussian_2p5nm", len(srf)
    except Exception as exc:      # 6S 波长范围/点数不接受时如实回退并留痕
        log.warning("波段 %.3f nm 自定义 SRF 失败，回退单色近似: %s", center_nm, exc)
        s.wavelength = Wavelength(center_um)

    alt = Altitudes()
    alt.set_target_custom_altitude(float(altitude_km))
    alt.set_sensor_satellite_level()
    s.altitudes = alt
    s.atmos_corr = AtmosCorr.AtmosCorrLambertianFromRadiance(REFERENCE_RADIANCE)
    s.run()
    o = s.outputs
    xa, xb, xc = float(o.coef_xa), float(o.coef_xb), float(o.coef_xc)
    finite_check = "ok"
    # 6S 在部分宽 SRF / 强吸收波段组合下会输出 NaN 系数（实测 2018.996 nm, 2xFWHM）。
    # 绝不静默把 NaN 写进 LUT：先回退单色重算并留痕，仍非有限则直接失败。
    if not _coefficients_finite(xa, xb, xc):
        if srf_mode != "gaussian_2p5nm":
            raise RuntimeError(
                f"6S 在 {center_nm:.3f} nm 返回非有限系数且已无 SRF 可回退；"
                "拒绝把 NaN 写入 LUT")
        log.warning("波段 %.3f nm 自定义 SRF 下 6S 系数非有限，回退单色重算", center_nm)
        s.wavelength = Wavelength(center_um)
        srf_mode, n_points = "monochromatic_fallback_nonfinite", 0
        start_um = end_um = None
        s.run()
        o = s.outputs
        xa, xb, xc = float(o.coef_xa), float(o.coef_xb), float(o.coef_xc)
        finite_check = "gaussian_srf_nonfinite_fell_back_to_monochromatic"
        if not _coefficients_finite(xa, xb, xc):
            raise RuntimeError(
                f"6S 在 {center_nm:.3f} nm 单色回退后仍返回非有限系数；拒绝伪造")
    try:
        reported = float(o.atmos_corrected_reflectance_lambertian)
    except (TypeError, ValueError):
        reported = None
    if reported is not None and not math.isfinite(reported):
        reported = None
    formula = float(surface_reflectance_from_radiance(
        REFERENCE_RADIANCE, xa, xb, xc))
    if not math.isfinite(formula):
        raise RuntimeError(f"6S 系数在 {center_nm:.3f} nm 未通过有限反演检查")
    return {
        "xa": xa, "xb": xb, "xc": xc,
        "sixs_reported_reflectance_at_L0": reported,
        "formula_reflectance_at_L0": formula,
        "formula_abs_diff": (abs(reported - formula) if reported is not None else None),
        "srf_mode": srf_mode, "srf_points": n_points,
        "finite_check": finite_check,
        "srf_start_um": start_um, "srf_end_um": end_um,
        "srf_fwhm_nm": fwhm_nm * srf_scale,
        "direct_solar_irradiance": float(o.direct_solar_irradiance),
        "diffuse_solar_irradiance": float(o.diffuse_solar_irradiance),
    }


# --------------------------------------------------------------------------
# 第一段：LUT 生成（atmos.run_backend 调用）
# --------------------------------------------------------------------------
def build_lut(cfg: dict, scene, cube_result: dict, work_dir: Path) -> dict:
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    sconf = cfg["atmospheric"]["sixs"]
    exe = Path(sconf["executable"])
    if not exe.is_file():
        return {"successful": False, "backend_method": "sixs_lut",
                "reason": "sixs_executable_not_found", "executable": str(exe)}
    os.environ["PATH"] = str(exe.parent) + os.pathsep + os.environ.get("PATH", "")
    try:
        import Py6S  # noqa: F401
    except ImportError as exc:
        return {"successful": False, "backend_method": "sixs_lut",
                "reason": "py6s_not_importable", "detail": str(exc)}

    check = self_check(exe)
    if not check["passed"]:
        return {"successful": False, "backend_method": "sixs_lut",
                "reason": "sixs_selftest_failed", "self_test": check}

    geom = read_geometry(scene)
    alt = read_target_altitude_km(cfg)
    scenarios = _scenarios(cfg)
    rows = [r for r in cube_result["band_records"] if r["analysis_candidate"]]
    if not rows:
        raise ValueError("没有候选波段，无法生成 6S LUT")
    nominal = next(s for s in scenarios if s.get("nominal"))

    lut_path = work_dir / "sixs_lut.csv"
    lut_rows, srf_rows = [], []
    total = len(rows) * len(scenarios)
    done = 0
    for k, row in enumerate(rows, 1):
        wl, fw = float(row["wavelength_nm"]), float(row["fwhm_nm"])
        for sc in scenarios:
            rec = _run_one(exe, geom, alt["target_altitude_km"], sc, wl, fw)
            lut_rows.append({
                "output_band": k, "wavelength_nm": wl, "fwhm_nm": fw,
                "source_cube_band": int(row["smoke_output_band"]),
                "scenario": sc["id"], "nominal": bool(sc.get("nominal", False)),
                "aot550": sc["aot550"], "water_g_cm2": sc["water_g_cm2"],
                "ozone_atm_cm": sc["ozone_atm_cm"],
                **rec})
            done += 1
            if done % 100 == 0:
                log.info("6S LUT 进度 %d/%d", done, total)
        # SRF 展宽敏感性：按波长分位选代表波段，避免 2 倍运行量
        if k in _sensitivity_band_indices(len(rows), int(sconf.get("srf_sensitivity_bands", 6))):
            wide = _run_one(exe, geom, alt["target_altitude_km"],
                            nominal, wl, fw, float(sconf.get("srf_sensitivity_scale", 2.0)))
            srf_rows.append({"output_band": k, "wavelength_nm": wl,
                             "nominal_fwhm_nm": fw,
                             "wide_fwhm_nm": wide["srf_fwhm_nm"],
                             "scenario": nominal["id"],
                             "xa_nominal": next(r["xa"] for r in lut_rows
                                                if r["output_band"] == k and r["nominal"]),
                             "xa_wide": wide["xa"],
                             "xb_nominal": next(r["xb"] for r in lut_rows
                                                if r["output_band"] == k and r["nominal"]),
                             "xb_wide": wide["xb"],
                             "xc_nominal": next(r["xc"] for r in lut_rows
                                                if r["output_band"] == k and r["nominal"]),
                             "xc_wide": wide["xc"],
                             "srf_mode_wide": wide["srf_mode"]})
    write_csv(lut_path, lut_rows)
    srf_path = work_dir / "srf_sensitivity.csv"
    if srf_rows:
        write_csv(srf_path, srf_rows)

    diffs = [r["formula_abs_diff"] for r in lut_rows
             if r["formula_abs_diff"] is not None]
    fallbacks = [r for r in lut_rows if r["srf_mode"] != "gaussian_2p5nm"]
    meta = {
        "lut_schema_version": LUT_SCHEMA_VERSION,
        "created_by": "ahsi_pipeline.sixs_lut.build_lut",
        "sixs_executable": str(exe),
        "sixs_self_test": check,
        "py6s_wrapper": "Py6S (pip); 6S binary from conda-forge 'sixs'",
        "geometry": geom, "target_altitude": alt,
        "band_count": len(rows), "scenario_count": len(scenarios),
        "nominal_scenario": nominal["id"], "scenarios": scenarios,
        "srf": {"model": "gaussian from center wavelength + FWHM",
                "step_um": SRF_STEP_UM, "grid_origin_um": .25,
                "half_width_fwhm": SRF_HALF_WIDTH_FWHM,
                "sensitivity_scale": sconf.get("srf_sensitivity_scale", 2.0),
                "caveat": "本景只有中心波长与 FWHM，无实测光谱响应函数；"
                          "高斯为显式近似，2xFWHM 敏感性见 srf_sensitivity.csv"},
        "reference_radiance": REFERENCE_RADIANCE,
        "inversion": "y = xa*L - xb ; rho = y/(1 + xc*y) ; L in W m-2 sr-1 um-1",
        "formula_reconciliation": {
            "max_abs_diff": max(diffs) if diffs else None,
            "mean_abs_diff": float(np.mean(diffs)) if diffs else None,
            "note": "6S 输出 xa/xb/xc 仅 5 位小数，与 6S 自报反射率的差异来自该打印精度；"
                    "不是模型分歧"},
        "srf_fallback_bands": len(fallbacks),
        "missing_inputs": [
            "无同步水汽观测（用假设情景区间代替）",
            "无同步气溶胶 AOT / AERONET（用假设情景区间代替）",
            "无同步臭氧观测（用假设情景区间代替）",
            "无实测逐波段光谱响应函数（用高斯近似）",
            "无实测地表反射率用于独立精度验证",
        ],
        "not_applied": ["地形光照校正", "邻近效应校正", "BRDF 校正",
                        "云/雪/阴影掩膜", "逐像元高程"],
        "deterministic": True,
        "lut_sha256": _sha256(lut_path),
        "srf_sensitivity_file": str(srf_path) if srf_rows else None,
    }
    save_json(work_dir / "sixs_lut_meta.json", meta)

    record = {
        "successful": True, "backend_method": "sixs_lut", "method": "sixs_lut",
        "stage": "build_lut", "lut_file": str(lut_path.resolve()),
        "meta_file": str((work_dir / "sixs_lut_meta.json").resolve()),
        "lut_sha256": meta["lut_sha256"],
        "band_count": len(rows), "lut_rows": len(lut_rows),
        "scenario_count": len(scenarios), "nominal_scenario": nominal["id"],
        "sixs_version": check["sixs_version"],
        "self_test": {k: v for k, v in check.items() if k != "stdout"},
        "formula_max_abs_diff": meta["formula_reconciliation"]["max_abs_diff"],
        "srf_fallback_bands": len(fallbacks),
        "srf_sensitivity_bands": len(srf_rows),
        "band_key": [[float(r["wavelength_nm"]), r["segment"],
                      int(r["source_band"]), int(r["smoke_output_band"])] for r in rows],
        "product_grade": PRODUCT_GRADE,
        "external_dependencies": [
            "6SV1.1 (conda-forge sixs, open-source local binary, no license)",
            "Py6S (pip, LGPL)"],
        "requires_external_binary": True,
        "deterministic": True,
        "atmospheric_correction_applied": True,
        "surface_reflectance_available": True,
        "quantitative_accuracy_independently_validated": False,
        "warning": ("6S model-estimated surface reflectance with assumed atmospheric "
                    "scenarios; NOT independently validated; do NOT treat as ground truth"),
    }
    save_json(work_dir / "sixs_stage_record.json", record)
    log.info("6S LUT 完成: %d 行, %d 情景, 公式最大偏差 %.3g",
             len(lut_rows), len(scenarios),
             meta["formula_reconciliation"]["max_abs_diff"] or float("nan"))
    return record


def _sensitivity_band_indices(n: int, count: int) -> set[int]:
    if n <= 0 or count <= 0:
        return set()
    k = min(count, n)
    return {int(round(i * (n - 1) / max(k - 1, 1))) + 1 for i in range(k)}


# --------------------------------------------------------------------------
# 第二段：反射率立方体写出（atmos.finalize 调用）
# --------------------------------------------------------------------------
def finalize_sixs_reflectance(scene, cfg: dict, scene_dir: Path,
                              cube_result: dict, stage_record: dict) -> dict:
    import csv

    if not stage_record.get("successful") or stage_record.get("backend_method") != "sixs_lut":
        raise ValueError("6S 阶段记录无效或方法不匹配")
    lut_path = Path(stage_record["lut_file"])
    if not lut_path.is_file():
        raise ValueError(f"6S LUT 文件不存在: {lut_path}")
    if _sha256(lut_path) != stage_record["lut_sha256"]:
        raise ValueError("6S LUT 内容已变化（疑似配置漂移）；禁止复用旧表，请用新 run_id 重算")
    with lut_path.open("r", encoding="utf-8-sig", newline="") as f:
        lut = list(csv.DictReader(f))
    rows = [r for r in cube_result["band_records"] if r["analysis_candidate"]]
    key = [[float(r["wavelength_nm"]), r["segment"],
            int(r["source_band"]), int(r["smoke_output_band"])] for r in rows]
    if key != stage_record["band_key"]:
        raise ValueError("候选波段表与 6S LUT 不一致；禁止复用旧表")
    nominal = stage_record["nominal_scenario"]
    coeff = {}
    scenarios = sorted(s["id"] for s in _scenarios(cfg))
    for r in lut:
        k = int(r["output_band"])
        scenario = r["scenario"]
        if not 1 <= k <= len(rows) or scenario not in scenarios:
            raise ValueError("6S LUT 包含意外的波段或情景")
        row = rows[k - 1]
        values = tuple(float(r[name]) for name in ("xa", "xb", "xc"))
        if (int(r["source_cube_band"]) != int(row["smoke_output_band"])
                or float(r["wavelength_nm"]) != float(row["wavelength_nm"])
                or float(r["fwhm_nm"]) != float(row["fwhm_nm"])
                or not _coefficients_finite(*values)):
            raise ValueError(f"6S LUT 波段 {k} 的来源、波长、FWHM 或系数无效")
        by_scenario = coeff.setdefault(k, {})
        if scenario in by_scenario:
            raise ValueError(f"6S LUT 波段 {k} 重复出现情景 {scenario}")
        by_scenario[scenario] = values
    if nominal not in scenarios or any(set(coeff.get(k, {})) != set(scenarios)
                                       for k in range(1, len(rows) + 1)):
        raise ValueError("6S LUT 缺少必要的波段或大气情景系数")

    out_dir = Path(scene_dir) / "sixs_reflectance_final"
    out_dir.mkdir(parents=True, exist_ok=True)
    nodata = float(cfg["grid"]["nodata"])
    with (rasterio.open(cube_result["coverage_mask"]) as d,
          rasterio.open(cube_result["radiance_file"]) as src):
        if (d.shape, d.crs, d.transform) != (src.shape, src.crs, src.transform):
            raise ValueError("6S 覆盖掩膜与辐亮度格网不一致")
        coverage = d.read(1, masked=True).filled(0) > 0
        units = src.tags().get("units", "")
        profile = {k: src.profile[k] for k in ("width", "height", "crs", "transform")}
    if not coverage.any():
        raise ValueError("6S 覆盖掩膜无有效观测")
    if " ".join(units.split()) != "W m-2 sr-1 um-1":
        raise ValueError(f"辐亮度单位不是 W m-2 sr-1 um-1（实际 {units!r}），"
                         "6S 系数不可直接使用")
    profile.update(driver="GTiff", count=len(rows), dtype="float32", nodata=nodata,
                   tiled=bool(cfg["runtime"]["tiled"]),
                   blockxsize=int(cfg["runtime"]["block_size"]),
                   blockysize=int(cfg["runtime"]["block_size"]),
                   interleave="band", compress=cfg["runtime"]["compression"],
                   predictor=3, zlevel=int(cfg["runtime"]["compress_level"]),
                   BIGTIFF="IF_SAFER")
    path = out_dir / f"{scene.scene_id}_6SV_surface_reflectance.tif"
    flags = np.zeros(coverage.shape, "uint8")
    available = coverage.copy()
    stats = []

    # 代表像元：覆盖区内按行列等间隔确定性取样，用于情景敏感性（非精度评估）
    yy, xx = np.where(coverage)
    sample_idx = (np.linspace(0, len(yy) - 1, 5).astype(int)
                  if len(yy) >= 5 else np.arange(len(yy)))
    sample_points = [(int(yy[i]), int(xx[i])) for i in sample_idx]
    sens_rows = []

    with rasterio.open(cube_result["radiance_file"]) as src, \
         rasterio.open(path, "w", **profile) as dst:
        dst.update_tags(
            product="6SV (6S V1.1) model-estimated surface reflectance; NOT independently validated",
            units="dimensionless",
            atmospheric_method="sixs_lut",
            nominal_scenario=nominal,
            aerosol_water_vapor_ozone_source="ASSUMED scenarios; NO synchronous measurement",
            spectral_response="Gaussian SRF from center wavelength + FWHM; 2.5 nm grid; NOT measured SRF",
            absolute_reflectance="MODEL_ESTIMATE",
            atmospheric_correction="APPLIED_6SV_RADIATIVE_TRANSFER",
            quantitative_accuracy_independently_validated="NO",
            terrain_illumination_correction="NOT_APPLIED",
            adjacency_correction="NOT_APPLIED",
            brdf_correction="NOT_APPLIED",
            usable_for="spectral shape, continuum removal, relative mineral diagnostics, "
                       "comparison under the SAME atmospheric scenarios",
            NOT_usable_for="absolute accuracy claims, cross-scene quantitative change detection "
                           "without re-deriving LUT and validation",
            pipeline="ahsi_pipeline")
        dst.update_tags(ns="ENVI", **{
            "wavelength": "{" + ",".join(str(r["wavelength_nm"]) for r in rows) + "}",
            "fwhm": "{" + ",".join(str(r["fwhm_nm"]) for r in rows) + "}",
            "wavelength units": "Nanometers",
            "data ignore value": str(nodata),
            "description": "6SV model-estimated surface reflectance; NOT independently validated"})
        for k, row in enumerate(rows, 1):
            L = src.read(int(row["smoke_output_band"]), masked=True).astype("float64").filled(np.nan)
            xa, xb, xc = coeff[k][nominal]
            rho = surface_reflectance_from_radiance(L, xa, xb, xc)
            valid = (coverage & np.isfinite(L) & (L > 0)
                     & np.isfinite(rho) & (np.abs(rho) <= np.finfo("float32").max))
            out = np.where(valid, rho, nodata).astype("float32")
            dst.write(out, k)
            dst.set_band_description(k, f"{row['wavelength_nm']:.3f} nm 6SV surface reflectance")
            dst.update_tags(k, wavelength_nm=row["wavelength_nm"], fwhm_nm=row["fwhm_nm"],
                            source_camera=row["segment"], source_native_band=row["source_band"],
                            sixs_xa=xa, sixs_xb=xb, sixs_xc=xc,
                            sixs_scenario=nominal,
                            product=PRODUCT_GRADE)
            flags[valid & (rho > 1)] |= 2
            flags[valid & (rho > 2)] |= 4
            sample = rho[::4, ::4][valid[::4, ::4]]
            p = (np.percentile(sample, [2, 50, 98]).tolist()
                 if sample.size else [None] * 3)
            stats.append({"output_band": k, "wavelength_nm": row["wavelength_nm"],
                          "sixs_xa": xa, "sixs_xb": xb, "sixs_xc": xc,
                          "valid_pixels": int(valid.sum()),
                          "reflectance_min": float(rho[valid].min()) if valid.any() else None,
                          "reflectance_max": float(rho[valid].max()) if valid.any() else None,
                          "sample_p02": p[0], "sample_median": p[1], "sample_p98": p[2],
                          "fraction_above_1": float(np.mean(rho[valid] > 1)) if valid.any() else None,
                          "fraction_above_2": float(np.mean(rho[valid] > 2)) if valid.any() else None})
            for j, (py, px) in enumerate(sample_points):
                if not valid[py, px]:
                    continue
                for sc in scenarios:
                    a, b, c = coeff[k][sc]
                    r_sc = float(surface_reflectance_from_radiance(L[py, px], a, b, c))
                    if not math.isfinite(r_sc):
                        continue
                    r_nom = float(rho[py, px])
                    sens_rows.append({"point": j, "row": py, "col": px,
                                      "output_band": k, "wavelength_nm": row["wavelength_nm"],
                                      "scenario": sc, "radiance": float(L[py, px]),
                                      "reflectance": r_sc,
                                      "delta_vs_nominal": r_sc - r_nom,
                                      "relative_delta": (r_sc - r_nom) / abs(r_nom)
                                      if abs(r_nom) > 1e-9 else None})
            available &= valid
            flags[coverage & ~valid] |= 1
        dst.update_tags(ns="ENVI", bbl="{" + ",".join(str(int(s["valid_pixels"] > 0)) for s in stats) + "}")
    write_csv(out_dir / "sixs_reflectance_band_quality.csv", stats)
    write_csv(out_dir / "scenario_sensitivity.csv", sens_rows)
    rel = [r["relative_delta"] for r in sens_rows if r["relative_delta"] is not None]
    sens = {"nominal_scenario": nominal, "scenarios": scenarios,
            "sample_points": len(sample_points), "rows": len(sens_rows),
            "relative_delta_p05": float(np.percentile(rel, 5)) if rel else None,
            "relative_delta_median": float(np.median(rel)) if rel else None,
            "relative_delta_p95": float(np.percentile(rel, 95)) if rel else None,
            "note": "假设情景导致的反射率离散，不是精度/误差评估；无同步大气观测"}
    save_json(out_dir / "scenario_sensitivity_summary.json", sens)

    # RGB 预览（与既有后端口径一致）
    with rasterio.open(path) as d:
        idx = [int(np.argmin([abs(r["wavelength_nm"] - w) for r in rows])) + 1
               for w in (650, 550, 470)]
        rgb = d.read(idx)
        valid = np.all(rgb != nodata, axis=0)
        img = np.zeros((*valid.shape, 4), "uint8")
        for i in range(3):
            v = rgb[i][valid]
            lo, hi = np.percentile(v, [2, 98]) if v.size else (0, 1)
            img[:, :, i] = np.clip((rgb[i] - lo) / max(hi - lo, 1e-8) * 255, 0, 255).astype("uint8")
        img[:, :, 3] = valid.astype("uint8") * 255
        im = Image.fromarray(img)
        im.thumbnail((1400, 1400))
        im.save(out_dir / "6SV_RGB_preview.png")

    mask_profile = {**profile, "count": 1, "dtype": "uint8", "nodata": 255, "predictor": 2}
    with rasterio.open(out_dir / "quality_flags.tif", "w", **mask_profile) as d:
        d.write(np.where(coverage, flags, 255).astype("uint8"), 1)
        d.update_tags(bit0="invalid / non-positive radiance / non-finite inversion",
                      bit1="6SV reflectance >1", bit2="6SV reflectance >2",
                      interpretation="screening flags only; >1 values are preserved")
    with rasterio.open(out_dir / "all_selected_bands_available.tif", "w", **mask_profile) as d:
        d.write(np.where(coverage, available, 255).astype("uint8"), 1)

    res_m = float(cfg["grid"]["resolution_m"])
    result = {
        "product_grade": PRODUCT_GRADE,
        "scene_id": scene.scene_id,
        "atmospheric_method": "sixs_lut",
        "nominal_scenario": nominal, "scenarios": scenarios,
        "reflectance_bands": len(rows), "reflectance_file": str(path.resolve()),
        "valid_reflectance_bands": sum(s["valid_pixels"] > 0 for s in stats),
        "lut_file": str(lut_path.resolve()),
        "sixs_version": stage_record["sixs_version"],
        "integer_scale_factor": 1.0,
        "coverage_area_km2": float(coverage.sum() * res_m ** 2 / 1e6),
        "all_bands_available_area_km2": float(available.sum() * res_m ** 2 / 1e6),
        "pixels_with_any_reflectance_above_1_fraction": float(np.mean((flags[coverage] & 2) > 0)),
        "pixels_with_any_reflectance_above_2_fraction": float(np.mean((flags[coverage] & 4) > 0)),
        "scenario_sensitivity": sens,
        "external_dependencies": stage_record["external_dependencies"],
        "deterministic": True,
        "atmospheric_correction_applied": True,
        "surface_reflectance_available": True,
        "quantitative_accuracy_independently_validated": False,
        "warning": ("6SV model-estimated surface reflectance; atmospheric parameters are "
                    "ASSUMED scenarios, not measurements; no independent accuracy validation; "
                    "no terrain/adjacency/BRDF correction"),
    }
    save_json(out_dir / "sixs_reflectance_quality_summary.json", result)
    return result
