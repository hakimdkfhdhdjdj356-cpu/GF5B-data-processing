"""无外部依赖的场景内辐射归一化：IARR 相对反射率。

IARR (Internal Average Relative Reflectance) = 逐波段把辐亮度除以该波段在
场景有效像元上的稳健均值。作用是压平太阳辐照度包络与大气透过率的一阶影响，
使剩余差异主要反映地表光谱形状，适合单景内做连续统去除、诊断吸收特征、
场景端元选择、MNF/PCA 与相对异常检测。

严格边界（必须在产物 tag、QC 报告、交付说明中同时声明）：
- 产物命名为 relative_reflectance，绝不叫 surface_reflectance；
- 不含任何大气校正，不是地表反射率，也没有绝对物理量纲；
- 逐波段缩放会改变光谱向量方向，不能与实验室绝对反射率光谱库直接做 SAM；
- 不能用于绝对反射率数值、跨景定量对比、时序变化检测。

实现只用 numpy + rasterio；无外部软件、无外部数据、无经验大气参数、无随机性，
同一输入 + 同一配置必然得到逐位相同的结果。

分两段执行：compute_factors 是耗时的第一遍全图扫描，结果落盘；
finalize_relative_reflectance 是第二遍写出，可独立断点续跑。
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import rasterio
from PIL import Image

from .core import save_json, write_csv

log = logging.getLogger(__name__)

MAX_STAT_SAMPLES = 200_000


def _candidate_rows(cube_result: dict) -> list[dict]:
    rows = [r for r in cube_result["band_records"] if r["analysis_candidate"]]
    if not rows:
        raise ValueError("没有候选波段，无法做 IARR 归一化；检查 bands.* 筛选规则")
    return rows


def _band_key(rows: list[dict]) -> list[list]:
    """因子与波段的绑定键；断点续跑时用它防止"换了配置却复用旧因子"。"""
    return [[float(r["wavelength_nm"]), r["segment"], int(r["source_band"]),
             int(r["smoke_output_band"])] for r in rows]


def _valid_like(raw: np.ndarray, coverage: np.ndarray, nodata: float,
                require_positive: bool) -> np.ndarray:
    valid = coverage & np.isfinite(raw) & (raw != nodata)
    if require_positive:
        valid &= raw > 0
    return valid


def _robust_factor(values: np.ndarray, iarr: dict) -> float:
    """等间隔确定性子采样求分位，避免整图排序，也避免任何随机性。"""
    lo, hi = iarr.get("robust_pct", [2, 98])
    if lo > 0 or hi < 100:
        step = max(1, values.size // MAX_STAT_SAMPLES)
        sample = values[::step]
        low_v, high_v = np.percentile(sample, [lo, hi])
        clipped = values[(values >= low_v) & (values <= high_v)]
        if clipped.size:
            values = clipped
    if iarr.get("normalization", "scene_mean") == "scene_median":
        return float(np.median(values))
    return float(np.mean(values))


def compute_factors(scene, cfg: dict, cube_result: dict, work_dir: Path) -> dict:
    """第一遍：逐波段求稳健归一化因子，落盘 iarr_factors.json。"""
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    iarr = cfg["atmospheric"]["iarr"]
    nodata = float(cfg["grid"]["nodata"])
    require_positive = bool(iarr.get("require_positive_radiance", True))
    min_valid = int(iarr["min_valid_pixels"])
    rows = _candidate_rows(cube_result)

    with rasterio.open(cube_result["coverage_mask"]) as d:
        coverage = d.read(1) > 0

    factors, stats = [], []
    with rasterio.open(cube_result["radiance_file"]) as src:
        for k, row in enumerate(rows, 1):
            raw = src.read(int(row["smoke_output_band"]))
            valid = _valid_like(raw, coverage, nodata, require_positive)
            n = int(valid.sum())
            rec = {"output_band": k, "source_cube_band": int(row["smoke_output_band"]),
                   "wavelength_nm": row["wavelength_nm"], "camera": row["segment"],
                   "native_band": row["source_band"], "valid_pixels": n}
            if n < min_valid:
                rec.update(normalization_factor=None, good=False,
                           reason=f"valid_pixels < min_valid_pixels ({min_valid})")
                factors.append(None)
            else:
                mu = _robust_factor(raw[valid], iarr)
                if not np.isfinite(mu) or mu <= 0:
                    rec.update(normalization_factor=mu, good=False,
                               reason="normalization factor not finite or <= 0")
                    factors.append(None)
                else:
                    rec.update(normalization_factor=mu, good=True, reason="")
                    factors.append(mu)
            stats.append(rec)

    record = {
        "successful": True, "method": "iarr", "stage": "compute_factors",
        "normalization": iarr.get("normalization", "scene_mean"),
        "robust_pct": list(iarr.get("robust_pct", [2, 98])),
        "min_valid_pixels": min_valid, "require_positive_radiance": require_positive,
        "radiance_file": str(Path(cube_result["radiance_file"]).resolve()),
        "coverage_mask": str(Path(cube_result["coverage_mask"]).resolve()),
        "band_key": _band_key(rows), "band_count": len(rows),
        "good_bands": sum(f is not None for f in factors),
        "factors": factors, "band_stats": stats,
        "external_dependencies": [], "deterministic": True,
        "product_grade": "relative_reflectance_NOT_surface_reflectance",
        "atmospheric_correction_applied": False,
        "surface_reflectance_available": False,
    }
    save_json(work_dir / "iarr_factors.json", record)
    write_csv(work_dir / "iarr_band_validity.csv", stats)
    log.info("IARR 归一化因子: %d/%d 波段可用", record["good_bands"], len(rows))
    return record


def finalize_relative_reflectance(scene, cfg: dict, out_dir: Path,
                                  cube_result: dict, factors_record: dict) -> dict:
    """第二遍：辐亮度 / 逐波段因子 -> 相对反射率立方体 + 质量表 + 旗标 + 预览。"""
    if not factors_record.get("successful") or factors_record.get("method") != "iarr":
        raise ValueError("IARR 因子记录无效或方法不匹配")
    rows = _candidate_rows(cube_result)
    if _band_key(rows) != factors_record["band_key"]:
        raise ValueError("候选波段表与已保存的 IARR 因子不一致（配置或输入已改变）；"
                         "请用新的 run_id 或 --force 重算，禁止复用旧因子")
    factors = factors_record["factors"]
    stats = [dict(s) for s in factors_record["band_stats"]]
    iarr = cfg["atmospheric"]["iarr"]
    nodata = float(cfg["grid"]["nodata"])
    require_positive = bool(factors_record["require_positive_radiance"])

    out_dir = Path(out_dir) / "relative_reflectance_final"
    out_dir.mkdir(parents=True, exist_ok=True)
    with rasterio.open(cube_result["coverage_mask"]) as d:
        coverage = d.read(1) > 0

    with rasterio.open(cube_result["radiance_file"]) as src:
        profile = {k: src.profile[k] for k in ("width", "height", "crs", "transform")}
    profile.update(driver="GTiff", count=len(rows), dtype="float32", nodata=nodata,
                   tiled=bool(cfg["runtime"]["tiled"]),
                   blockxsize=int(cfg["runtime"]["block_size"]),
                   blockysize=int(cfg["runtime"]["block_size"]),
                   interleave="band", compress=cfg["runtime"]["compression"],
                   predictor=3, zlevel=int(cfg["runtime"]["compress_level"]),
                   BIGTIFF="IF_SAFER")
    path = out_dir / f"{scene.scene_id}_relative_reflectance.tif"
    flags = np.zeros(coverage.shape, "uint8")
    available = coverage.copy()

    with rasterio.open(cube_result["radiance_file"]) as src, rasterio.open(path, "w", **profile) as dst:
        dst.update_tags(
            product="IARR relative reflectance; NOT surface reflectance; NOT atmospherically corrected",
            units="dimensionless (scene-mean normalized per band)",
            atmospheric_method="iarr",
            normalization=str(factors_record["normalization"]),
            robust_percentile=str(factors_record["robust_pct"]),
            absolute_reflectance="NO",
            atmospheric_correction="NOT_APPLIED",
            quantitative_accuracy_independently_validated="NO",
            usable_for="within-scene spectral shape: continuum removal, absorption features, scene endmembers",
            NOT_usable_for="absolute reflectance, cross-scene quantitative comparison, direct SAM against lab spectra",
            pipeline="ahsi_pipeline")
        dst.update_tags(ns="ENVI", **{
            "wavelength": "{" + ",".join(str(r["wavelength_nm"]) for r in rows) + "}",
            "fwhm": "{" + ",".join(str(r["fwhm_nm"]) for r in rows) + "}",
            "wavelength units": "Nanometers",
            "bbl": "{" + ",".join("1" if f is not None else "0" for f in factors) + "}",
            "data ignore value": str(nodata),
            "description": "IARR relative reflectance; NOT surface reflectance"})
        for k, (row, mu) in enumerate(zip(rows, factors), 1):
            raw = src.read(int(row["smoke_output_band"]))
            out = np.full(raw.shape, nodata, "float32")
            valid = (_valid_like(raw, coverage, nodata, require_positive)
                     if mu is not None else np.zeros_like(coverage))
            if mu is not None:
                rel = raw.astype("float32") / np.float32(mu)
                out[valid] = rel[valid]
                flags[valid & (rel > 2)] |= 2
                flags[valid & (rel > 3)] |= 4
                sample = rel[::4, ::4][valid[::4, ::4]]
                p = np.percentile(sample, [2, 50, 98]).tolist() if sample.size else [None] * 3
                stats[k - 1].update(relative_min=float(rel[valid].min()) if valid.any() else None,
                                    relative_max=float(rel[valid].max()) if valid.any() else None,
                                    relative_sample_p02=p[0], relative_sample_p50=p[1],
                                    relative_sample_p98=p[2],
                                    fraction_above_2=float(np.mean(rel[valid] > 2)) if valid.any() else None,
                                    fraction_above_3=float(np.mean(rel[valid] > 3)) if valid.any() else None)
            else:
                stats[k - 1].update(relative_min=None, relative_max=None,
                                    relative_sample_p02=None, relative_sample_p50=None,
                                    relative_sample_p98=None, fraction_above_2=None,
                                    fraction_above_3=None)
            dst.write(out, k)
            dst.set_band_description(k, f"{row['wavelength_nm']:.3f} nm IARR relative reflectance")
            dst.update_tags(k, wavelength_nm=row["wavelength_nm"], fwhm_nm=row["fwhm_nm"],
                            source_camera=row["segment"], source_native_band=row["source_band"],
                            normalization_factor=mu if mu is not None else "BAD_BAND",
                            product="relative_reflectance_NOT_surface_reflectance")
            available &= valid
            flags[coverage & ~valid] |= 1
    write_csv(out_dir / "relative_reflectance_band_quality.csv", stats)

    with rasterio.open(path) as d:
        idx = [int(np.argmin([abs(r["wavelength_nm"] - w) for r in rows])) + 1 for w in (650, 550, 470)]
        rgb = d.read(idx)
        valid = np.all(rgb != nodata, axis=0)
        with rasterio.open(out_dir / f"{scene.scene_id}_RGB_relative_reflectance.tif", "w",
                           **{**profile, "count": 3}) as o:
            o.write(rgb)
        img = np.zeros((*valid.shape, 4), "uint8")
        for i in range(3):
            v = rgb[i][valid]
            lo, hi = np.percentile(v, [2, 98]) if v.size else (0, 1)
            img[:, :, i] = np.clip((rgb[i] - lo) / max(hi - lo, 1e-8) * 255, 0, 255).astype("uint8")
        img[:, :, 3] = valid.astype("uint8") * 255
        im = Image.fromarray(img)
        im.thumbnail((1400, 1400))
        im.save(out_dir / "IARR_RGB_preview.png")

    mask_profile = {**profile, "count": 1, "dtype": "uint8", "nodata": 255, "predictor": 2}
    with rasterio.open(out_dir / "quality_flags.tif", "w", **mask_profile) as d:
        d.write(np.where(coverage, flags, 255).astype("uint8"), 1)
        d.update_tags(bit0="some candidate band invalid or bad", bit1="relative value >2",
                      bit2="relative value >3",
                      interpretation="screening flags only; values are NOT surface reflectance")
    with rasterio.open(out_dir / "all_selected_bands_available.tif", "w", **mask_profile) as d:
        d.write(np.where(coverage, available, 255).astype("uint8"), 1)

    res_m = float(cfg["grid"]["resolution_m"])
    result = {
        "product_grade": "relative_reflectance_NOT_surface_reflectance",
        "atmospheric_method": "iarr",
        "normalization": factors_record["normalization"],
        "robust_pct": factors_record["robust_pct"],
        "scene_id": scene.scene_id,
        "reflectance_bands": len(rows), "good_bands": factors_record["good_bands"],
        "bad_band_wavelengths_nm": [s["wavelength_nm"] for s in stats if not s["good"]],
        "reflectance_file": str(path.resolve()),
        "integer_scale_factor": 1.0,   # 无整数缩放：直接 float32 相对值，禁止任何二次缩放
        "coverage_area_km2": float(coverage.sum() * res_m ** 2 / 1e6),
        "all_bands_available_area_km2": float(available.sum() * res_m ** 2 / 1e6),
        "value_above_2_fraction": float(np.mean((flags[coverage] & 2) > 0)),
        "value_above_3_fraction": float(np.mean((flags[coverage] & 4) > 0)),
        "external_dependencies": [],
        "deterministic": True,
        "atmospheric_correction_applied": False,
        "surface_reflectance_available": False,
        "warning": "IARR relative reflectance is NOT surface reflectance; "
                   "per-band normalization changes spectral vector direction, "
                   "do NOT match directly against absolute laboratory spectra",
    }
    save_json(out_dir / "relative_reflectance_quality_summary.json", result)
    return result
