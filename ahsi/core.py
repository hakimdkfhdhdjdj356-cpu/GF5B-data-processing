"""数值核心：逐波段辐射定标 + 显式 RPC/DEM 正射 + 成果写出。

数值逻辑与 2026-11-06 已验证的 process.py / finalize_reflectance.py 保持一致，
差别只在：全部参数来自配置、输出按 run_id 隔离、可重入、不硬编码卫星与结论。
禁止用普通仿射变换替代 RPC 处理这些无地理参考的 L1 TIFF。
"""
from __future__ import annotations

import csv
import json
import logging
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.transform import RPCTransformer, from_origin
from rasterio.warp import calculate_default_transform, reproject, transform as warp_coordinates
from PIL import Image

from .config import aux_path
from .scene import Scene

log = logging.getLogger(__name__)

RESAMPLING = {"nearest": Resampling.nearest, "bilinear": Resampling.bilinear,
              "cubic": Resampling.cubic, "average": Resampling.average}


def save_json(path: Path, value) -> None:
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with Path(path).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader()
        w.writerows(rows)


def _target_crs(scene: Scene, cfg: dict) -> str:
    if not cfg["grid"]["auto_utm"]:
        return cfg["grid"]["crs"]
    lon = float(scene.meta.get("center_lon") or 0)
    lat = float(scene.meta.get("center_lat") or 0)
    zone = int((lon + 180) // 6) + 1
    return f"EPSG:{32600 + zone if lat >= 0 else 32700 + zone}"


def _rpc_options(cfg: dict) -> dict:
    dem = aux_path(cfg, "rpc_dem")
    geoid = aux_path(cfg, "geoid_grid")
    if not dem.is_file():
        raise FileNotFoundError(f"RPC 正射所需 DEM 不存在: {dem}")
    opts = {"RPC_DEM": str(dem)}
    if cfg["grid"]["apply_vdatum_shift"]:
        if not geoid.is_file():
            raise FileNotFoundError(f"EGM2008 高程网格不存在: {geoid}")
        opts["RPC_DEM_SRS"] = f"+proj=longlat +datum=WGS84 +geoidgrids={geoid} +type=crs"
        opts["RPC_DEM_APPLY_VDATUM_SHIFT"] = "YES"
    return opts


def _grid(scene: Scene, cfg: dict, target_crs: str, opts: dict):
    res = float(cfg["grid"]["resolution_m"])
    requested = cfg["grid"].get("bounds_m")
    if requested is not None:
        if len(requested) != 4 or not all(math.isfinite(float(value)) for value in requested):
            raise ValueError("grid.bounds_m must contain four finite target CRS coordinates")
        left, bottom, right, top = map(float, requested)
        if left >= right or bottom >= top:
            raise ValueError("grid.bounds_m must be ordered west, south, east, north")
        bounds = [(left, bottom, right, top)]
    else:
        bounds = []
        for seg in scene.segments:
            with rasterio.open(seg.image) as ds:
                tr, w, h = calculate_default_transform("EPSG:4326", target_crs, ds.width, ds.height,
                                                       rpcs=ds.rpcs, resolution=res, **opts)
                bounds.append((tr.c, tr.f + h * tr.e, tr.c + w * tr.a, tr.f))
    left = math.floor(min(b[0] for b in bounds) / res) * res
    bottom = math.floor(min(b[1] for b in bounds) / res) * res
    right = math.ceil(max(b[2] for b in bounds) / res) * res
    top = math.ceil(max(b[3] for b in bounds) / res) * res
    width = int(round((right - left) / res))
    height = int(round((top - bottom) / res))
    minimum = 0 if requested is not None else 100
    if not (minimum < width < 20000 and minimum < height < 20000):
        raise ValueError(f"目标网格尺寸异常 {width}x{height}，检查 RPC / CRS / 分辨率配置")
    return from_origin(left, top, res, res), width, height, [left, bottom, right, top]


def _validity_mask(ds, chunk: int) -> np.ndarray:
    """整相机光谱全为 0 才算无效；单波段 0 值可能是真实低信号，不能当 NoData。"""
    valid = np.zeros(ds.shape, bool)
    for start in range(1, ds.count + 1, chunk):
        idx = list(range(start, min(start + chunk, ds.count + 1)))
        valid |= np.any(ds.read(idx) != 0, axis=0)
    return valid


def _build_vrt(path: Path, cube: Path, rows: list[dict], width: int, height: int,
               transform, crs_wkt: str, nodata: float) -> None:
    root = ET.Element("VRTDataset", rasterXSize=str(width), rasterYSize=str(height))
    ET.SubElement(root, "SRS").text = crs_wkt
    ET.SubElement(root, "GeoTransform").text = ",".join(
        map(str, [transform.c, transform.a, transform.b, transform.f, transform.d, transform.e]))
    md = ET.SubElement(root, "Metadata", domain="ENVI")
    for key, val in [("wavelength", "{" + ",".join(str(r["wavelength_nm"]) for r in rows) + "}"),
                     ("fwhm", "{" + ",".join(str(r["fwhm_nm"]) for r in rows) + "}"),
                     ("wavelength units", "Nanometers")]:
        ET.SubElement(md, "MDI", key=key).text = val
    for index, row in enumerate(rows, 1):
        band = ET.SubElement(root, "VRTRasterBand", dataType="Float32", band=str(index))
        ET.SubElement(band, "NoDataValue").text = str(nodata)
        ET.SubElement(band, "Description").text = f"{row['wavelength_nm']:.3f} nm; {row['segment']}; radiance"
        src = ET.SubElement(band, "SimpleSource")
        ET.SubElement(src, "SourceFilename", relativeToVRT="1").text = cube.name
        ET.SubElement(src, "SourceBand").text = str(row["output_band"])
        ET.SubElement(src, "SrcRect", xOff="0", yOff="0", xSize=str(width), ySize=str(height))
        ET.SubElement(src, "DstRect", xOff="0", yOff="0", xSize=str(width), ySize=str(height))
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def orthorectify(scene: Scene, cfg: dict, out_dir: Path,
                 max_bands: int | None = None, progress=log.info) -> dict:
    """L1 DN -> 辐亮度立方体（GeoTIFF）+ 波段质量表 + 覆盖掩膜 + RGB + 候选波段 VRT。"""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    target_crs = _target_crs(scene, cfg)
    opts = _rpc_options(cfg)
    transform, width, height, bounds = _grid(scene, cfg, target_crs, opts)
    nodata = float(cfg["grid"]["nodata"])
    chunk = int(cfg["runtime"]["band_chunk"])
    resampling = RESAMPLING[cfg["grid"]["resampling"]]

    records = [dict(r) for r in scene.band_table]
    if max_bands:
        records = records[:max_bands]
        progress(f"SMOKE 模式：仅输出前 {len(records)} 个波段")
    lookup = {(r["segment"], r["source_band"]): r for r in records}
    for i, r in enumerate(records, 1):
        r["smoke_output_band"] = i

    profile = {"driver": "GTiff", "width": width, "height": height, "count": len(records),
               "crs": target_crs, "transform": transform, "dtype": "float32", "nodata": nodata,
               "tiled": bool(cfg["runtime"]["tiled"]),
               "blockxsize": int(cfg["runtime"]["block_size"]),
               "blockysize": int(cfg["runtime"]["block_size"]),
               "interleave": "band", "compress": cfg["runtime"]["compression"],
               "predictor": 3, "zlevel": int(cfg["runtime"]["compress_level"]),
               "BIGTIFF": "IF_SAFER"}
    cube = out_dir / f"{scene.scene_id}_radiance.tif"
    common = np.ones((height, width), bool)
    geometry_checks, calibration_checks = [], []

    with rasterio.open(cube, "w", **profile) as dst:
        dst.update_tags(product="At-sensor spectral radiance; NOT surface reflectance",
                        units=cfg["sensors"]["calibration"]["units"],
                        calibration=cfg["sensors"]["calibration"]["formula"],
                        geometry=f"Explicit RPC + {cfg['grid']['rpc_dem']}; {cfg['grid']['resampling']}",
                        source_scene=scene.scene_id,
                        spectral_band_order="ascending wavelength; camera identity retained",
                        atmospheric_correction="NOT_APPLIED",
                        terrain_radiometric_correction="NOT_APPLIED",
                        pipeline="ahsi_pipeline")
        dst.update_tags(ns="ENVI",
                        wavelength="{" + ",".join(str(r["wavelength_nm"]) for r in records) + "}",
                        fwhm="{" + ",".join(str(r["fwhm_nm"]) for r in records) + "}",
                        wavelength_units="Nanometers")
        for r in records:
            i = r["smoke_output_band"]
            dst.set_band_description(i, f"{r['wavelength_nm']:.3f} nm {r['segment']}")
            dst.update_tags(i, wavelength_nm=r["wavelength_nm"], fwhm_nm=r["fwhm_nm"],
                            source_band=r["source_band"], source_segment=r["segment"],
                            source_file=r["source_file"], gain=r["gain"], offset=r["offset"],
                            analysis_candidate=str(r["analysis_candidate"]),
                            units=cfg["sensors"]["calibration"]["units"])

        for seg in scene.segments:
            with rasterio.open(seg.image) as ds:
                valid = _validity_mask(ds, chunk)
                mask = np.zeros((height, width), "uint8")
                reproject(valid.astype("uint8"), mask, rpcs=ds.rpcs, src_crs="EPSG:4326",
                          dst_crs=target_crs, dst_transform=transform, src_nodata=0, dst_nodata=0,
                          resampling=Resampling.nearest, **opts)
                common &= mask > 0
                with RPCTransformer(ds.rpcs, **opts) as tr:
                    if cfg["grid"].get("bounds_m") is not None:
                        left, bottom, right, top = bounds
                        lon, lat = warp_coordinates(target_crs, "EPSG:4326",
                                                    [(left + right) / 2], [(bottom + top) / 2])
                        src_row, src_col = tr.rowcol(lon, lat, op=lambda v: v)
                    else:
                        src_row = [ds.height // 2 + .5]
                        src_col = [ds.width // 2 + .5]
                    lon, lat = tr.xy(src_row, src_col, offset="ul")
                    row, col = tr.rowcol(lon, lat, op=lambda v: v)
                    residual = float(math.hypot(row[0] - src_row[0], col[0] - src_col[0]))
                    if not math.isfinite(residual):
                        raise ValueError(f"{seg.name}: RPC/DEM roundtrip failed inside target grid")
                geometry_checks.append({"segment": seg.name, "rpc_roundtrip_pixels": residual,
                                        "independent_ground_accuracy_measured": False})

                for start in range(1, seg.count + 1, chunk):
                    indexes = [b for b in range(start, min(start + chunk, seg.count + 1))
                               if (seg.name, b) in lookup]
                    if not indexes:
                        continue
                    native = ds.read(indexes).astype("float32")
                    native[:, ~valid] = np.nan
                    warped = np.full((len(indexes), height, width), np.nan, "float32")
                    reproject(native, warped, rpcs=ds.rpcs, src_crs="EPSG:4326", dst_crs=target_crs,
                              dst_transform=transform, src_nodata=np.nan, dst_nodata=np.nan,
                              resampling=resampling, num_threads=2, warp_mem_limit=256, **opts)
                    for i, band in enumerate(indexes):
                        rec = lookup[(seg.name, band)]
                        vals = native[i][valid]
                        rec.update(native_min_DN=float(vals.min()), native_max_DN=float(vals.max()),
                                   native_mean_DN=float(vals.mean()),
                                   zero_DN_fraction=float(np.mean(vals == 0)),
                                   negative_DN_fraction=float(np.mean(vals < 0)))
                        finite = np.isfinite(warped[i])
                        rad = warped[i] * rec["gain"] + rec["offset"]
                        if finite.any():
                            rec["orthorectified_mean_radiance"] = float(rad[finite].mean())
                            rec["orthorectified_valid_pixels"] = int(finite.sum())
                        dst.write(np.where(finite, rad, nodata).astype("float32"), rec["smoke_output_band"])
                        if band == indexes[0]:
                            loc = (min(ds.height // 2, ds.height - 1), min(ds.width // 2, ds.width - 1))
                            dn = float(native[i, loc[0], loc[1]])
                            calibration_checks.append({"segment": seg.name, "source_band": band,
                                                       "native_DN": dn,
                                                       "computed_radiance": dn * rec["gain"] + rec["offset"]})
                    progress(f"{scene.scene_id} {seg.name} {indexes[-1]}/{seg.count}")
                    del native, warped

    write_csv(out_dir / "band_parameters_and_quality.csv", records)
    selected = [r for r in records if r["analysis_candidate"]]
    if selected:
        _build_vrt(out_dir / f"{scene.scene_id}_analysis_candidates.vrt", cube,
                   [{**r, "output_band": r["smoke_output_band"]} for r in selected],
                   width, height, transform, CRS.from_string(target_crs).to_wkt(), nodata)

    with rasterio.open(cube) as ds:
        rgb_rows = [min(records, key=lambda r: abs(r["wavelength_nm"] - w)) for w in (650, 550, 470)]
        rgb = ds.read([r["smoke_output_band"] for r in rgb_rows])
        rgb_valid = np.all(rgb != nodata, axis=0) & common
        rgb[:, ~rgb_valid] = nodata
        with rasterio.open(out_dir / f"{scene.scene_id}_RGB_radiance.tif", "w",
                           **{**profile, "count": 3}) as d:
            d.write(rgb)
            d.update_tags(interpretation="Three actual radiance bands; no reflectance correction")
        display = np.zeros((height, width, 4), "uint8")
        stretches = []
        for i in range(3):
            lo, hi = np.percentile(rgb[i][rgb_valid], [2, 98]) if rgb_valid.any() else (0, 1)
            stretches.append([float(lo), float(hi)])
            display[:, :, i] = np.clip((rgb[i] - lo) / max(hi - lo, 1e-8) * 255, 0, 255).astype("uint8")
        display[:, :, 3] = rgb_valid.astype("uint8") * 255
        img = Image.fromarray(display)
        img.thumbnail((1600, 1600))
        img.save(out_dir / "RGB_preview.png")

        readback_ok = True
        if cfg["qc"]["require_complete_readback"]:
            for start in range(1, ds.count + 1, chunk):
                vals = ds.read(list(range(start, min(start + chunk, ds.count + 1))))
                if not np.isfinite(vals).all():
                    readback_ok = False
                    break
        if selected:
            with rasterio.open(out_dir / f"{scene.scene_id}_analysis_candidates.vrt") as v:
                assert v.count == len(selected)

    with rasterio.open(out_dir / "joint_camera_coverage_mask.tif", "w",
                       **{**profile, "count": 1, "dtype": "uint8", "nodata": 0, "predictor": 2}) as d:
        d.write(common.astype("uint8"), 1)
        d.update_tags(meaning="Both cameras spatially covered; clouds/snow/shadow NOT verified")

    result = {
        "scene_id": scene.scene_id, "satellite": scene.satellite,
        "source_level": scene.meta.get("product_level"), "start_time": scene.meta.get("start_time"),
        "source_bands": scene.total_bands, "output_bands": len(records),
        "candidate_bands": len(selected), "smoke_mode": bool(max_bands),
        "wavelength_min_nm": records[0]["wavelength_nm"], "wavelength_max_nm": records[-1]["wavelength_nm"],
        "crs": target_crs, "resolution_m": float(cfg["grid"]["resolution_m"]),
        "width": width, "height": height, "bounds_m": bounds,
        "joint_camera_area_km2": float(common.sum() * float(cfg["grid"]["resolution_m"]) ** 2 / 1e6),
        "radiance_file": str(cube.resolve()), "coverage_mask": str((out_dir / "joint_camera_coverage_mask.tif").resolve()),
        "vrt_file": str((out_dir / f"{scene.scene_id}_analysis_candidates.vrt").resolve()) if selected else None,
        "geometry_checks": geometry_checks, "calibration_samples": calibration_checks,
        "complete_raster_readback_passed": readback_ok,
        "radiometric_calibration_applied": True, "rpc_dem_orthorectification_applied": True,
        "atmospheric_correction_applied": False, "surface_reflectance_available": False,
        "ground_control_accuracy_measured": False, "destriping_applied": False,
        "rgb_actual_wavelengths_nm": [r["wavelength_nm"] for r in rgb_rows],
        "rgb_display_stretches": stretches,
        "band_records": records,
    }
    save_json(out_dir / "radiance_summary.json", result)
    return result


def export_envi(scene: Scene, cfg: dict, out_dir: Path, cube_result: dict) -> Path | None:
    """从候选波段 VRT 导出标准 ENVI 辐亮度，作为 ENVI QUAC 的输入。"""
    vrt = Path(cube_result["vrt_file"]) if cube_result.get("vrt_file") else None
    if not vrt or not vrt.is_file():
        return None
    out_dir = Path(out_dir) / "ENVI_input"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{scene.scene_id}_calibrated_analysis_radiance.dat"
    rows = [r for r in cube_result["band_records"] if r["analysis_candidate"]]
    with rasterio.open(vrt) as src:
        profile = {k: src.profile[k] for k in ("width", "height", "crs", "transform", "count", "dtype", "nodata")}
        profile.update(driver="ENVI", interleave="bsq")
        with rasterio.open(path, "w", **profile) as dst:
            chunk = int(cfg["runtime"]["band_chunk"])
            for start in range(1, src.count + 1, chunk):
                bands = list(range(start, min(start + chunk, src.count + 1)))
                dst.write(src.read(bands), indexes=bands)
            dst.update_tags(ns="ENVI", **{
                "wavelength": "{" + ",".join(str(r["wavelength_nm"]) for r in rows) + "}",
                "fwhm": "{" + ",".join(str(r["fwhm_nm"]) for r in rows) + "}",
                "wavelength units": "Nanometers",
                "data ignore value": str(cfg["grid"]["nodata"]),
                "description": "Calibrated RPC/DEM orthorectified radiance; NOT surface reflectance"})
    with rasterio.open(path) as d:
        assert d.count == len(rows)
        assert len(d.tags(ns="ENVI")["wavelength"].strip("{}").split(",")) == d.count
    return path


def finalize_reflectance(scene: Scene, cfg: dict, out_dir: Path, quac_dat: Path,
                         cube_result: dict) -> dict:
    """QUAC uint16 -> float32 反射率立方体；只缩放一次；保留 >1 值并标记。"""
    out_dir = Path(out_dir) / "reflectance_final"
    out_dir.mkdir(parents=True, exist_ok=True)
    nodata = float(cfg["grid"]["nodata"])
    with rasterio.open(cube_result["coverage_mask"]) as d:
        coverage = d.read(1) > 0
    rows = [r for r in cube_result["band_records"] if r["analysis_candidate"]]
    stats, flags = [], np.zeros(coverage.shape, "uint8")
    available = coverage.copy()

    with rasterio.open(quac_dat) as src:
        envi = src.tags(ns="ENVI")
        scale = float(envi.get("reflectance_scale_factor", cfg["atmospheric"]["reflectance_scale"]))
        bbl = [int(x) for x in envi.get("bbl", "").strip("{} ").replace(",", " ").split()] or [1] * src.count
        if src.count != len(rows) or len(bbl) != src.count:
            raise ValueError("QUAC 输出波段数与候选波段表不一致")
        profile = {k: src.profile[k] for k in ("width", "height", "crs", "transform", "count")}
        profile.update(driver="GTiff", dtype="float32", nodata=nodata, tiled=True,
                       blockxsize=256, blockysize=256, interleave="band",
                       compress=cfg["runtime"]["compression"], predictor=3, BIGTIFF="IF_SAFER")
        path = out_dir / f"{scene.scene_id}_surface_reflectance.tif"
        with rasterio.open(path, "w", **profile) as dst:
            dst.update_tags(product="ENVI QUAC estimated surface reflectance", units="dimensionless",
                            scale_applied=f"integer / {scale:g}; ALREADY float reflectance, DO NOT scale again",
                            atmospheric_method=f"{cfg['atmospheric']['method']}, {cfg['atmospheric']['sensor']}",
                            quantitative_accuracy_independently_validated="NO",
                            terrain_illumination_correction="NOT_APPLIED")
            dst.update_tags(ns="ENVI", **{k: envi[k] for k in ("wavelength", "fwhm", "wavelength_units", "bbl") if k in envi})
            for i, row in enumerate(rows, 1):
                raw = src.read(i)
                valid = coverage & (raw != src.nodata) & bool(bbl[i - 1])
                refl = raw.astype("float32") / scale
                dst.write(np.where(valid, refl, nodata).astype("float32"), i)
                dst.set_band_description(i, f"{row['wavelength_nm']:.3f} nm QUAC reflectance")
                dst.update_tags(i, wavelength_nm=row["wavelength_nm"], fwhm_nm=row["fwhm_nm"],
                                source_camera=row["segment"], source_native_band=row["source_band"],
                                quac_bbl=bbl[i - 1])
                available &= valid
                flags[coverage & ~valid] |= 1
                flags[valid & (refl > 1)] |= 2
                flags[valid & (refl > 2)] |= 4
                sample = refl[::4, ::4][valid[::4, ::4]]
                p = np.percentile(sample, [2, 50, 98]).tolist() if sample.size else [None] * 3
                stats.append({"output_band": i, "wavelength_nm": row["wavelength_nm"],
                              "camera": row["segment"], "native_band": row["source_band"],
                              "quac_bbl": bbl[i - 1], "valid_pixels": int(valid.sum()),
                              "reflectance_min": float(refl[valid].min()) if valid.any() else None,
                              "reflectance_max": float(refl[valid].max()) if valid.any() else None,
                              "sample_p02": p[0], "sample_median": p[1], "sample_p98": p[2],
                              "fraction_above_1": float(np.mean(refl[valid] > 1)) if valid.any() else None,
                              "fraction_above_2": float(np.mean(refl[valid] > 2)) if valid.any() else None})
    write_csv(out_dir / "reflectance_band_quality.csv", stats)

    with rasterio.open(path) as d:
        idx = [int(np.argmin([abs(r["wavelength_nm"] - w) for r in rows])) + 1 for w in (650, 550, 470)]
        rgb = d.read(idx)
        valid = np.all(rgb != nodata, axis=0)
        with rasterio.open(out_dir / f"{scene.scene_id}_RGB_reflectance.tif", "w",
                           **{**profile, "count": 3}) as o:
            o.write(rgb)
        img = np.zeros((*valid.shape, 4), "uint8")
        for i in range(3):
            lo, hi = np.percentile(rgb[i][valid], [2, 98]) if valid.any() else (0, 1)
            img[:, :, i] = np.clip((rgb[i] - lo) / max(hi - lo, 1e-8) * 255, 0, 255).astype("uint8")
        img[:, :, 3] = valid.astype("uint8") * 255
        im = Image.fromarray(img)
        im.thumbnail((1400, 1400))
        im.save(out_dir / "QUAC_RGB_preview.png")

    mask_profile = {**profile, "count": 1, "dtype": "uint8", "nodata": 255, "predictor": 2}
    with rasterio.open(out_dir / "quality_flags.tif", "w", **mask_profile) as d:
        d.write(np.where(coverage, flags, 255).astype("uint8"), 1)
        d.update_tags(bit0="some selected band has QUAC nodata/bad band", bit1="some reflectance >1",
                      bit2="some reflectance >2",
                      interpretation="screening flags, not cloud/snow/shadow labels; >1 preserved")
    with rasterio.open(out_dir / "all_selected_bands_available.tif", "w", **mask_profile) as d:
        d.write(np.where(coverage, available, 255).astype("uint8"), 1)

    result = {"scene_id": scene.scene_id, "reflectance_bands": len(rows), "quac_good_bands": int(sum(bbl)),
              "integer_scale_factor": scale, "nodata": nodata,
              "reflectance_file": str(path.resolve()),
              "coverage_area_km2": float(coverage.sum() * float(cfg["grid"]["resolution_m"]) ** 2 / 1e6),
              "all_bands_available_area_km2": float(available.sum() * float(cfg["grid"]["resolution_m"]) ** 2 / 1e6),
              "pixels_with_any_reflectance_above_1_fraction": float(np.mean((flags[coverage] & 2) > 0)),
              "pixels_with_any_reflectance_above_2_fraction": float(np.mean((flags[coverage] & 4) > 0)),
              "warning": "QUAC scene-based estimated reflectance; no independent accuracy assessment"}
    save_json(out_dir / "reflectance_quality_summary.json", result)
    return result
