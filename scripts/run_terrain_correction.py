"""GF-5B 反射率产品的地形光照辐射校正阶段（可独立重跑，零商业软件依赖）。

阶段定位：接在 6SV 大气校正之后、下游光谱分类之前。不修改 run_ahsi.py 的
阶段调度，因此不会破坏既有 run 的断点续跑标记与冻结产物；用独立输出目录。

示例：
  python scripts/run_terrain_correction.py \
      --reflectance <.../GF5B_..._6SV_surface_reflectance.tif> \
      --dem <assets/rpc_dem_ellipsoidal.tif> \
      --xml <原始产品 .xml> \
      --coverage <.../joint_camera_coverage_mask.tif> \
      --out <.../terrain_corrected_final>
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import rasterio

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ahsi import terrain  # noqa: E402


def read_sun(xml_path: Path) -> terrain.SunGeometry:
    root = ET.parse(xml_path).getroot()
    flat = {e.tag: (e.text or "").strip() for e in root.iter()}
    aliases = {"zenith": ("SolarZenith", "SunZenith", "SolarZenithAngle"),
               "azimuth": ("SolarAzimuth", "SunAzimuth", "SolarAzimuthAngle")}
    def pick(keys):
        for k in keys:
            if k in flat and flat[k]:
                return float(flat[k])
        raise ValueError(f"XML 缺少太阳角字段 {keys}: {xml_path}")
    return terrain.SunGeometry(pick(aliases["zenith"]),
                               pick(aliases["azimuth"]),
                               source=f"xml_center:{xml_path.name}")


def read_wavelengths(ds) -> tuple[np.ndarray, np.ndarray, dict]:
    envi = ds.tags(ns="ENVI")
    wl = envi.get("wavelength", "").strip("{} \n")
    fw = envi.get("fwhm", "").strip("{} \n")
    if not wl:
        raise ValueError("反射率立方体缺 ENVI wavelength 标签，拒绝按假设波段号校正")
    w = np.asarray([float(v) for v in wl.replace(",", " ").split()], float)
    f = (np.asarray([float(v) for v in fw.replace(",", " ").split()], float)
         if fw else np.full_like(w, np.nan))
    if (len(w) != ds.count or not np.isfinite(w).all() or not np.all(np.diff(w) > 0)):
        raise ValueError("波长标签须有限、严格递增且与波段数一致")
    if fw and (f.shape != w.shape or not (np.isfinite(f) & (f > 0)).all()):
        raise ValueError("FWHM 标签与波段数不一致或含无效值")
    return w, f, envi


def write_tif(path: Path, arr: np.ndarray, profile: dict, tags: dict,
              envi_tags: dict | None = None, descriptions=None, band_tags=None):
    p = profile.copy()
    p.update(count=arr.shape[0], dtype="float32",
             compress="deflate", tiled=True, blockxsize=256, blockysize=256)
    with rasterio.open(path, "w", **p) as d:
        d.write(np.asarray(arr, "float32"))
        d.update_tags(**{k: str(v) for k, v in tags.items()})
        if envi_tags:
            d.update_tags(ns="ENVI", **envi_tags)
        if descriptions:
            for i, t in enumerate(descriptions, start=1):
                d.set_band_description(i, str(t))
        for i, metadata in enumerate(band_tags or (), start=1):
            d.update_tags(i, **metadata)
    return path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reflectance", required=True, type=Path)
    ap.add_argument("--dem", required=True, type=Path)
    ap.add_argument("--xml", required=True, type=Path)
    ap.add_argument("--coverage", type=Path, default=None)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--method", choices=("scs_c", "c_correction", "cosine"),
                    default="scs_c")
    ap.add_argument("--min-cos-i", type=float, default=0.10)
    ap.add_argument("--min-regression-pixels", type=int, default=1000)
    ap.add_argument("--force", action="store_true", help="允许覆盖已有输出（默认拒绝）")
    a = ap.parse_args()

    out = a.out
    stem = a.reflectance.stem.removesuffix("_surface_reflectance")
    main_tif = out / f"{stem}_terrain_corrected_reflectance.tif"
    if out.exists() and any(out.iterdir()) and not a.force:
        print(f"拒绝覆盖已有输出目录: {out}（如确需重算请显式 --force）", file=sys.stderr)
        return 1
    sun = read_sun(a.xml)
    sun.validate()
    with rasterio.open(a.reflectance) as ds:
        profile, tags = ds.profile.copy(), dict(ds.tags())
        wl, fwhm, envi = read_wavelengths(ds)
        band_tags = [ds.tags(i) for i in ds.indexes]
        rho = ds.read(masked=True).astype("float32").filled(np.nan)
        nodata = ds.nodata if ds.nodata is not None else -9999.0
        crs, transform = ds.crs, ds.transform
        h, w = ds.height, ds.width
        grade = tags.get("product", "unknown")
        prev_atm = tags.get("atmospheric_correction", "unknown")
        if (tags.get("atmospheric_method") != "sixs_lut"
                or prev_atm != "APPLIED_6SV_RADIATIVE_TRANSFER"):
            raise ValueError("此入口需要带 6SV 大气校正标记的反射率产品")
        if tags.get("terrain_illumination_correction", "NOT_APPLIED") != "NOT_APPLIED":
            print("输入产品已做过地形校正，拒绝二次校正", file=sys.stderr)
            return 1
    print(f"输入产品: {grade}")
    print(f"大气校正: {prev_atm} | 波段 {rho.shape[0]} | 格网 {h}x{w} | {crs}")
    print(f"太阳几何: zenith={sun.zenith_deg} azimuth={sun.azimuth_deg} ({sun.source})")

    if (not crs or not crs.is_projected or crs.linear_units_factor[1] != 1.0
            or transform.a <= 0 or transform.e >= 0 or transform.b or transform.d
            or not np.isclose(transform.a, -transform.e)):
        raise ValueError("地形梯度要求北向上的米制方格投影")
    dem = terrain.reproject_dem(a.dem, crs, transform, w, h)
    print(f"DEM 重投影完成，有效像元 {int(np.isfinite(dem).sum())}/{dem.size}")

    cell_m = float(abs(transform.a))
    # 缺 DEM 邻点的梯度一并失效，避免填充值在 AOI 边界制造人工陡坡。
    slope, aspect = terrain.slope_aspect(dem, cell_m)
    slope[~np.isfinite(dem)] = np.nan
    aspect[~np.isfinite(dem)] = np.nan
    if not np.isfinite(slope).any():
        raise ValueError("DEM 与反射率格网没有可计算坡度的重叠区域")
    cos_i = terrain.illumination_cosine(slope, aspect, sun)

    coverage = None
    if a.coverage:
        with rasterio.open(a.coverage) as d:
            if (d.shape, d.crs, d.transform) != ((h, w), crs, transform):
                raise ValueError("覆盖掩膜与反射率格网不一致")
            coverage = d.read(1, masked=True).filled(0) == 1
        print(f"共同覆盖掩膜有效像元 {int(coverage.sum())}")

    result = terrain.apply_correction(
        rho, cos_i, slope, dem, sun, wl, method=a.method,
        min_cos_i=a.min_cos_i,
        min_regression_pixels=a.min_regression_pixels,
        nodata=nodata, coverage=coverage, aspect=aspect)
    m = result.meta
    print(f"方法 {a.method} | C 回归有效波段 {m['c_fit_valid_bands']}/{rho.shape[0]}")
    print(f"坡度中位 {m['slope_deg']['median']:.2f}° | cos(i) 中位 {m['cos_i']['median']:.4f}")
    print(f"低光照像元比例 {m['low_illumination_fraction']:.4f}")
    print("校正相对变化 p05/med/p95 =", [m[f"relative_change_{q}"]
                                            for q in ("p05", "median", "p95")])
    out.mkdir(parents=True, exist_ok=True)
    profile["nodata"] = nodata

    new_tags = dict(tags)
    new_tags.update({
        "product": ("6SV model-estimated surface reflectance with first-order "
                    "terrain illumination correction; NOT independently validated"),
        "atmospheric_method": tags.get("atmospheric_method", "sixs_lut"),
        "terrain_illumination_correction": f"APPLIED_{a.method.upper()}",
        "terrain_correction_sun_geometry": f"zenith={sun.zenith_deg}; azimuth={sun.azimuth_deg}; {sun.source}",
        "terrain_correction_dem": str(a.dem.name),
        "terrain_pixelwise_ogp_angles": "NOT_APPLIED",
        "terrain_diffuse_adjacency_brdf": "NOT_APPLIED",
        "absolute_reflectance": "MODEL_ESTIMATE",
        "quantitative_accuracy_independently_validated": "NO",
        "usable_for": ("within-scene spectral shape, continuum removal, relative mineral "
                       "diagnostics under the SAME atmospheric+terrain assumptions"),
        "NOT_usable_for": ("absolute accuracy claims, direct lab-library SAM without "
                           "re-evaluation, cross-scene quantitative comparison"),
        "caveat_spectral_shape": ("per-band C changes spectral vector direction; "
                                  "re-evaluate before SAM/library matching"),
        "upstream_product": str(a.reflectance.resolve()),
    })
    write_tif(main_tif, result.corrected, profile, new_tags, envi,
              descriptions=[f"{x:.3f} nm" for x in wl], band_tags=band_tags)
    print(f"写出 {main_tif}")

    with rasterio.open(out / "terrain_qa.tif", "w", **{
            **profile, "count": 1, "dtype": "uint8", "nodata": None,
            "compress": "deflate", "predictor": 2}) as d:
        d.write(result.qa[np.newaxis])
        d.update_tags(**{f"qa_bit_{k}": str(v) for k, v in terrain.QA_BITS.items()},
                      method=a.method)

    diag = out / "terrain_geometry_diagnostics.tif"
    with rasterio.open(diag, "w", **{**profile, "count": 4, "dtype": "float32",
                                     "nodata": -9999.0, "compress": "deflate"}) as d:
        for i, (arr, name) in enumerate(zip(
                [result.dem_on_grid, result.slope_deg, result.aspect_deg, result.cos_i],
                ["ellipsoidal_height_m", "slope_deg", "aspect_deg_from_north_cw",
                 "cos_incidence_angle"]), start=1):
            d.write(np.where(np.isfinite(arr), arr, -9999.0).astype("float32"), i)
            d.set_band_description(i, name)

    with (out / "terrain_band_parameters.csv").open("w", encoding="utf-8-sig",
                                                    newline="") as f:
        cols = ["band", "wavelength_nm", "c_fit_valid", "c_parameter", "c_source",
                "slope_m", "intercept_b", "r_squared", "n_regression_pixels", "reason"]
        wr = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        wr.writeheader()
        for r in result.band_records:
            wr.writerow(r)

    meta = dict(m)
    meta.update({
        "stage": "terrain_illumination_correction",
        "input_reflectance": str(a.reflectance.resolve()),
        "input_product_grade": grade,
        "dem_source": str(a.dem.resolve()),
        "coverage_mask": str(a.coverage.resolve()) if a.coverage else None,
        "grid": {"crs": str(crs), "transform": list(transform)[:6],
                 "shape": [h, w], "resolution_m": cell_m},
        "wavelength_nm_first_last": [float(wl[0]), float(wl[-1])],
        "qa_bits": terrain.QA_BITS,
        "output_files": {"reflectance": main_tif.name, "qa": "terrain_qa.tif",
                         "geometry": diag.name,
                         "band_parameters": "terrain_band_parameters.csv"},
    })
    (out / "terrain_correction_meta.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    print(f"元数据与 QA 已写出到 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
