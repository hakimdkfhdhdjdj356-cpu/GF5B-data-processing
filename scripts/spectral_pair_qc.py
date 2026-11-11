"""同点光谱质检：校正前辐亮度 vs 校正后光谱产品，按同一位置与空间支持配对。

纯 Python（numpy / rasterio / pandas / matplotlib / scipy），不依赖 ENVI 或任何商业软件。

设计约束（与 ahsi_pipeline 的产品等级纪律一致）：
- 每个产品必须自带 ENVI wavelength / fwhm 标签；波长轴按各自真实标签读取，绝不按波段序号硬对齐。
- 辐亮度 330 波段与校正后 254 波段按「最近波长 + 显式记录 delta_nm」配对，
  delta 超过 FWHM/2 时标记 wavelength_mismatch，不伪造同波长。
- 单位、产品等级、是否绝对反射率全部从 GeoTIFF tag 读取并写入产物；
  相对反射率绝不称为地表反射率，也不与绝对反射率做数值等同。
- 交付件只有中心波长与 FWHM，没有完整 SRF。脚本把高斯近似、其假设与
  对 SWIR 吸收特征的敏感性检查全部显式记录，不把近似说成真实响应。

用法：
  python scripts/spectral_pair_qc.py \
      --radiance <330-band radiance.tif> \
      --product iarr=<relative_reflectance.tif> \
      --product quac=<surface_reflectance.tif> \
      --points <representative_spectra.csv with site,row,col> \
      --out <output dir>
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from scipy.ndimage import gaussian_filter1d

SRF_NOTE = (
    "Product delivery provides per-band center wavelength and FWHM only; "
    "no full spectral response function (SRF) and no detector-level wavelength drift. "
    "The Gaussian smoothing is applied to interpolated, already sampled spectra. "
    "It probes sensitivity to additional smoothing, rather than modeling the native SRF."
)

WINDOWS = {
    "AlOH_kaolinite_2130_2240": (2130.0, 2240.0),
    "carbonate_FeMgOH_2290_2400": (2290.0, 2400.0),
}


def _envi_list(tags: dict, key: str, n: int) -> np.ndarray:
    raw = tags.get(key) or tags.get(key.capitalize()) or tags.get(key.upper()) or ""
    parts = [p for p in raw.strip("{} \t").replace(",", " ").split() if p]
    if len(parts) != n:
        raise ValueError(f"ENVI {key} 数量 {len(parts)} != 波段数 {n}，拒绝按序号推断波长")
    return np.asarray([float(p) for p in parts], float)


def _unit_and_grade(tags: dict, path: Path) -> tuple[str, str, str]:
    product = tags.get("product", "") or ""
    atm = tags.get("atmospheric_method", "") or ""
    absolute = tags.get("absolute_reflectance", "") or ""
    low = product.lower()
    if "radiance" in low:
        return (tags.get("units") or "UNSPECIFIED", "at_sensor_radiance", "NO")
    if "iarr" in atm.lower() or "relative reflectance" in low:
        return ("dimensionless (per-band scene-mean normalized)",
                "relative_reflectance_NOT_surface_reflectance", "NO")
    if "sixs_lut" in atm.lower():
        return ("dimensionless 6SV model-estimated surface reflectance",
                "sixs_model_estimated_surface_reflectance_NOT_independently_validated",
                absolute or "MODEL_ESTIMATE")
    if "reflectance" in low:
        return ("dimensionless reflectance estimate (nominal 0-1)",
                "estimated_reflectance_NOT_independently_validated",
                absolute or "NOT_STATED")
    return (f"UNSPECIFIED (file={path.name})", "UNSPECIFIED", "NOT_STATED")


def load_cube(path: Path) -> dict:
    with rasterio.open(path) as d:
        tags = d.tags()
        envi = d.tags(ns="ENVI")
        n = d.count
        wl = _envi_list(envi, "wavelength", n)
        fwhm = _envi_list(envi, "fwhm", n)
        unit, grade, absolute = _unit_and_grade(tags, path)
        if not (np.isfinite(wl).all() and np.all(np.diff(wl) >= 0)
                and (np.isfinite(fwhm) & (fwhm > 0)).all()):
            raise ValueError(f"{path.name}: 波长顺序或 FWHM 无效")
        origins = []
        for i in d.indexes:
            band = d.tags(i)
            camera = band.get("source_camera") or band.get("source_segment")
            source_band = band.get("source_native_band") or band.get("source_band")
            origins.append((camera, int(source_band)) if camera and source_band else None)
        return dict(origins=origins, path=str(path.resolve()), name=path.name, count=n,
                    width=d.width, height=d.height, crs=str(d.crs),
                    transform=list(d.transform)[:6], nodata=d.nodata,
                    wavelength_nm=wl, fwhm_nm=fwhm, unit=unit, grade=grade,
                    absolute_reflectance=absolute, product_tag=tags.get("product", ""),
                    atmospheric_method=tags.get("atmospheric_method", ""))


def read_spectrum(cube: dict, row: int, col: int) -> np.ndarray:
    with rasterio.open(cube["path"]) as d:
        v = d.read(window=((row, row + 1), (col, col + 1)),
                   masked=True)[:, 0, 0].astype(float).filled(np.nan)
    v[~np.isfinite(v)] = np.nan
    return v


def load_sites(points_csv: Path | None, ref: dict, n_default: int = 5) -> list[dict]:
    if points_csv:
        df = pd.read_csv(points_csv, encoding="utf-8-sig")
        need = {"site", "row", "col"}
        if not need.issubset(df.columns):
            raise ValueError(f"点位文件缺列 {need - set(df.columns)}")
        sites = (df[["site", "row", "col"]].drop_duplicates()
                 .sort_values("site").to_dict("records"))
    else:
        # 确定性空间分位选点，不引入随机性
        with rasterio.open(ref["path"]) as d:
            a = d.read(1, masked=True).astype(float).filled(np.nan)
        ys, xs = np.nonzero(np.isfinite(a))
        if not len(ys):
            raise ValueError("辐亮度第一波段没有有效选点")
        order = np.argsort(ys * ref["width"] + xs)
        picks = np.linspace(0, len(order) - 1, min(n_default, len(order))).astype(int)
        sites = [{"site": f"Q{i + 1}", "row": int(ys[order[k]]), "col": int(xs[order[k]])}
                 for i, k in enumerate(picks)]
    if not sites or len({s["site"] for s in sites}) != len(sites) or any(pd.isna(s["site"]) for s in sites):
        raise ValueError("点位表须包含非空且唯一的点名")
    for s in sites:
        if any(not np.isfinite(float(s[k])) or not float(s[k]).is_integer() for k in ("row", "col")):
            raise ValueError(f"行列号须为整数: {s}")
        s["row"], s["col"] = int(s["row"]), int(s["col"])
        if not (0 <= s["row"] < ref["height"] and 0 <= s["col"] < ref["width"]):
            raise ValueError(f"点位越界: {s}")
    return sites


def _upper_hull(w: np.ndarray, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """从左到右的上凸包：所有样本点位于该折线之下，作为连续统。"""
    pts = sorted(zip(w.tolist(), x.tolist()))
    up: list[tuple[float, float]] = []
    for p in reversed(pts):
        while len(up) >= 2:
            (ox, oy), (ax, ay) = up[-2], up[-1]
            cross = (ax - ox) * (p[1] - oy) - (ay - oy) * (p[0] - ox)
            if cross <= 0:
                up.pop()
            else:
                break
        up.append(p)
    up.reverse()
    return (np.asarray([p[0] for p in up], float),
            np.asarray([p[1] for p in up], float))


def continuum_removal(wl: np.ndarray, v: np.ndarray, lo: float, hi: float,
                      min_depth: float = 0.01) -> dict:
    """窗口内上凸包连续统去除；返回吸收中心、深度与有效波段数。

    v1（已废弃）用窗口两端点连线，对单调下降谱会把右端点误判为吸收中心；
    v2 改用上凸包，并显式给出 absorption_detected，无显著吸收时不报中心。
    """
    m = (wl >= lo) & (wl <= hi) & np.isfinite(v) & (v > 0)
    n = int(m.sum())
    if n < 5:
        return {"valid": False, "n_bands": n, "reason": "fewer than 5 valid bands"}
    w, x = wl[m], v[m]
    hw, hx = _upper_hull(w, x)
    if len(hw) < 2 or np.any(hx <= 0):
        return {"valid": False, "n_bands": n, "reason": "invalid upper hull"}
    cont = np.interp(w, hw, hx)
    if np.any(cont <= 0):
        return {"valid": False, "n_bands": n, "reason": "non-positive continuum"}
    depth = 1.0 - x / cont
    j = int(np.argmax(depth))
    detected = bool(depth[j] >= min_depth)
    return {"valid": True, "n_bands": n,
            "window_nm": [lo, hi],
            "continuum_method": "upper_convex_hull_v2",
            "hull_points": int(len(hw)),
            "absorption_detected": detected,
            "absorption_center_nm": round(float(w[j]), 3) if detected else None,
            "max_depth": round(float(depth[j]), 6),
            "min_depth_for_detection": min_depth,
            "continuum_left_nm": round(float(hw[0]), 3),
            "continuum_right_nm": round(float(hw[-1]), 3),
            "median_fwhm_nm": None}


def srf_sensitivity(wl: np.ndarray, v: np.ndarray, fwhm: np.ndarray,
                    lo: float, hi: float) -> dict:
    """假设高斯 SRF：比较标称 FWHM 与 2×FWHM 展宽后的吸收中心漂移。"""
    m = (wl >= lo - 60) & (wl <= hi + 60) & np.isfinite(v) & (v > 0)
    if int(m.sum()) < 5:
        return {"valid": False, "reason": "fewer than 5 valid bands"}
    w, x, f = wl[m], v[m], fwhm[m]
    sigma = float(np.median(f)) / 2.3548
    if sigma <= 0:
        return {"valid": False, "reason": "non-positive FWHM"}
    step = 1.0
    grid = np.arange(math.floor(w.min()) - 10, math.ceil(w.max()) + 11, step)
    spec = np.interp(grid, w, x)
    if w.min() > lo or w.max() < hi:
        return {"valid": False, "reason": "observed wavelengths do not bracket the window"}
    out = {}
    for label, factor in (("nominal_fwhm", 1.0), ("doubled_fwhm", 2.0)):
        sm = gaussian_filter1d(spec, max(sigma * factor / step, 1e-6))
        r = continuum_removal(grid, sm, lo, hi)
        out[label] = r
    a, b = out["nominal_fwhm"], out["doubled_fwhm"]
    out["depth_change"] = (round(b["max_depth"] - a["max_depth"], 6)
                           if a.get("valid") and b.get("valid") else None)
    ca, cb = a.get("absorption_center_nm"), b.get("absorption_center_nm")
    out["center_shift_nm"] = (round(cb - ca, 3)
                              if ca is not None and cb is not None else None)
    out["assumed_gaussian_sigma_nm"] = round(sigma, 4)
    out["valid"] = True
    return out


def flags_for(grade: str, v: float) -> list[str]:
    out = []
    if not np.isfinite(v):
        return ["missing"]
    if v <= 0:
        out.append("non_positive")
    if "relative" in grade and v > 2:
        out.append("above_2_relative")
    if "relative" in grade and v > 3:
        out.append("above_3_relative")
    if "estimated" in grade and v > 1:
        out.append("above_1_reflectance")
    if "estimated" in grade and v > 2:
        out.append("above_2_reflectance")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="同点光谱质检（纯 Python）")
    ap.add_argument("--radiance", required=True, help="校正前定标辐亮度立方体")
    ap.add_argument("--product", action="append", default=[], metavar="LABEL=PATH",
                    help="校正后产品，可重复；标签用于输出命名")
    ap.add_argument("--points", help="同点 CSV，需含 site,row,col；缺省按空间分位确定性选 5 点")
    ap.add_argument("--out", required=True, help="输出目录")
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"请指定新的质检输出目录: {out}")
    rad = load_cube(Path(args.radiance))
    products = {}
    for spec in args.product:
        if "=" not in spec:
            raise SystemExit(f"--product 需为 LABEL=PATH: {spec}")
        label, path = spec.split("=", 1)
        if not label or label == "radiance" or label in products:
            raise ValueError(f"产品标签为空、重复或占用保留名 radiance: {label!r}")
        products[label] = load_cube(Path(path))
    if not products:
        raise SystemExit("至少需要一个 --product")

    # 同一空间支持：要求所有立方体网格完全一致，不做隐式重采样
    grid_keys = ("width", "height", "crs", "transform")
    grid_match = rad["crs"] != "None" and all(
        all(p[k] == rad[k] for k in grid_keys) for p in products.values())
    if not grid_match:
        raise SystemExit("输入立方体网格不一致；请先配准到同一格网，禁止隐式重采样伪造同点")

    sites = load_sites(Path(args.points) if args.points else None, rad)
    cubes = {"radiance": rad, **products}

    spectra = {(label, s["site"]): read_spectrum(c, s["row"], s["col"])
               for label, c in cubes.items() for s in sites}

    rows, stats = [], {label: {"observations": 0, "missing": 0,
                               "missing_radiance": 0, "missing_product": 0, "non_positive": 0,
                               "above_1_reflectance": 0, "above_2_reflectance": 0,
                               "above_2_relative": 0, "above_3_relative": 0,
                               "wavelength_mismatch": 0}
                       for label in products}
    swir = {label: {} for label in cubes}

    source_lookup = {origin: i for i, origin in enumerate(rad["origins"]) if origin is not None}
    pairs = {}
    for label, cube in products.items():
        pairs[label] = [(source_lookup[origin], "camera_source_band") if origin in source_lookup
                        else (int(np.argmin(np.abs(rad["wavelength_nm"] - wl))), "nearest_wavelength")
                        for origin, wl in zip(cube["origins"], cube["wavelength_nm"])]
    for s in sites:
        site = s["site"]
        rv = spectra[("radiance", site)]
        for label, c in products.items():
            pv = spectra[(label, site)]
            for j, wl in enumerate(c["wavelength_nm"]):
                k, pair_method = pairs[label][j]
                dnm = float(abs(wl - rad["wavelength_nm"][k]))
                fwhm = float(c["fwhm_nm"][j])
                mismatch = bool(dnm > max(fwhm, 1.0) / 2.0)
                r, p = float(rv[k]), float(pv[j])
                fl = flags_for(c["grade"], p)
                st = stats[label]
                st["observations"] += 1
                missing_r, missing_p = not np.isfinite(r), not np.isfinite(p)
                st["missing_radiance"] += int(missing_r)
                st["missing_product"] += int(missing_p)
                st["missing"] += int(missing_r or missing_p)
                if missing_r:
                    fl.append("missing_radiance")
                for f in ("non_positive", "above_1_reflectance", "above_2_reflectance",
                          "above_2_relative", "above_3_relative"):
                    st[f] += int(f in fl)
                st["wavelength_mismatch"] += int(mismatch)
                rows.append(dict(site=site, row=s["row"], col=s["col"], product=label,
                                 product_wavelength_nm=round(float(wl), 4), fwhm_nm=fwhm,
                                 radiance_wavelength_nm=round(float(rad["wavelength_nm"][k]), 4),
                                 radiance_band_index=k + 1, delta_wavelength_nm=round(dnm, 4),
                                 wavelength_mismatch=mismatch, pair_method=pair_method,
                                 radiance_value=None if not np.isfinite(r) else r,
                                 radiance_unit=rad["unit"],
                                 product_value=None if not np.isfinite(p) else p,
                                 product_unit=c["unit"], product_grade=c["grade"],
                                 flags="|".join(fl + (["wavelength_mismatch"] if mismatch else []))))
        for label, c in cubes.items():
            v = spectra[(label, site)]
            swir[label][site] = {}
            for wname, (lo, hi) in WINDOWS.items():
                r = continuum_removal(c["wavelength_nm"], v, lo, hi)
                r["median_fwhm_nm"] = (round(float(np.median(c["fwhm_nm"][
                    (c["wavelength_nm"] >= lo) & (c["wavelength_nm"] <= hi)])), 4)
                    if np.any((c["wavelength_nm"] >= lo) & (c["wavelength_nm"] <= hi)) else None)
                r["srf_sensitivity"] = srf_sensitivity(c["wavelength_nm"], v,
                                                       c["fwhm_nm"], lo, hi)
                swir[label][site][wname] = r

    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(out / "paired_spectra.csv", index=False, encoding="utf-8-sig")

    report = {
        "purpose": "same-pixel, same-spatial-support paired spectral QA; NO accuracy claim",
        "srf_note": SRF_NOTE,
        "grid_match_all_products": grid_match,
        "grid": {k: rad[k] for k in ("width", "height", "crs", "resolution_m_note", "nodata")
                 if k in rad} | {"transform": rad["transform"]},
        "sites": sites,
        "inputs": {label: {k: c[k] for k in ("path", "count", "unit", "grade",
                                             "absolute_reflectance", "product_tag",
                                             "atmospheric_method")}
                   for label, c in cubes.items()},
        "wavelength_ranges_nm": {label: [float(c["wavelength_nm"][0]),
                                         float(c["wavelength_nm"][-1])]
                                 for label, c in cubes.items()},
        "per_product_stats": stats,
        "swir_absorption": swir,
        "continuum_removal_version": "upper_convex_hull_v2 (v1 endpoint-line retired; see pair_qc_endpoint_v1_SUPERSEDED.json)",
        "caveats": [
            "Radiance has units W m^-2 sr^-1 um^-1; relative and estimated reflectance are dimensionless and NOT interchangeable.",
            "v1 continuum used window endpoint line and mis-reported the right window edge as the absorption center on monotonically decreasing spectra; v2 uses an upper convex hull and reports absorption_detected=False when max depth < 0.01.",
            "IARR per-band normalization changes the spectral vector direction; its absorption depth is shape-only, not absolute reflectance.",
            "Continuum removal on at-sensor radiance is dominated by the solar irradiance envelope and atmospheric transmission; compare shapes only.",
            "Gaussian SRF is an assumption. Delivery has no full SRF, so center/depth shifts above are sensitivity diagnostics, not sensor calibration.",
            "Sites are single pixels; no spatial support aggregation, no cloud/shadow screening, no independent ground truth.",
        ],
    }
    (out / "pair_qc.json").write_text(json.dumps(report, ensure_ascii=False, indent=2),
                                      encoding="utf-8")

    # 图：每个产品一列，上排全谱、下排 SWIR；所有点同图，英文标签避免字体缺失
    labels = ["radiance"] + list(products)
    fig, axes = plt.subplots(2, len(labels), figsize=(5.2 * len(labels), 7.4),
                             squeeze=False)
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(sites), 1)))
    for ci, label in enumerate(labels):
        c = cubes[label]
        for s, color in zip(sites, colors):
            v = spectra[(label, s["site"])]
            for ri, (lo, hi) in enumerate([(380, 2500), (2100, 2450)]):
                m = (c["wavelength_nm"] >= lo) & (c["wavelength_nm"] <= hi)
                axes[ri][ci].plot(c["wavelength_nm"][m], v[m], lw=0.9,
                                  color=color, label=s["site"] if ri == 0 else None)
        axes[0][ci].set_title(f"{label}\n{c['unit']}", fontsize=9)
        axes[1][ci].set_title(f"{label} SWIR", fontsize=9)
        for ri in (0, 1):
            axes[ri][ci].set_xlabel("wavelength (nm)", fontsize=8)
            axes[ri][ci].grid(alpha=0.3, lw=0.4)
            axes[ri][ci].tick_params(labelsize=7)
        axes[0][ci].set_ylabel("at-sensor radiance" if label == "radiance"
                               else "dimensionless", fontsize=8)
        axes[1][ci].set_ylabel("value in product units", fontsize=8)
    axes[0][0].legend(fontsize=7, ncol=2)
    fig.suptitle("Same-pixel paired spectra (before / after correction) - units differ, do NOT compare magnitudes",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out / "paired_spectra.png", dpi=140)
    plt.close(fig)

    print(json.dumps({
        "out": str(out.resolve()), "sites": len(sites),
        "products": list(products), "paired_rows": len(df),
        "grid_match": grid_match,
        "per_product_stats": stats,
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
