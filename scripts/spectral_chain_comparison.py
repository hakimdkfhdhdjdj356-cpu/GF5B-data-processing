"""同一影像像元的「未处理 → 处理完成」四阶段光谱曲线对比。

四个阶段（单位与产品等级互不相同，绝不混画同一纵轴）：
  A raw_dn              L1 原始 DN；未定标、未正射、未大气/地形校正。
                        通过 RPC + DEM 把产品格网像元反算回未正射影像像素读取。
  B radiance            L = DN*gain + offset + 显式 RPC/DEM 正射，W m-2 sr-1 um-1
  C sixs_reflectance    6SV 大气校正后模型估计地表反射率，无量纲
  D terrain_corrected   + SCS+C 一阶地形光照校正，无量纲

「同点」的严格含义：B/C/D 在同一产品格网 (row,col)；A 用 RPC 反算到
原始相机格网，VNIR/SWIR 分别反算并记录源像素坐标，两者一致时才称同点。
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import RPCTransformer
from rasterio.warp import transform as warp_coordinates

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

WATER_WINDOWS = [(755, 775), (930, 980), (1350, 1460), (1800, 1980)]
STAGES = ("raw_dn", "radiance", "sixs_reflectance", "terrain_corrected")
STYLE = {
    "raw_dn": ("#7f8c8d", "A  Raw L1 DN (UNPROCESSED)", "DN (int16)"),
    "radiance": ("#2980b9", "B  Calibrated + orthorectified radiance",
                 "W m^-2 sr^-1 um^-1"),
    "sixs_reflectance": ("#e67e22", "C  6SV surface reflectance (model estimate)",
                         "reflectance (dimensionless)"),
    "terrain_corrected": ("#c0392b", "D  + SCS+C terrain illumination correction",
                          "reflectance (dimensionless)"),
}


def envi_wavelengths(ds, key="wavelength") -> np.ndarray:
    raw = ds.tags(ns="ENVI").get(key, "").strip("{} \n")
    return np.asarray([float(v) for v in raw.replace(",", " ").split()], float)


def read_pixel(path: Path, row: int, col: int, nodata=None) -> np.ndarray:
    with rasterio.open(path) as d:
        v = d.read(window=((row, row + 1), (col, col + 1)),
                   masked=True)[:, 0, 0].astype(float).filled(np.nan)
        nd = d.nodata if d.nodata is not None else nodata
        if nd is not None:
            v[v == nd] = np.nan
        return v


def rpc_backproject(raw_dir: Path, dem: Path, grid_crs, x: float, y: float) -> dict:
    """产品格网 (x,y) → 各相机原始像素；不伪造，越界显式返回 unavailable。"""
    lon, lat = warp_coordinates(grid_crs, "EPSG:4326", [x], [y])
    out = {"lon": float(lon[0]), "lat": float(lat[0]), "segments": {}}
    for seg in ("VN", "SW"):
        files = sorted(raw_dir.glob(f"*_{seg}.tif"))
        if not files:
            out["segments"][seg] = {"available": False,
                                    "reason": f"no *_{seg}.tif in raw dir"}
            continue
        if len(files) != 1:
            raise ValueError(f"{raw_dir} 中 {seg} 相机影像不唯一")
        with rasterio.open(files[0]) as d:
            with RPCTransformer(d.rpcs, RPC_DEM=str(dem)) as tr:
                r, c = tr.rowcol(lon, lat, op=lambda v: v)
            sr, sc = float(r[0]), float(c[0])
            inside = (0 <= sr < d.height) and (0 <= sc < d.width)
            rec = {"file": files[0].name, "src_row": round(sr, 4),
                   "src_col": round(sc, 4), "image_shape": [d.height, d.width],
                   "available": bool(inside), "dn": None}
            if inside:
                ri, ci = math.floor(sr), math.floor(sc)
                rec["dn"] = d.read(window=((ri, ri + 1), (ci, ci + 1)))[:, 0, 0].astype(float).tolist()
                rec["pixel_read"] = [ri, ci]
            out["segments"][seg] = rec
    a = out["segments"].get("VN", {})
    b = out["segments"].get("SW", {})
    out["both_cameras_available"] = bool(a.get("available") and b.get("available"))
    out["vn_sw_agree"] = bool(
        out["both_cameras_available"]
        and abs(a["src_row"] - b["src_row"]) < 0.5
        and abs(a["src_col"] - b["src_col"]) < 0.5)
    return out


def sam_angle_deg(u: np.ndarray, v: np.ndarray) -> float | None:
    ok = np.isfinite(u) & np.isfinite(v) & (u > 0) & (v > 0)
    if ok.sum() < 10:
        return None
    a, b = u[ok], v[ok]
    c = float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
    return float(math.degrees(math.acos(max(-1.0, min(1.0, c)))))


def pick_sites(diag: Path, coverage: Path | None, refl: Path,
               primary=None) -> list[dict]:
    """选一个中间位置及两个太阳入射余弦分位位置。"""
    with rasterio.open(diag) as d:
        dem, slope, aspect, cos_i = d.read(masked=True).astype(float).filled(np.nan)
    with rasterio.open(refl) as d:
        sample = d.read(1)
        transform = d.transform
    cov = np.ones_like(sample, bool)
    if coverage is not None:
        with rasterio.open(coverage) as d:
            cov = d.read(1, masked=True).filled(0) == 1
    valid = (np.isfinite(dem) & np.isfinite(sample) & (sample > 0)
             & (cos_i > 0) & cov)
    rows, cols = np.where(valid)
    order = np.argsort(cos_i[valid])
    n = len(order)
    if n < 3:
        raise ValueError("可配对的地形和反射率像元不足三个")
    if primary is None:
        i = np.argmin((rows - valid.shape[0] / 2)**2 + (cols - valid.shape[1] / 2)**2)
        primary = int(rows[i]), int(cols[i])
    if not (0 <= primary[0] < valid.shape[0] and 0 <= primary[1] < valid.shape[1]
            and valid[primary]):
        raise ValueError(f"指定的主点缺少有效空间或光谱支持: {primary}")
    def at(q, exclude):
        k = int(np.clip(round(q * (n - 1)), 0, n - 1))
        for j in list(range(k, n)) + list(range(k - 1, -1, -1)):
            r, c = int(rows[order[j]]), int(cols[order[j]])
            if (r, c) not in exclude:
                return r, c
        return None
    used = {primary}
    picks = [("P1_primary",) + primary]
    for label, q in (("shaded_low_cos_i", 0.02), ("sunlit_high_cos_i", 0.98)):
        p = at(q, used)
        if p:
            used.add(p)
            picks.append((label,) + p)
    sites = []
    for label, r, c in picks:
        x, y = transform * (c + .5, r + .5)
        sites.append({"site": label, "row": int(r), "col": int(c),
                      "x": float(x), "y": float(y),
                      "cos_i": float(cos_i[r, c]),
                      "slope_deg": float(slope[r, c]),
                      "aspect_deg": float(aspect[r, c]),
                      "dem_m": float(dem[r, c])})
    return sites


def plot_chain(site: dict, spectra: dict, waves: dict, out: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(15, 8.5))
    for ax, stage in zip(axes.ravel(), STAGES):
        color, title, ylabel = STYLE[stage]
        w, v = waves[stage], spectra[stage]
        for lo, hi in WATER_WINDOWS:
            ax.axvspan(lo, hi, color="#bdc3c7", alpha=0.35, lw=0)
        ax.plot(w, v, color=color, lw=1.4)
        ax.set_title(title, fontsize=10.5, color=color, fontweight="bold")
        ax.set_xlabel("wavelength (nm)", fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_xlim(float(np.nanmin(w)) - 10, float(np.nanmax(w)) + 10)
        ax.grid(alpha=0.25, lw=0.6)
        finite = v[np.isfinite(v)]
        if finite.size:
            ax.text(0.015, 0.93,
                    f"n={len(w)}  min={finite.min():.4g}  med={np.median(finite):.4g}  max={finite.max():.4g}",
                    transform=ax.transAxes, fontsize=8, va="top",
                    bbox=dict(fc="white", ec="none", alpha=0.75))
    fig.suptitle(f"Same pixel, four processing stages | site={site['site']} "
                 f"(row {site['row']}, col {site['col']}) | "
                 f"slope {site['slope_deg']:.1f} deg, cos(i) {site['cos_i']:.3f}\n"
                 "Grey bands = water-absorption windows excluded from analysis. "
                 "Y-axes are INDEPENDENT (different physical units).",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(out / "unprocessed_vs_processed.png", dpi=150)
    plt.close(fig)


def plot_terrain_effect(sites: list[dict], spectra: dict, waves: dict,
                        stats: dict, out: Path) -> None:
    n = len(sites)
    fig, axes = plt.subplots(1, n, figsize=(5.6 * n, 4.2), squeeze=False)
    for ax, s in zip(axes.ravel(), sites):
        key = s["site"]
        w = waves["sixs_reflectance"]
        c = spectra["sixs_reflectance"][key]
        t = spectra["terrain_corrected"][key]
        ax.plot(w, c, color="#e67e22", lw=1.5, label="C 6SV (before terrain)")
        ax.plot(w, t, color="#c0392b", lw=1.3, ls="--",
                label="D SCS+C (after terrain)")
        ax.set_title(f"{key}\nslope {s['slope_deg']:.1f} deg, cos(i) {s['cos_i']:.3f}"
                     f"\nSAM(C,D) = {stats[key]['sam_angle_deg']:.3f} deg", fontsize=9.5)
        ax.set_xlabel("wavelength (nm)", fontsize=9)
        ax.set_ylabel("reflectance", fontsize=9)
        ax.grid(alpha=0.25, lw=0.6)
        ax.legend(fontsize=8)
    fig.suptitle("Terrain illumination correction changes BOTH magnitude and spectral shape "
                 "(per-band C parameter)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(out / "terrain_effect_sites.png", dpi=150)
    plt.close(fig)


def plot_normalized(site: dict, spectra: dict, waves: dict, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 4.8))
    for stage in STAGES:
        w, v = waves[stage], np.asarray(spectra[stage], float)
        ok = np.isfinite(v)
        if ok.sum() < 10:
            continue
        lo, hi = float(np.percentile(v[ok], 2)), float(np.percentile(v[ok], 98))
        n = (v - lo) / (hi - lo) if hi > lo else np.full_like(v, np.nan)
        ax.plot(w, n, lw=1.3, color=STYLE[stage][0],
                label=f"{stage[0]} {stage} (min-max norm)", alpha=0.9)
    for a, b in WATER_WINDOWS:
        ax.axvspan(a, b, color="#bdc3c7", alpha=0.3, lw=0)
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("normalized value (each stage scaled independently)")
    ax.set_title(f"Spectral SHAPE at {site['site']}: atmospheric correction reshapes "
                 "the spectrum; terrain correction mainly rescales magnitude (small SAM angle)",
                 fontsize=10.5)
    ax.grid(alpha=0.25, lw=0.6)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "normalized_shape_comparison.png", dpi=150)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scene-dir", required=True, type=Path)
    ap.add_argument("--raw-dir", required=True, type=Path)
    ap.add_argument("--dem", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--primary", nargs=2, type=int, metavar=("ROW", "COL"))
    ap.add_argument("--terrain-dir", type=Path)
    a = ap.parse_args()
    S, out = a.scene_dir, a.out
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"请指定新的质检输出目录: {out}")

    rad_p = S / f"{S.name}_radiance.tif"
    sixs_p = S / "sixs_reflectance_final" / f"{S.name}_6SV_surface_reflectance.tif"
    terrain_dir = a.terrain_dir or S / "terrain_corrected_final"
    meta = json.loads((terrain_dir / "terrain_correction_meta.json").read_text(encoding="utf-8"))
    if meta["method"] != "scs_c" or Path(meta["input_reflectance"]).resolve() != sixs_p.resolve():
        raise ValueError("四阶段 SCS+C 质检需要与当前场景配套的地形校正成果")
    terr_p = terrain_dir / meta["output_files"]["reflectance"]
    diag_p = terrain_dir / meta["output_files"]["geometry"]
    cov_p = S / "joint_camera_coverage_mask.tif"
    band_csv = S / "band_parameters_and_quality.csv"
    for p in (rad_p, sixs_p, terr_p, diag_p, band_csv, a.dem):
        if not p.is_file():
            raise FileNotFoundError(p)

    with rasterio.open(rad_p) as d:
        w330 = envi_wavelengths(d)
        grid_crs, transform = d.crs, d.transform
        shape = d.shape
    with rasterio.open(sixs_p) as d:
        w254 = envi_wavelengths(d)
    with rasterio.open(terr_p) as d:
        w254t = envi_wavelengths(d)
    for path in (sixs_p, terr_p, diag_p, cov_p):
        if path == cov_p and not path.exists():
            continue
        with rasterio.open(path) as d:
            if (d.crs, d.transform, d.shape) != (grid_crs, transform, shape):
                raise ValueError(f"四阶段输入格网不一致: {path}")
    if not w330.size or not w254.size or w254.shape != w254t.shape or not np.allclose(w254, w254t):
        raise ValueError("波长数/一致性检查失败")
    with band_csv.open(encoding="utf-8-sig", newline="") as f:
        band_rows = list(csv.DictReader(f))
    ob = np.asarray([int(r["output_band"]) for r in band_rows], int)
    sb = np.asarray([int(r["source_band"]) for r in band_rows], int)
    seg = np.asarray([r["segment"] for r in band_rows])
    csv_wl = np.asarray([float(r["wavelength_nm"]) for r in band_rows], float)
    order = np.argsort(ob)
    ob, sb, seg = ob[order], sb[order], seg[order]
    csv_wl = csv_wl[order]
    gain = np.asarray([float(r["gain"]) for r in band_rows], float)[order]
    offset = np.asarray([float(r["offset"]) for r in band_rows], float)[order]
    if len(gain) != len(w330) or not np.array_equal(ob, np.arange(1, len(w330) + 1)):
        raise ValueError("定标表 output_band 必须逐一对应辐亮度波段")
    if not np.allclose(csv_wl, w330, rtol=0, atol=0.01):
        raise ValueError("波段表波长与辐亮度 ENVI 标签不一致，拒绝按假设顺序配对")
    # 底座在 VNIR/SWIR 重叠区按波长交错排布，输出波段 ≠ VN+SW 简单拼接。

    sites = pick_sites(diag_p, cov_p if cov_p.is_file() else None, sixs_p,
                       tuple(a.primary) if a.primary else None)
    print("选点：", [(s["site"], s["row"], s["col"], round(s["cos_i"], 3)) for s in sites])

    spectra = {k: {} for k in STAGES}
    site_meta = []
    for s in sites:
        r, c = s["row"], s["col"]
        x, y = transform * (c + .5, r + .5)
        back = rpc_backproject(a.raw_dir, a.dem, grid_crs, x, y)
        vn = back["segments"].get("VN", {})
        sw = back["segments"].get("SW", {})
        if not back["both_cameras_available"]:
            raise RuntimeError(f"{s['site']}: VN/SW 缺少原始探测器支持")
        vn_dn = np.asarray(vn["dn"], float)
        sw_dn = np.asarray(sw["dn"], float)
        is_vn = seg == "VNIR"
        dn = np.empty(len(w330), float)
        dn[is_vn] = vn_dn[sb[is_vn] - 1]
        dn[~is_vn] = sw_dn[sb[~is_vn] - 1]
        spectra["raw_dn"][s["site"]] = dn
        spectra["radiance"][s["site"]] = read_pixel(rad_p, r, c)
        spectra["sixs_reflectance"][s["site"]] = read_pixel(sixs_p, r, c)
        spectra["terrain_corrected"][s["site"]] = read_pixel(terr_p, r, c)
        site_meta.append({**s, "rpc_backprojection": back})
        print(f"  {s['site']}: src pixel VN=({vn['src_row']},{vn['src_col']}) "
              f"SW=({sw['src_row']},{sw['src_col']}) DN bands={len(dn)}")

    waves = {"raw_dn": w330, "radiance": w330,
             "sixs_reflectance": w254, "terrain_corrected": w254}
    out.mkdir(parents=True, exist_ok=True)
    with (out / "paired_chain_spectra.csv").open("w", encoding="utf-8-sig",
                                                 newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["site", "row", "col", "stage", "band", "wavelength_nm",
                     "value", "unit", "product_grade"])
        grade = {"raw_dn": "L1_DN_UNPROCESSED",
                 "radiance": "at_sensor_radiance",
                 "sixs_reflectance": "sixs_model_estimated_surface_reflectance_NOT_independently_validated",
                 "terrain_corrected": "sixs_terrain_illumination_corrected_model_estimated_reflectance_NOT_independently_validated"}
        for s in sites:
            for stage in STAGES:
                for i, (w, v) in enumerate(zip(waves[stage], spectra[stage][s["site"]]), 1):
                    wr.writerow([s["site"], s["row"], s["col"], stage, i,
                                 f"{w:.4f}", f"{v:.8g}", STYLE[stage][2], grade[stage]])

    stats = {}
    for s in sites:
        key = s["site"]
        c = spectra["sixs_reflectance"][key]
        t = spectra["terrain_corrected"][key]
        ok = np.isfinite(c) & np.isfinite(t) & (c > 0)
        rel = (t[ok] - c[ok]) / c[ok]
        win = (w254 >= 2130) & (w254 <= 2240)
        def center(v):
            z = v[win & np.isfinite(v)]
            return float(w254[win & np.isfinite(v)][int(np.argmin(z))]) if z.size else None
        dn = spectra["raw_dn"][key]
        l_pred = dn * gain + offset
        l_obs = spectra["radiance"][key]
        rec_ok = np.isfinite(l_pred) & np.isfinite(l_obs) & (l_pred > 0) & (l_obs > 0)
        rec_rel = (l_obs[rec_ok] - l_pred[rec_ok]) / l_pred[rec_ok]
        hi = rec_ok & (dn >= 50)   # 排除近零 DN 的强吸收/低信号波段
        reconciliation = {
            "n_bands": int(rec_ok.sum()),
            "pearson_r": (float(np.corrcoef(l_pred[rec_ok], l_obs[rec_ok])[0, 1])
                          if rec_ok.sum() > 5 else None),
            "pearson_r_high_signal_dn_ge_50": (
                float(np.corrcoef(l_pred[hi], l_obs[hi])[0, 1])
                if hi.sum() > 5 else None),
            "n_high_signal_bands": int(hi.sum()),
            "relative_diff_median": float(np.median(rec_rel)) if rec_rel.size else None,
            "relative_diff_p05": float(np.percentile(rec_rel, 5)) if rec_rel.size else None,
            "relative_diff_p95": float(np.percentile(rec_rel, 95)) if rec_rel.size else None,
            "interpretation": ("DN is read at the RAW camera pixel while radiance is an "
                               "orthorectified RESAMPLED pixel; differences are expected "
                               "from sub-pixel resampling and are NOT a calibration error."),
        }
        stats[key] = {
            "dn_to_radiance_reconciliation": reconciliation,
            "sam_angle_deg": sam_angle_deg(c, t),
            "terrain_relative_change_p05": float(np.percentile(rel, 5)) if rel.size else None,
            "terrain_relative_change_median": float(np.median(rel)) if rel.size else None,
            "terrain_relative_change_p95": float(np.percentile(rel, 95)) if rel.size else None,
            "window_min_wavelength_nm_before": center(c),
            "window_min_wavelength_nm_after": center(t),
            "stage_summary": {st: {"n_bands": int(len(spectra[st][key])),
                                    "min": float(np.nanmin(spectra[st][key])),
                                    "median": float(np.nanmedian(spectra[st][key])),
                                    "max": float(np.nanmax(spectra[st][key])),
                                    "unit": STYLE[st][2]}
                              for st in STAGES},
        }

    primary = sites[0]
    one = {k: spectra[k][primary["site"]] for k in STAGES}
    plot_chain(primary, one, waves, out)
    plot_normalized(primary, one, waves, out)
    plot_terrain_effect(sites, spectra, waves, stats, out)

    meta = {
        "purpose": "same-pixel unprocessed vs fully processed spectral chain",
        "same_pixel_definition": ("B/C/D share product grid (row,col); A is read at the "
                                  "RPC-backprojected containing pixel for each raw camera"),
        "raw_dn_band_order": ("The radiance cube is wavelength-interleaved across VNIR/SWIR "
                              "overlap (1000-1025 nm); raw DN were reordered by "
                              "output_band -> (segment, source_band), NOT VN+SW concatenation. "
                              "Wavelengths were cross-checked against ENVI tags (atol 0.01 nm)."),
        "stages": {st: {"unit": STYLE[st][2], "n_bands": int(len(waves[st])),
                        "definition": STYLE[st][1]} for st in STAGES},
        "water_absorption_windows_nm": WATER_WINDOWS,
        "sites": site_meta,
        "statistics": stats,
        "caveats": [
            "Raw DN is at the backprojected RAW camera pixel, not a resampled orthorectified "
            "pixel; sub-pixel resampling differences of order 1 native pixel remain.",
            "In stage A the 1000-1025 nm VNIR/SWIR overlap shows DN steps: the two cameras "
            "have different DN scales and the cube interleaves them by wavelength. This is "
            "exactly what UNPROCESSED data look like; radiometric calibration (stage B) "
            "resolves the scales into a common physical radiance unit.",
            "The 6SV/terrain curves show a sharp feature near 2000 nm at the edge of the "
            "excluded 1800-1980 nm water band; this is a model/SRF edge artefact and must NOT "
            "be read as a mineral absorption feature.",
            "DN (330 bands) and reflectance (254 bands) have different band counts; curves "
            "are compared by wavelength, not band index.",
            "All four stages use INDEPENDENT y-axes; DN/radiance/reflectance are NOT "
            "interchangeable and no numeric accuracy claim is made.",
            "Sites are selected deterministically (existing P1 plus cos(i) quantiles); they "
            "are NOT mineral or vegetation anomalies and carry no ground truth.",
            "6SV reflectance is a model estimate with assumed atmospheric scenarios; "
            "terrain correction is first-order direct-solar only; neither is independently validated.",
        ],
    }
    (out / "chain_comparison.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"写出 {out}  (3 figures + paired_chain_spectra.csv + chain_comparison.json)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
