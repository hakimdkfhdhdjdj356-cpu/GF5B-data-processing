"""Build AOI QA, analysis grids and the local asset catalog after the AHSI run."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import fiona
import numpy as np
import rasterio as rio
from rasterio.enums import Resampling
from rasterio.features import geometry_mask
from rasterio.transform import RPCTransformer
from rasterio.warp import reproject, transform, transform_geom
from rasterio.windows import Window
from scipy.ndimage import binary_erosion, median_filter
from shapely.geometry import shape

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ahsi.atmos import METHOD_INFO
from ahsi.config import aux_path, load_config, load_yaml
from ahsi.core import save_json, write_csv
from ahsi.scene import parse_scene

WINDOWS = {
    "clay_kaolinite_doublet_context": [2135, 2225],
    "white_mica_illite_AlOH": [2175, 2240],
    "chlorite_epidote_FeMgOH_context": [2225, 2370],
    "carbonates_CO3": [2290, 2405],
    "longwave_context": [2405, 2450],
}
NODATA = -9999.0
MASK_BITS = {"inside_aoi": 1, "observed": 2, "all_retained_bands_valid": 4,
             "some_reflectance_gt_1": 8, "some_reflectance_gt_2": 16, "source_dn4095_suspect": 32}


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def aoi_on(dataset, feature, crs):
    geometry = transform_geom(crs, dataset.crs, feature)
    return geometry_mask([geometry], dataset.shape, dataset.transform, invert=True)


def output_profile(grid, count, dtype="float32", nodata=NODATA):
    return dict(driver="GTiff", width=grid.width, height=grid.height, crs=grid.crs,
                transform=grid.transform, count=count, dtype=dtype, nodata=nodata,
                tiled=True, blockxsize=128, blockysize=128, compress="deflate",
                predictor=3 if dtype == "float32" else 2)


def write_mask(path, grid, inside, observed, complete, gt1=None, gt2=None, suspect=None):
    flags = inside.astype("uint8") | (observed.astype("uint8") << 1) | (complete.astype("uint8") << 2)
    if gt1 is not None:
        flags |= (gt1.astype("uint8") << 3) | (gt2.astype("uint8") << 4)
    if suspect is not None:
        flags |= (suspect.astype("uint8") << 5)
    with rio.open(path, "w", **output_profile(grid, 1, "uint8", 255)) as dst:
        dst.write(flags, 1)
        dst.update_tags(**MASK_BITS, note="Bits stored as integer masks; no GF5B cloud/snow/shadow product")
    return {"aoi_pixels": int(inside.sum()), "observed": int(observed.sum()),
            "all_retained_bands_valid": int(complete.sum()),
            "source_dn4095_suspect_pixels": int(suspect.sum()) if suspect is not None else 0}


def highpass_noise(a, valid):
    ok = valid[:-1, :-1] & valid[1:, 1:] & valid[:-1, 1:] & valid[1:, :-1]
    residual = ((a[:-1, :-1] + a[1:, 1:] - a[:-1, 1:] - a[1:, :-1]) / 2)[ok]
    return float(1.4826 * np.median(np.abs(residual - np.median(residual)))) if residual.size else None


def native_detector_qa(scene, grid, inside, dem):
    """Use unique detector cells sampled by the center grid for local QC."""
    y, x = np.where(inside)
    xs, ys = grid.transform * (x + .5, y + .5)
    lons, lats = transform(grid.crs, "EPSG:4326", xs.tolist(), ys.tolist())
    output = {}
    suspect = np.zeros(inside.shape, bool)
    for segment in scene.segments:
        with rio.open(segment.image) as src, RPCTransformer(src.rpcs, RPC_DEM=str(dem)) as rpc:
            rows, cols = map(np.asarray, rpc.rowcol(lons, lats, op=lambda v: v))
            valid = np.isfinite(rows) & np.isfinite(cols)
            centers = np.flatnonzero(valid)
            rows = np.floor(rows[valid]).astype(int)
            cols = np.floor(cols[valid]).astype(int)
            valid = (rows >= 0) & (cols >= 0) & (rows < src.height) & (cols < src.width)
            centers = centers[valid]
            rows, cols = rows[valid], cols[valid]
            row0, col0 = int(rows.min()), int(cols.min())
            window = Window(col0, row0, int(cols.max()) - col0 + 1, int(rows.max()) - row0 + 1)
            support = np.zeros((int(window.height), int(window.width)), bool)
            support[rows - row0, cols - col0] = True
            data = src.read(window=window)
            support &= np.any(data != 0, axis=0)
            cap = np.any(data == 4095, axis=0)[rows - row0, cols - col0]
            suspect[y[centers[cap]], x[centers[cap]]] = True
            # Remove broad terrain gradients before aggregating residuals along detector columns.
            counts = support.sum(axis=0)
            core = counts >= max(30, .5 * counts.max())
            container_max = np.iinfo(data.dtype).max
            for i, band in enumerate(data, 1):
                finite = band[support]
                a = band.astype("float32")
                hp = a - median_filter(a, size=(1, 7), mode="nearest")
                # Use only valid cells in the along-track average of each detector column.
                profile = np.nanmean(np.where(support[:, core], hp[:, core], np.nan), axis=0)
                level = max(abs(float(np.median(finite))), 1.0)
                output[segment.name, i] = dict(
                    detector_support=int(finite.size), detector_max_DN=int(finite.max()),
                    detector_zero_fraction=float(np.mean(finite == 0)),
                    detector_negative_fraction=float(np.mean(finite < 0)),
                    detector_container_cap_count=int(np.sum(finite == container_max)),
                    detector_4095_count=int(np.sum(finite == 4095)),
                    detector_stripe_column_p95_fraction=float(np.quantile(np.abs(profile), .95) / level),
                )
    return output, suspect


def read_band(ds, band):
    return ds.read(band, masked=True).astype("float32").filled(np.nan)


def band_qa(scene, rad, refl, inside, observed, dem):
    rows = [dict(row) for row in scene.band_table]
    detector, suspect = native_detector_qa(scene, rad, observed, dem)
    analysis = 0
    for row in rows:
        i = row["output_band"]
        a = read_band(rad, i)
        ok = inside & np.isfinite(a)
        positive = ok & (a > 0)
        noise = highpass_noise(a, positive)
        median = float(np.median(a[positive])) if positive.any() else None
        row.update(detector[row["segment"], row["source_band"]])
        row.update(aoi_finite_radiance_pixels=int(ok.sum()),
                   aoi_positive_radiance_fraction=float(positive.sum() / inside.sum()),
                   aoi_mean_radiance=float(np.mean(a[ok])) if ok.any() else None,
                   aoi_median_positive_radiance=median, spatial_residual_mad_radiance=noise,
                   spatial_residual_to_median=(noise / median if noise is not None and median else None),
                   analysis_band=None, aoi_valid_analysis_pixels=0, effective_valid=False,
                   aoi_mean_analysis=None, aoi_median_analysis=None)
        if row["analysis_candidate"]:
            analysis += 1
            row["analysis_band"] = analysis
            a = read_band(refl, analysis)
            ok = observed & np.isfinite(a)
            tags = refl.tags(analysis)
            if (abs(float(tags["wavelength_nm"]) - row["wavelength_nm"]) > 1e-4
                    or int(tags["source_native_band"]) != row["source_band"]
                    or tags["source_camera"] != row["segment"]):
                raise ValueError(f"Band map mismatch at QUAC band {analysis}")
            row["quac_bbl"] = int(tags["quac_bbl"])
            row["aoi_valid_analysis_pixels"] = int(ok.sum())
            row["effective_valid"] = bool(row["quac_bbl"] and ok.sum() >= .98 * observed.sum())
            row["aoi_mean_analysis"] = float(np.mean(a[ok])) if ok.any() else None
            row["aoi_median_analysis"] = float(np.median(a[ok])) if ok.any() else None
            if not row["effective_valid"]:
                row["screening_reason"] = "quac_bad_or_aoi_valid_fraction_below_0.98"
    envi = refl.tags(ns="ENVI")
    for field, csv_key in (("wavelength", "wavelength_nm"), ("fwhm", "fwhm_nm")):
        expected = [r[csv_key] for r in rows if r["analysis_candidate"]]
        actual = np.fromstring(envi[field].strip("{} "), sep=",")
        if len(expected) != len(actual) or not np.allclose(expected, actual, rtol=0, atol=.001):
            raise ValueError(f"QUAC ENVI {field} differs from the source mapping")
    if analysis != refl.count:
        raise ValueError("QUAC band count does not match the source band table")
    return rows, suspect


def resample(a, src, dst, method, nodata=np.nan):
    out = np.full(dst.shape, nodata, "float32")
    reproject(a, out, src_transform=src.transform, src_crs=src.crs, src_nodata=nodata,
              dst_transform=dst.transform, dst_crs=dst.crs, dst_nodata=nodata, resampling=method)
    return out


def write_analysis_grid(src, grid, inside, native_observed, rows, path, qa,
                        native_suspect, fractions=None, min_fraction=.75):
    """Average directly to EMIT support or reproject once to project 30 m."""
    average = fractions is not None
    method = Resampling.average if average else Resampling.bilinear
    obs = resample(native_observed.astype("float32"), src, grid,
                   Resampling.average if average else Resampling.nearest, -1.0)
    observed = inside & (obs >= (min_fraction if average else 1.0))
    complete, gt1, gt2 = observed.copy(), np.zeros(grid.shape, bool), np.zeros(grid.shape, bool)
    selected = [r for r in rows if r["analysis_band"] is not None]
    profile = output_profile(grid, len(selected))
    with rio.open(path, "w", **profile) as dst:
        support_dst = rio.open(fractions, "w", **profile) if average else None
        try:
            dst.update_tags(**src.tags(), resampling=method.name, source_file=str(src.name))
            dst.update_tags(ns="ENVI", **src.tags(ns="ENVI"))
            for r in selected:
                i = r["analysis_band"]
                a = read_band(src, i)
                a[~native_observed] = np.nan
                out = resample(a, src, grid, method)
                if average:
                    fraction = resample(np.isfinite(a).astype("float32"), src, grid, method, -1.0)
                    support_dst.write(np.where(inside & (fraction >= 0), fraction, NODATA), i)
                    support_dst.set_band_description(i, src.descriptions[i - 1])
                    out[fraction < min_fraction] = np.nan
                valid = observed & np.isfinite(out)
                if r["effective_valid"]:
                    complete &= valid
                    gt1 |= valid & (out > 1)
                    gt2 |= valid & (out > 2)
                dst.write(np.where(valid, out, NODATA), i)
                dst.set_band_description(i, src.descriptions[i - 1])
                dst.update_tags(i, **src.tags(i))
            if support_dst is not None:
                support_dst.update_tags(meaning="Area-weighted fraction of valid native support by GF5B band",
                                        minimum_fraction=min_fraction)
        finally:
            if support_dst is not None:
                support_dst.close()
    suspect = observed & (resample(native_suspect.astype("float32"), src, grid,
                                   Resampling.average if average else Resampling.nearest, -1.0) > 0)
    return write_mask(qa, grid, inside, observed, complete, gt1, gt2, suspect)


def alignment_check(a, b, valid, pixel_size_m, radius):
    """Relative image similarity with fixed domain for every tested offset."""
    a = np.hypot(*np.gradient(a))
    b = np.hypot(*np.gradient(b))
    finite = valid & np.isfinite(a) & np.isfinite(b)
    support = binary_erosion(finite, np.ones((2 * radius + 1, 2 * radius + 1), bool))
    rr, cc = np.where(support)
    if len(rr) < 100:
        return {"status": "insufficient_support", "pixels": len(rr)}
    qa = a[rr, cc]
    mid_r, mid_c = np.median(rr), np.median(cc)
    quadrants = [((rr >= mid_r) == south) & ((cc >= mid_c) == east)
                 for south in (False, True) for east in (False, True)]
    results = []
    for dr in range(-radius, radius + 1):
        for dc in range(-radius, radius + 1):
            q = b[rr + dr, cc + dc]
            results.append((dr, dc, float(np.corrcoef(qa, q)[0, 1]),
                            [float(np.corrcoef(qa[m], q[m])[0, 1]) for m in quadrants]))
    best = max(results, key=lambda x: x[2])
    zero = next(row for row in results if row[:2] == (0, 0))
    return {"pixels": len(rr), "step_m": pixel_size_m, "search_radius_pixels": radius,
            "offset_convention": "sample reference at (row+dr, col+dc) versus the first image",
            "zero_score": zero[2], "best_score": best[2], "best_dr_dc": list(best[:2]),
            "quadrant_best_dr_dc": [list(max(results, key=lambda x: x[3][i])[:2]) for i in range(4)]}


def emit_clear(qa):
    names = ("cloud_flag", "cirrus_flag", "water_flag", "spacecraft_flag",
             "dilated_cloud_flag", "aggregate_flag")
    clean = np.ones(qa.shape, bool)
    for name in names:
        clean &= read_band(qa, qa.descriptions.index(name) + 1) == 0
    return clean


def prepare_emit_support(assets, native, native_obs, rows, aoi, aoi_crs, emit, emit_qa, emit_props, suspect):
    data_path = assets / "GF5B_on_EMIT_60m.tif"
    qa_path = assets / "GF5B_on_EMIT_60m_QA.tif"
    fractions = assets / "GF5B_on_EMIT_60m_support_fraction.tif"
    inside = aoi_on(emit, aoi, aoi_crs)
    if (emit.shape, emit.transform, emit.crs) != (emit_qa.shape, emit_qa.transform, emit_qa.crs):
        raise ValueError("EMIT QA grid does not match its reflectance raster")
    counts = write_analysis_grid(native, emit, inside, native_obs, rows, data_path, qa_path, suspect, fractions)
    clean = inside & emit_clear(emit_qa)
    wavelengths = np.asarray(emit_props["reflectance_wavelengths"])
    if emit.count != len(wavelengths):
        raise ValueError("EMIT wavelength and raster counts differ")
    common_path = assets / "EMIT_GF5B_common_observation_mask_60m.tif"
    alignment = []
    with rio.open(qa_path) as q, rio.open(data_path) as gf:
        gf_complete = (q.read(1) & MASK_BITS["all_retained_bands_valid"]) != 0
        emit_spectral_available = np.ones(emit.shape, bool)
        for i in emit.indexes:
            emit_spectral_available &= np.isfinite(read_band(emit, i))
        common = clean & emit_spectral_available & gf_complete
        with rio.open(common_path, "w", **output_profile(emit, 1, "uint8", 255)) as dst:
            dst.write(np.where(inside, common, 255).astype("uint8"), 1)
            dst.update_tags(meaning="1=joint spectral coverage and EMIT flag-clear; 0=unsupported; 255=outside AOI",
                            GF5B_cloud_status="unknown")
        for w in (650, 860):
            r = min((r for r in rows if r["effective_valid"]), key=lambda r: abs(r["wavelength_nm"] - w))
            i = int(np.argmin(abs(wavelengths - w))) + 1
            result = alignment_check(read_band(gf, r["analysis_band"]), read_band(emit, i), common, 60, 3)
            alignment.append(dict(gf5b_nm=r["wavelength_nm"], emit_nm=float(wavelengths[i - 1]), **result))
    counts.update(emit_qa_clear=int(clean.sum()), emit_all_285_finite=int((inside & emit_spectral_available).sum()),
                  common=int(common.sum()), alignment=alignment)
    return data_path, qa_path, fractions, common_path, counts


def sentinel_alignment(root, s2, gf_path, common_path, emit, rows, emit_props):
    """Check whether the GF5B/EMIT relative offset also appears against local S2A."""
    wavelengths = np.asarray(emit_props["reflectance_wavelengths"])
    results = []
    with (rio.open(root / s2["qa"]) as s2qa, rio.open(gf_path) as gf,
          rio.open(common_path) as common_ds):
        common = common_ds.read(1) == 1
        s2_clean = read_band(s2qa, 1) == 1
        for wave, s2_key, band in ((650, "sr10", "B4"), (860, "sr20", "B8A")):
            with rio.open(root / s2[s2_key]) as src:
                mask = resample(s2_clean.astype("float32"), s2qa, src, Resampling.nearest, -1.0) == 1
                values = read_band(src, src.descriptions.index(band) + 1) * .0001
                values[~mask] = np.nan
                support = resample(np.isfinite(values).astype("float32"), src, emit, Resampling.average, -1.0)
                ref = resample(values, src, emit, Resampling.average)
                valid = common & (support >= .75)
            gf_band = min((r for r in rows if r["effective_valid"]), key=lambda r: abs(r["wavelength_nm"] - wave))
            emit_band = int(np.argmin(abs(wavelengths - wave))) + 1
            for sensor, a in (("GF5B", read_band(gf, gf_band["analysis_band"])),
                              ("EMIT", read_band(emit, emit_band))):
                results.append(dict(sensor=sensor, reference=f"Sentinel-2A {band}",
                                    **alignment_check(a, ref, valid, 60, 3)))
    return results


def representative_spectra(path, ds, inside, rows):
    rr, cc = np.where(inside)
    records = []
    for fx, fy in ((.2, .2), (.8, .2), (.5, .5), (.2, .8), (.8, .8)):
        r = np.quantile(rr, fy)
        c = np.quantile(cc, fx)
        j = np.argmin((rr - r)**2 + (cc - c)**2)
        row, col = int(rr[j]), int(cc[j])
        vals = ds.read(window=Window(col, row, 1, 1), masked=True).astype("float32")[:, 0, 0]
        for band in rows:
            i = band["analysis_band"]
            if i is not None:
                records.append(dict(site=f"P{len(records) // ds.count + 1}", row=row, col=col,
                                    wavelength_nm=band["wavelength_nm"], analysis_band=i,
                                    value=None if np.ma.is_masked(vals[i - 1]) else float(vals[i - 1])))
    write_csv(path, records)


def raster_info(path):
    with rio.open(path) as ds:
        return {"crs": str(ds.crs) if ds.crs else None, "transform": list(ds.transform)[:6],
                "shape": list(ds.shape), "bands": ds.count, "dtype": ds.dtypes[0],
                "nodata": ds.nodata if ds.nodata is None or np.isfinite(ds.nodata) else None,
                "nodata_encoding": "NaN" if ds.nodata is not None and np.isnan(ds.nodata) else None,
                "resolution": list(ds.res) if ds.crs else None, "has_rpc": bool(ds.rpcs)}


def build(cfg, run):
    assets = run / "assets"
    xml, = Path(cfg["input_root"]).glob("*.xml")
    scene = parse_scene(xml, cfg)
    sd = run / "scenes" / scene.scene_id
    rad_result = read_json(sd / "radiance_summary.json")
    refl_result = read_json(sd / "stage_reflectance.json")["result"]
    quac = read_json(sd / "stage_quac.json")["result"]
    if not quac["successful"] or quac["backend_method"] != "envi_quac":
        raise ValueError("This asset configuration expects completed QUAC")
    grade = METHOD_INFO["envi_quac"]["product_grade"]
    root = Path(cfg["project_root"])
    local = cfg["scene_assets"]
    aoi_path = root / local["aoi"]
    emit_path = root / local["emit_reflectance"]
    emit_qa_path = root / local["emit_qa"]
    emit_meta = root / local["emit_metadata"]
    aster = root / local["aster_native_dir"]
    aoi_qa_path = assets / "GF5B_native30m_AOI_QA.tif"
    project_path = assets / "GF5B_project30m_QUAC.tif"
    project_qa = assets / "GF5B_project30m_AOI_QA.tif"
    band_path = assets / "GF5B_wavelength_and_valid_bands.csv"
    with fiona.open(aoi_path) as f:
        aoi = next(iter(f)).geometry.__geo_interface__
        aoi_crs = f.crs
    paths = {"radiance": Path(rad_result["radiance_file"]), "analysis_native": Path(refl_result["reflectance_file"])}
    with (rio.open(paths["radiance"]) as rad, rio.open(paths["analysis_native"]) as native,
          rio.open(rad_result["coverage_mask"]) as coverage, rio.open(emit_path) as emit,
          rio.open(emit_qa_path) as emit_qa, rio.open(aux_path(cfg, "rpc_dem")) as project):
        if (rad.shape, rad.transform, rad.crs) != (native.shape, native.transform, native.crs):
            raise ValueError("QUAC and radiance grids differ")
        if project.crs != local["project_crs"]:
            raise ValueError("Project grid and configured CRS differ")
        inside = aoi_on(native, aoi, aoi_crs)
        observed = inside & (coverage.read(1) == 1)
        rows, suspect = band_qa(scene, rad, native, inside, observed, aux_path(cfg, "rpc_dem"))
        retained = [r["analysis_band"] for r in rows if r["effective_valid"]]
        cube = native.read(retained, masked=True)
        complete = observed & ~np.ma.getmaskarray(cube).any(axis=0)
        counts = {"native": write_mask(aoi_qa_path, native, inside, observed, complete,
                                       observed & np.any(cube.filled(0) > 1, axis=0),
                                       observed & np.any(cube.filled(0) > 2, axis=0), suspect)}
        project_inside = aoi_on(project, aoi, aoi_crs)
        counts["project"] = write_analysis_grid(
            native, project, project_inside, coverage.read(1) == 1, rows, project_path, project_qa, suspect)
        emit_props = read_json(emit_meta)["properties"]
        common_data, common_qa, fractions, common_mask, counts["emit_support"] = prepare_emit_support(
            assets, native, observed, rows, aoi, aoi_crs, emit, emit_qa, emit_props, suspect)
        counts["sentinel_alignment"] = sentinel_alignment(
            root, local["sentinel2"], common_data, common_mask, emit, rows, emit_props)
        representative_spectra(assets / "GF5B_representative_spectra.csv", native, complete, rows)
        by_camera = [[r for r in rows if r["segment"] == segment] for segment in ("VNIR", "SWIR")]
        vn, sw = min(((v, s) for v in by_camera[0] for s in by_camera[1]),
                     key=lambda pair: abs(pair[0]["wavelength_nm"] - pair[1]["wavelength_nm"]))
        a, b = [read_band(rad, r["output_band"]) for r in (vn, sw)]
        counts["camera_alignment"] = alignment_check(a, b, observed, 30, 3)
        counts["camera_alignment"]["wavelengths_nm"] = [vn["wavelength_nm"], sw["wavelength_nm"]]
        aoi_geom = shape(transform_geom(aoi_crs, project.crs, aoi))
        counts["aoi_area_m2"] = float(aoi_geom.area)
        counts["native_any_above_1"] = int((observed & np.any(cube.filled(0) > 1, axis=0)).sum())
        counts["native_any_above_2"] = int((observed & np.any(cube.filled(0) > 2, axis=0)).sum())

    write_csv(band_path, rows)
    save_json(assets / "GF5B_valid_band_mask.json", {
        "order": "ascending wavelength, original camera and 1-based band retained",
        "radiance_330": [int(r["effective_valid"]) for r in rows],
        "quac_output": [int(r["effective_valid"]) for r in rows if r["analysis_band"] is not None],
        "definition": "wavelength screening, QUAC BBL=1, finite on >=98% of jointly observed AOI centers",
        "not_instrument_certification": True,
    })
    windows = {name: {"bounds_nm": pair, "bands": [r["analysis_band"] for r in rows if r["effective_valid"]
               and pair[0] <= r["wavelength_nm"] <= pair[1]]} for name, pair in WINDOWS.items()}
    summary = dict(counts, windows=windows, source_bands=len(rows), retained_bands=len(retained),
                   excluded_reasons=dict(Counter(r["screening_reason"] for r in rows if not r["effective_valid"])))
    save_json(assets / "GF5B_spatial_spectral_summary.json", summary)

    def relative(path):
        return Path(os.path.relpath(path, assets)).as_posix()

    def product(path, qa=None, mapping=None, **metadata):
        return dict(path=relative(path), grid=raster_info(path), qa=relative(qa) if qa else None,
                    band_column=mapping, **metadata)

    source = {"directory": str(xml.parent), "xml": relative(xml), "inventory": [
        {"file": x.name, "bytes": x.stat().st_size} for x in sorted(xml.parent.iterdir()) if x.is_file()],
        "xml_fields": {el.tag: el.text for el in ET.parse(xml).getroot()},
        "scene": scene.scene_id, "ogp": scene.ogp,
        "srfs": "Center wavelength and FWHM only; response curves and detector spectral smile unavailable",
        "footprint": relative(scene.files["footprint"]), "pixel_quality_from_supplier": None}
    for segment in scene.segments:
        source[segment.name] = {"raw": relative(segment.image), "rpc": relative(segment.rpc_file),
                                "wavelength_fwhm": relative(segment.spectral_file),
                                "calibration": relative(segment.radcal_file)}
    aster_products = {}
    for mode, count in (("VNIR", 3), ("SWIR", 6)):
        path = aster / f"ASTER_{mode}_native_AOI_reflectance.tif"
        with rio.open(path) as ds:
            band_ids = [description.split("_", 1)[0] for description in ds.descriptions]
        aster_products[mode] = product(
            path, aster / f"ASTER_{mode}_native_AOI_QA.tif", bands=count, scale=1.0,
            spectral_bands=[dict(id=b, passband_nm=local["aster_bands_nm"][b],
                                 nominal_midpoint_nm=sum(local["aster_bands_nm"][b]) / 2)
                            for b in band_ids],
            qa_status="First plane retains clear code 0; second plane stored as raw codes")
    s2 = local["sentinel2"]
    save_json(assets / "spectral_asset_catalog.json", {
        "schema_version": 1, "aoi": relative(aoi_path), "project_crs": local["project_crs"],
        "qa_bits": MASK_BITS, "sensors": {
            "GF5B": {"sensor": "AHSI", "acquired": scene.meta["start_time"], "time_zone": "unspecified in XML",
                "source_grade": scene.meta["product_level"], "nominal_gsd_m": 30, "source": source,
                "spectral_table": relative(band_path), "valid_band_mask": "GF5B_valid_band_mask.json",
                "coverage_summary": summary, "diagnostic_windows": windows,
                "surface_reflectance_independently_validated": False,
                "products": {
                    "radiance_native": product(paths["radiance"], aoi_qa_path, "output_band",
                        grade="at_sensor_radiance", units=cfg["sensors"]["calibration"]["units"]),
                    "analysis_native": product(paths["analysis_native"], aoi_qa_path, "analysis_band",
                        grade=grade, units="dimensionless"),
                    "analysis_project": product(project_path, project_qa, "analysis_band",
                        grade=grade, units="dimensionless"),
                    "emit_support": product(common_data, common_qa, "analysis_band",
                        grade=grade, units="dimensionless", support_fraction=relative(fractions),
                        minimum_support_fraction=.75, joint_observation_mask=relative(common_mask)),
                    **{segment.name.lower() + "_dn": product(segment.image, grade="LEVEL1_DN",
                          camera=segment.name) for segment in scene.segments},
                }},
            "EMIT": {"sensor": "EMIT", "acquired": "2023-02-04T04:04:48Z", "source_grade": "NASA L2A RFL",
                "grid_role": "60 m project delivery grid, distinct from native detector grid",
                "spectral_metadata": relative(emit_meta), "srf_status": "metadata center wavelength and FWHM",
                "local_v02_band_screening": relative(root / local["emit_v02_bands"]),
                "comparison_band_selection": "All finite delivery bands; v0.2 retained-band policy stored separately",
                "products": {"project": product(emit_path, emit_qa_path, grade="L2A_surface_reflectance",
                    uncertainty=relative(emit_path.with_name("EMIT_uncertainty_285bands_60m.tif")))}},
            "ASTER": {"sensor": "Terra ASTER", "acquired": "2006-03-31T03:55:10Z",
                      "source_grade": "AST_07XT.004 surface reflectance", "products": aster_products,
                      "native_asset_metadata": relative(root / local["aster_metadata"]),
                      "srf_reference": "https://asterweb.jpl.nasa.gov/content/01_mission/03_instrument/01_Characteristics/default.htm",
                      "role": "separate historical observation"},
            "Sentinel-2": {"sensor": "MSI", "platform": s2["platform"], "acquired": s2["date"], "source_grade": "L2A surface reflectance",
                           "spectral_metadata": relative(root / s2["metadata"]),
                           "spectral_bands": {k: dict(equivalent_wavelength_nm=v[0], bandwidth_nm=v[1])
                                              for k, v in s2["bands_nm"].items()},
                           "spectral_reference": "https://sentiwiki.copernicus.eu/web/s2-mission",
                           "srf_reference": "https://sentiwiki.copernicus.eu/web/s2-documents",
                           "products": {f"native_{res}m": product(root / item, root / s2["qa"], reflectance_scale=.0001,
                                                                    categorical_SCL_preserved=(res == 20))
                                        for res, item in ((10, s2["sr10"]), (20, s2["sr20"]))}},
        },
        "comparison": {"method": "GF5B area-weighted aggregation to the existing EMIT delivery grid",
                       "min_GF5B_band_support_fraction": .75, "band_stacking": False,
                       "cross_sensor_reflectance_calibration": False, "temporal_difference_allowed": False,
                       "co_registration_status": "relative offset observed; not corrected"},
    })
    print(assets / "spectral_asset_catalog.json")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("config", type=Path)
    cli.add_argument("run", type=Path)
    cli.add_argument("--asset-config", type=Path)
    args = cli.parse_args()
    cfg = load_config(args.config)
    cfg["scene_assets"].update(load_yaml(args.asset_config or args.config.with_name("scene_assets.yaml")))
    build(cfg, args.run.resolve())
