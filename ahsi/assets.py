"""Local spectral data access for registered AHSI scene products."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
from pathlib import Path

import fiona
import numpy as np
import rasterio
from rasterio.errors import WindowError
from rasterio.features import geometry_mask, geometry_window
from rasterio.warp import transform_geom
from rasterio.windows import Window


@dataclass
class Spectra:
    data: np.ndarray
    wavelength_nm: np.ndarray
    fwhm_nm: np.ndarray
    indexes: np.ndarray
    band_accepted: np.ndarray
    qa: np.ndarray
    support: np.ndarray
    transform: object
    crs: str | None
    window: Window
    grade: str
    terrain_qa: np.ndarray | None = None

    def payload(self):
        result = dict(raw=self.data, wave=self.wavelength_nm, fwhm_nm=self.fwhm_nm,
                    band_index=self.indexes, band_accepted=self.band_accepted, qa=self.qa,
                    mask=self.support, transform=np.array(self.transform)[:6],
                    crs=self.crs or "", window=self.window.flatten(), grade=self.grade)
        if self.terrain_qa is not None:
            result["terrain_qa"] = self.terrain_qa
        return result


class SpectralAsset:
    """The catalog owns all asset locations, band maps and QA semantics."""

    def __init__(self, config):
        config = Path(config).resolve()
        self.config = json.loads(config.read_text(encoding="utf-8-sig"))
        self.catalog_path = (config.parent / self.config["catalog"]).resolve()
        self.catalog = json.loads(self.catalog_path.read_text(encoding="utf-8-sig"))
        if self.catalog["schema_version"] != 1:
            raise ValueError("Unsupported spectral catalog schema")
        self.gf = self.catalog["sensors"]["GF5B"]
        with self.path(self.gf["spectral_table"]).open(encoding="utf-8-sig", newline="") as stream:
            self.bands = list(csv.DictReader(stream))

    def path(self, item):
        return (self.catalog_path.parent / item).resolve()

    def product(self, name):
        if name not in self.config["products"]:
            raise ValueError(f"{name!r} is unavailable for {self.config['role']}")
        return self.gf["products"][name]

    def _selected(self, product, diagnostic, retained):
        if diagnostic is not None:
            low, high = self.gf["diagnostic_windows"][diagnostic]["bounds_nm"]
        else:
            low, high = -np.inf, np.inf
        camera = product.get("camera")
        key = "source_band" if camera else product["band_column"]
        rows = [r for r in self.bands if r[key] and (not camera or r["segment"] == camera)
                and low <= float(r["wavelength_nm"]) <= high
                and (not retained or r["effective_valid"] == "True")]
        if not rows:
            raise ValueError("No requested bands are available in this product")
        return rows, [int(r[key]) for r in rows]

    @staticmethod
    def _window(ds, window):
        window = window or Window(0, 0, ds.width, ds.height)
        x, y, width, height = window.flatten()
        if (any(int(v) != v for v in (x, y, width, height)) or min(x, y) < 0
                or min(width, height) < 1 or x + width > ds.width or y + height > ds.height):
            raise ValueError("Window must have integer offsets and fit inside the raster")
        return window

    def describe(self):
        return {"catalog": str(self.catalog_path), "role": self.config["role"],
                "sensor": self.gf["sensor"], "acquired": self.gf["acquired"],
                "source_grade": self.gf["source_grade"], "nominal_gsd_m": self.gf["nominal_gsd_m"],
                "coverage": self.gf["coverage_summary"], "diagnostic_windows": self.gf["diagnostic_windows"],
                "qa_bits": self.catalog["qa_bits"], "source": self.gf["source"],
                "available_products": self.gf["products"]}

    def read(self, name=None, window=None, diagnostic=None, retained=None):
        name = name or self.config.get("default_product")
        product = self.product(name)
        sixs = product["grade"].startswith("sixs_")
        if retained is None:
            retained = sixs or name in ("analysis_native", "analysis_project", "emit_support")
        if product.get("camera") and window is None:
            raise ValueError("Raw detector reads require a bounded pixel window")
        rows, indexes = self._selected(product, diagnostic, retained)
        with rasterio.open(self.path(product["path"])) as ds:
            win = self._window(ds, window)
            bbl = ds.tags(ns="ENVI").get("bbl")
            if sixs and bbl:
                accepted = [int(v) == 1 for v in bbl.strip("{} ").replace(",", " ").split()]
                if len(accepted) != ds.count:
                    raise ValueError("6SV BBL and data band counts differ")
                if retained:
                    rows = [r for r, i in zip(rows, indexes) if accepted[i - 1]]
                    indexes = [i for i in indexes if accepted[i - 1]]
            else:
                accepted = [True] * ds.count
            if not indexes:
                raise ValueError("No accepted bands are available in this product")
            data = ds.read(indexes, window=win, masked=True).astype("float32").filled(np.nan)
            if product["qa"]:
                with rasterio.open(self.path(product["qa"])) as q:
                    if (q.shape, q.crs, q.transform) != (ds.shape, ds.crs, ds.transform):
                        raise ValueError("QA and data grids differ")
                    qa = q.read(1, window=win, masked=True).filled(0)
                support = (qa & self.catalog["qa_bits"]["observed"]) != 0
            else:
                support = np.any(ds.read(window=win) != 0, axis=0)
                qa = np.zeros(support.shape, dtype="uint8")
                qa[support] = self.catalog["qa_bits"]["observed"]
            data[~np.isfinite(data)] = np.nan
            data[:, ~support] = np.nan
            terrain_qa = None
            if sixs:
                # The original catalog QA supplies shared AOI geometry, VN/SW coverage
                # and DN4095 only. Product-specific reflectance bits are recomputed.
                bits = self.catalog["qa_bits"]
                names = ("all_retained_bands_valid", "some_reflectance_gt_1", "some_reflectance_gt_2")
                qa &= np.uint8(255 ^ sum(bits[n] for n in names))
                _, eligible = self._selected(product, None, True)
                eligible = [i for i in eligible if accepted[i - 1]]
                if eligible:
                    vals = data if indexes == eligible else ds.read(
                        eligible, window=win, masked=True).astype("float32").filled(np.nan)
                    finite = np.isfinite(vals)
                    qa[support & finite.all(axis=0)] |= bits[names[0]]
                    qa[support & np.any(finite & (vals > 1), axis=0)] |= bits[names[1]]
                    qa[support & np.any(finite & (vals > 2), axis=0)] |= bits[names[2]]
                if product.get("terrain_qa"):
                    with rasterio.open(self.path(product["terrain_qa"])) as q:
                        if (q.shape, q.crs, q.transform) != (ds.shape, ds.crs, ds.transform):
                            raise ValueError("Terrain QA and data grids differ")
                        terrain_qa = q.read(1, window=win, masked=True).filled(255)
            return Spectra(data, np.array([float(r["wavelength_nm"]) for r in rows]),
                           np.array([float(r["fwhm_nm"]) for r in rows]), np.asarray(indexes),
                           np.array([r["effective_valid"] == "True" and accepted[i - 1]
                                     for r, i in zip(rows, indexes)]), qa, support,
                           ds.window_transform(win), str(ds.crs) if ds.crs else None, win, product["grade"],
                           terrain_qa)

    def read_emit_pair(self, window=None, diagnostic=None):
        gf = self.read("emit_support", window, diagnostic)
        product = self.product("emit_support")
        emit = self.catalog["sensors"]["EMIT"]
        meta = json.loads(self.path(emit["spectral_metadata"]).read_text(encoding="utf-8-sig"))["properties"]
        wave = np.asarray(meta["reflectance_wavelengths"])
        fwhm = np.asarray(meta["reflectance_fwhm"])
        indexes = np.arange(1, len(wave) + 1)
        if diagnostic:
            lo, hi = self.gf["diagnostic_windows"][diagnostic]["bounds_nm"]
            selected = (wave >= lo) & (wave <= hi)
            wave, fwhm, indexes = wave[selected], fwhm[selected], indexes[selected]
        reference = emit["products"]["project"]
        with (rasterio.open(self.path(reference["path"])) as ds,
              rasterio.open(self.path(product["joint_observation_mask"])) as q,
              rasterio.open(self.path(product["support_fraction"])) as sf,
              rasterio.open(self.path(reference["qa"])) as emit_q):
            if any((a.shape, a.crs, a.transform) != (ds.shape, ds.crs, ds.transform)
                   for a in (q, sf, emit_q)) or (str(ds.crs), ds.window_transform(gf.window)) != (gf.crs, gf.transform):
                raise ValueError("Common-support grids differ")
            mask = q.read(1, window=gf.window, masked=True).filled(0) == 1
            raw = ds.read(indexes.tolist(), window=gf.window, masked=True).astype("float32").filled(np.nan)
            fractions = sf.read(gf.indexes.tolist(), window=gf.window, masked=True).filled(np.nan)
            flags = emit_q.read(window=gf.window, masked=True).filled(np.nan)
            qa_names = np.array(emit_q.descriptions)
        return gf, dict(raw=raw, wave=wave, fwhm_nm=fwhm, band_index=indexes,
                        qa=flags, qa_names=qa_names, grade=reference["grade"]), mask, fractions

    def polygon_bag(self, geometry, crs):
        if self.config["role"] != "lithology":
            raise ValueError("Polygon bags require the lithology input configuration")
        name = self.config["default_product"]
        with rasterio.open(self.path(self.product(name)["path"])) as ds:
            geom = transform_geom(crs, ds.crs, geometry)
            try:
                win = geometry_window(ds, [geom])
            except WindowError:
                return None
            width = ds.width
        spectra = self.read(name, win)
        inside = geometry_mask([geom], spectra.support.shape, spectra.transform, invert=True)
        included = inside & spectra.support & np.isfinite(spectra.data).any(axis=0)
        rows, cols = np.where(included)
        if rows.size == 0:
            return None
        values = spectra.data[:, rows, cols].T
        x, y = spectra.transform * (cols + .5, rows + .5)
        bag = spectra.payload()
        band_counts = np.isfinite(values).sum(axis=0)
        usable = band_counts > 0
        quantiles = np.full((3, len(band_counts)), np.nan, "float32")
        mean = np.full(len(band_counts), np.nan, "float32")
        std = mean.copy()
        if usable.any():
            quantiles[:, usable] = np.nanquantile(values[:, usable], [.1, .5, .9], axis=0)
            mean[usable] = np.nanmean(values[:, usable], axis=0)
            std[usable] = np.nanstd(values[:, usable], axis=0)
        bag.update(raw=np.where(included[None], spectra.data, np.nan), mask=included,
                   pixel_spectra=values, parent_ids=(rows + int(win.row_off)) * width + cols + int(win.col_off),
                   pixel_x=x, pixel_y=y, band_counts=band_counts, band_quantiles=quantiles,
                   band_mean=mean, band_std=std, quantile_levels=np.array([.1, .5, .9]),
                   centers_in_polygon=int(inside.sum()))
        return bag


def cli(config):
    parser = argparse.ArgumentParser(description="Local GF-5B data access")
    parser.add_argument("--config", type=Path, default=config)
    sub = parser.add_subparsers(dest="operation", required=True)
    sub.add_parser("info")
    for name in ("read", "compare"):
        cmd = sub.add_parser(name)
        cmd.add_argument("--window", nargs=4, type=int, metavar=("COL", "ROW", "WIDTH", "HEIGHT"))
        cmd.add_argument("--diagnostic")
        cmd.add_argument("--out", type=Path)
        if name == "read":
            cmd.add_argument("--product")
            cmd.add_argument("--all-bands", action="store_true")
            cmd.add_argument("--spectral-role", choices=("unassigned", "scene_endmember_candidate",
                                                         "mixed_pixel_candidate"), default="unassigned")
    polygon = sub.add_parser("polygon")
    polygon.add_argument("vector", type=Path)
    polygon.add_argument("--layer")
    polygon.add_argument("--index", type=int, default=0, help="Zero-based feature position")
    polygon.add_argument("--out", type=Path)
    args = parser.parse_args()
    asset = SpectralAsset(args.config)
    if args.operation == "info":
        print(json.dumps(asset.describe(), ensure_ascii=False, indent=2))
        return
    if asset.config["role"] == "structure":
        parser.error("The structure configuration provides coverage and catalog metadata through info")
    if args.operation == "polygon":
        with fiona.open(args.vector, layer=args.layer) as source:
            geom = next((feat.geometry for i, feat in enumerate(source) if i == args.index), None)
            if geom is None:
                parser.error("Feature has no geometry or its position is unavailable")
            bag = asset.polygon_bag(geom.__geo_interface__, source.crs)
        if bag is None:
            print(json.dumps({"status": "no_observed_pixel_centers"}))
            return
        payload = bag
        report = {"grade": bag["grade"], "pixels": len(bag["parent_ids"]), "bands": len(bag["wave"]),
                  "centers_in_polygon": bag["centers_in_polygon"],
                  "min_band_support": int(bag["band_counts"].min())}
    else:
        window = Window(*args.window) if args.window else None
        if args.operation == "compare":
            gf, emit, common, fractions = asset.read_emit_pair(window, args.diagnostic)
            payload = {f"gf5b_{k}": v for k, v in gf.payload().items()}
            payload.update({f"emit_{k}": v for k, v in emit.items()})
            payload.update(common_mask=common, gf5b_valid_fraction=fractions,
                           co_registration_status=asset.catalog["comparison"]["co_registration_status"])
            report = {"gf5b_shape": gf.data.shape, "emit_shape": emit["raw"].shape,
                      "common_pixels": int(common.sum()), **asset.catalog["comparison"]}
        else:
            spectra = asset.read(args.product, window, args.diagnostic, False if args.all_bands else None)
            payload = dict(spectra.payload(), spectral_role=args.spectral_role)
            report = {"shape": spectra.data.shape, "grade": spectra.grade, "crs": spectra.crs,
                      "observed_pixels": int(spectra.support.sum()), "spectral_role": args.spectral_role}
    if args.out:
        with args.out.open("xb") as stream:
            np.savez_compressed(stream, **payload)
        report["output"] = str(args.out.resolve())
    print(json.dumps(report, ensure_ascii=False))
