"""Plot local spectral QA from the finished GF-5B asset run."""
from pathlib import Path
import argparse
import csv
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import numpy as np
import rasterio


def table(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def vector(rows, name):
    return np.array([float(row[name]) if row[name] else np.nan for row in rows])


def plot(assets):
    rows = table(assets / "GF5B_wavelength_and_valid_bands.csv")
    sites = table(assets / "GF5B_representative_spectra.csv")
    meta = json.loads((assets / "GF5B_spatial_spectral_summary.json").read_text(encoding="utf-8"))
    wave = vector(rows, "wavelength_nm")
    accepted = np.array([r["effective_valid"] == "True" for r in rows])
    blue, orange, gray = "#216680", "#a54c25", "#a3a5aa"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    fig, axs = plt.subplots(4, 2, figsize=(11, 10), layout="constrained")
    ax = axs.ravel()

    ax[0].scatter(np.arange(1, len(rows) + 1)[~accepted], wave[~accepted], s=6, c=gray, label="Screened")
    ax[0].scatter(np.arange(1, len(rows) + 1)[accepted], wave[accepted], s=6, c=blue, label="Retained")
    ax[0].set(xlabel="Radiance band (1-based)", ylabel="Center wavelength (nm)")
    ax[0].legend(frameon=False)

    for camera, color in (("VNIR", blue), ("SWIR", orange)):
        take = np.array([r["segment"] == camera for r in rows])
        ax[1].plot(wave[take], vector(rows, "aoi_mean_radiance")[take], color=color, label=camera)
    ax[1].set(xlabel="Center wavelength (nm)", ylabel="AOI mean radiance")
    ax[1].legend(frameon=False)

    noise = vector(rows, "spatial_residual_to_median") * 100
    ax[2].scatter(wave[~accepted], noise[~accepted], s=6, c=gray)
    ax[2].scatter(wave[accepted], noise[accepted], s=6, c=blue)
    ax[2].set(xlabel="Center wavelength (nm)", ylabel="Spatial residual / median (%)", yscale="log")

    stripe = vector(rows, "detector_stripe_column_p95_fraction") * 100
    ax[3].scatter(wave[~accepted], stripe[~accepted], s=6, c=gray)
    ax[3].scatter(wave[accepted], stripe[accepted], s=6, c=orange)
    ax[3].set(xlabel="Center wavelength (nm)", ylabel="Column residual P95 / median DN (%)")

    n = meta["native"]["observed"]
    ax[4].plot(wave[accepted], vector(rows, "aoi_valid_analysis_pixels")[accepted] / n * 100, ".", c=blue)
    ax[4].set(xlabel="Center wavelength (nm)", ylabel="Finite QUAC values in AOI (%)")

    for name, color in zip(sorted({r["site"] for r in sites}),
                           ("#004c6d", "#478268", "#933e57", "#a6591d", "#695c9d")):
        records = [r for r in sites if r["site"] == name]
        w, z = vector(records, "wavelength_nm"), vector(records, "value")
        gaps = np.flatnonzero(np.diff(w) > 14) + 1
        for i, indexes in enumerate(np.split(np.arange(len(w)), gaps)):
            ax[5].plot(w[indexes], z[indexes], lw=.9, color=color, label=name if i == 0 else None)
        valid = (w >= 2080) & (w <= 2450)
        ax[6].plot(w[valid], z[valid], lw=.9, color=color)
    ax[5].legend(frameon=False, ncol=3)
    ax[5].set(xlabel="Center wavelength (nm)", ylabel="QUAC estimated reflectance")
    for lo, hi in ((2155, 2175), (2190, 2220), (2240, 2270), (2310, 2360)):
        ax[6].axvspan(lo, hi, color="#dce4e8", zorder=0)
    ax[6].set(xlabel="SWIR wavelength (nm)", ylabel="QUAC estimated reflectance")

    with rasterio.open(assets / "GF5B_native30m_AOI_QA.tif") as ds:
        flags = ds.read(1)
    status = np.zeros(flags.shape, dtype="uint8")
    status[(flags & 1) != 0] = 1
    status[(flags & 2) != 0] = 2
    status[(flags & 4) != 0] = 3
    colors = ["#f3f2f1", "#98999e", "#d68d54", "#438393"]
    ax[7].imshow(status, cmap=ListedColormap(colors), norm=BoundaryNorm(np.arange(-.5, 4), 4),
                 interpolation="nearest")
    ax[7].set(xlabel="UTM 30 m grid column", ylabel="UTM 30 m grid row")
    ax[7].legend(handles=[Patch(facecolor=colors[i], label=label) for i, label in
                 ((0, "Outside AOI"), (2, "Partial band support"), (3, "All retained bands finite"))],
                 frameon=True, loc="lower left", fontsize=7)

    titles = ("Wavelength and band eligibility", "Full-spectrum mean radiance", "Local noise proxy",
              "Native detector striping proxy", "Per-band AOI support", "Five AOI spectra",
              "SWIR diagnostic regions", "AOI and spectral support")
    for i, (panel, title) in enumerate(zip(ax, titles)):
        panel.set_title(f"{chr(97 + i)}  {title}", loc="left", fontweight="bold", fontsize=9)
    fig.savefig(assets / "GF5B_spectral_and_spatial_QA.png", dpi=300)
    fig.savefig(assets / "GF5B_spectral_and_spatial_QA.pdf")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("assets", type=Path)
    plot(parser.parse_args().assets)
