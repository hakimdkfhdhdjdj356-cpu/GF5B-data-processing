"""Prepare the local ellipsoidal DSM for the bounded AHSI RPC run."""
import argparse
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject
from rasterio.windows import from_bounds
from rasterio.warp import transform_bounds

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ahsi.config import aux_path, load_config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    cfg = load_config(parser.parse_args().config)
    root = Path(cfg["project_root"])
    paths = cfg["scene_assets"]
    output = aux_path(cfg, "rpc_dem")
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(output)

    with rasterio.open(root / paths["dem"]) as dem, rasterio.open(root / paths["geoid"]) as geoid:
        if dem.descriptions[0] != "DEM" or "EGM2008" not in dem.tags().get("meaning", ""):
            raise ValueError("Expected the local EGM2008 DSM in band 1")
        bounds = transform_bounds(dem.crs, "EPSG:4326", *dem.bounds)
        grid_window = from_bounds(*bounds, transform=geoid.transform)
        grid_window = rasterio.windows.Window(
            np.floor(grid_window.col_off) - 2, np.floor(grid_window.row_off) - 2,
            np.ceil(grid_window.width) + 5, np.ceil(grid_window.height) + 5)
        undulation = np.full(dem.shape, np.nan, dtype="float32")
        reproject(geoid.read(1, window=grid_window), undulation,
                  src_crs="EPSG:4326", src_transform=geoid.window_transform(grid_window),
                  dst_crs=dem.crs, dst_transform=dem.transform, dst_nodata=np.nan,
                  resampling=Resampling.bilinear)
        heights = dem.read(1, masked=True).filled(np.nan)
        if not np.isfinite(undulation[np.isfinite(heights)]).all():
            raise ValueError("Geoid grid does not cover valid DSM pixels")
        h = heights + undulation
        profile = dem.profile.copy()
        profile.update(count=1, dtype="float32", nodata=-9999.0, compress="deflate")
        with rasterio.open(output, "w", **profile) as dst:
            dst.write(np.where(np.isfinite(h), h, -9999.0), 1)
            dst.update_tags(vertical_reference="ellipsoid; h = H(EGM2008) + N",
                            dsm_source=str(dem.name), geoid_source=str(geoid.name),
                            undulation_min_m=float(np.nanmin(undulation)),
                            undulation_max_m=float(np.nanmax(undulation)))
    print(output)


if __name__ == "__main__":
    main()
