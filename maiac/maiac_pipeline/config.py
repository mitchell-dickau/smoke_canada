"""Frozen constants for the MAIAC 25 km Canada pipeline.

Everything here is fixed for the life of the project. The target grid in
particular must never change: every monthly file has to share one CRS,
transform, and set of x/y coordinates so the months concatenate along `time`
without alignment games.
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------- product ---
SHORT_NAME = "MCD19A2"
VERSION = "061"

# Deliberately larger than Canada. The exact land mask is applied at native
# 1 km resolution later, so this only has to be generous enough not to clip.
CANADA_BBOX = (-141.0, 41.5, -52.0, 83.5)
CONUS_BBOX = CANADA_BBOX  # Backward compatibility alias

# AOD_QA bit fields (MCD19A2 user guide).
#   bits 3-4    land/water/snow/ice   0 = land
#   bits 8-11   AOD quality           0 = best quality, 11 = high-but-reliable
#   bits 13-14  aerosol model         1 = smoke
QA_SURFACE_SHIFT, QA_SURFACE_MASK = 3, 0b11
QA_QUALITY_SHIFT, QA_QUALITY_MASK = 8, 0b1111
QA_MODEL_SHIFT, QA_MODEL_MASK = 13, 0b11

SURFACE_LAND = 0
AEROSOL_MODEL_SMOKE = 1

# Primary product keeps only best-quality retrievals.
QUALITY_PRIMARY = (0,)
QUALITY_PERMISSIVE = (0, 11)

AOD_SCALE = 0.001
SUBDATASET_AOD = ":grid1km:Optical_Depth_055"
SUBDATASET_QA = ":grid1km:AOD_QA"

# ------------------------------------------------------------ target grid ---
# EPSG:3978 - NAD83 / Canada Atlas Lambert Equal Area
TARGET_CRS = "EPSG:3978"
TARGET_RES = 25_000.0
TARGET_XMIN = -2_400_000.0
TARGET_XMAX = 3_075_000.0
TARGET_YMIN = -800_000.0
TARGET_YMAX = 3_900_000.0

TARGET_NX = int(round((TARGET_XMAX - TARGET_XMIN) / TARGET_RES))  # 219
TARGET_NY = int(round((TARGET_YMAX - TARGET_YMIN) / TARGET_RES))  # 188
TARGET_NCELLS = TARGET_NX * TARGET_NY


def target_coords() -> tuple[np.ndarray, np.ndarray]:
    """Cell-centre x and y coordinate vectors. y descends, north-up raster."""
    x = TARGET_XMIN + (np.arange(TARGET_NX) + 0.5) * TARGET_RES
    y = TARGET_YMAX - (np.arange(TARGET_NY) + 0.5) * TARGET_RES
    return x.astype("float64"), y.astype("float64")


# ------------------------------------------- MODIS sinusoidal tile geometry --
# The sphere MODIS land products are gridded on, not WGS84.
MODIS_SINU_PROJ = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"
MODIS_TILE_SIZE = 1_111_950.5196666666  # metres, 10 degrees at the equator
MODIS_X_MIN = -20_015_109.354
MODIS_Y_MAX = 10_007_554.677
MODIS_TILE_PIXELS = 1200
MODIS_PIXEL_SIZE = MODIS_TILE_SIZE / MODIS_TILE_PIXELS  # 926.625433 m


def tile_pixel_centres(h: int, v: int) -> tuple[np.ndarray, np.ndarray]:
    """Sinusoidal x/y of every 1 km pixel centre in tile hHHvVV.

    Returns two (1200, 1200) float64 arrays.
    """
    x0 = MODIS_X_MIN + h * MODIS_TILE_SIZE
    y0 = MODIS_Y_MAX - v * MODIS_TILE_SIZE
    cols = x0 + (np.arange(MODIS_TILE_PIXELS) + 0.5) * MODIS_PIXEL_SIZE
    rows = y0 - (np.arange(MODIS_TILE_PIXELS) + 0.5) * MODIS_PIXEL_SIZE
    return np.meshgrid(cols, rows)


# --------------------------------------------------------------- run range ---
# MCD19A2.061 begins 2000-02-24 (Terra only). Aqua joins 2002-07, so months
# before that are single-platform and have roughly half the observations per
# pixel-day. That is a real feature of the record, not a defect -- but it is
# why `valid_pixel_day_weight` is written to every file.
ARCHIVE_START = "2000-02"

# Canada boundary polygon source (Natural Earth 10m Admin-1 Provinces/Territories)
CANADA_ADMIN1_URL = "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_admin_1_states_provinces.geojson"
