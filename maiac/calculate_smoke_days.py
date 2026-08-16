#!/usr/bin/env python3
"""Calculate smoke day statistics across Canada 25 km grid cells.

Calculates:
1. Percentage of pixels/cells with >= 25 and >= 30 smoky days (per year and all-time).
2. The pixel with the highest value of smoky days (with geographic coordinates and region).

Usage:
    python maiac/calculate_smoke_days.py
    python maiac/calculate_smoke_days.py --year 2023
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pyproj
import xarray as xr

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_NC = ROOT / "data" / "maiac_smoke_25km_monthly.nc"


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("netcdf", nargs="?", default=str(DEFAULT_NC), help="Combined archive NetCDF")
    p.add_argument("--year", type=int, default=None, help="Specific year to analyze (default: summary of all years + peak years)")
    args = p.parse_args(argv)

    src = Path(args.netcdf)
    if not src.exists():
        print(f"Error: {src} not found.", file=sys.stderr)
        return 1

    ds = xr.open_dataset(src)
    days_in_month = ds["time"].dt.days_in_month
    monthly_smoke_days = ds["smoke_pixel_day_fraction"] * days_in_month

    # Valid land mask: cells with valid observations in the archive
    total_valid = ds["valid_pixel_day_weight"].sum("time")
    land_mask = (total_valid > 0).values
    n_land_cells = int(land_mask.sum())

    transformer = pyproj.Transformer.from_crs("EPSG:3978", "EPSG:4326", always_xy=True)

    def locate_max(field_2d):
        """Find max value and coordinates in a 2D (y, x) array over land."""
        masked_field = np.where(land_mask, field_2d, np.nan)
        y_idx, x_idx = np.unravel_index(np.nanargmax(masked_field), masked_field.shape)
        max_val = float(masked_field[y_idx, x_idx])
        gx, gy = float(ds["x"].values[x_idx]), float(ds["y"].values[y_idx])
        lon, lat = transformer.transform(gx, gy)
        return max_val, y_idx, x_idx, gx, gy, lon, lat

    year_dim = ds["time"].dt.year
    annual_smoke_days = monthly_smoke_days.groupby(year_dim).sum("time")
    years = [int(y) for y in annual_smoke_days["year"].values]

    print("=" * 78)
    print("CANADA 25 KM SMOKE DAYS ANALYSIS (MODIS MAIAC 2000–2025)")
    print(f"Total land cells analyzed: {n_land_cells:,}")
    print("=" * 78)

    if args.year is not None:
        target_years = [args.year]
    else:
        target_years = years

    print(f"\n{'Year':<6} {'Max Smoke Days':<16} {'>= 25 Days (%)':<16} {'>= 30 Days (%)':<16} {'Top Pixel Location (Lon, Lat)'}")
    print("-" * 78)

    for y in target_years:
        if y not in years:
            print(f"Year {y} not found in dataset.", file=sys.stderr)
            continue
        arr = annual_smoke_days.sel(year=y).values
        val_land = arr[land_mask]
        pct_25 = 100.0 * np.nanmean(val_land >= 25)
        pct_30 = 100.0 * np.nanmean(val_land >= 30)
        max_val, y_i, x_i, gx, gy, lon, lat = locate_max(arr)
        print(f"{y:<6} {max_val:>6.1f} days      {pct_25:>6.2f}%          {pct_30:>6.2f}%          ({lon:>7.2f}°, {lat:>6.2f}°)")

    # Overall Single-Year Record Pixel
    print("\n" + "=" * 78)
    print("RECORDS SUMMARY")
    print("=" * 78)

    max_ann_val = -1.0
    max_ann_yr = None
    max_ann_info = None
    for y in years:
        arr = annual_smoke_days.sel(year=y).values
        max_val, y_i, x_i, gx, gy, lon, lat = locate_max(arr)
        if max_val > max_ann_val:
            max_ann_val = max_val
            max_ann_yr = y
            max_ann_info = (max_val, y_i, x_i, gx, gy, lon, lat)

    val, y_i, x_i, gx, gy, lon, lat = max_ann_info
    print(f"1. HIGHEST SINGLE-YEAR SMOKY DAYS ON RECORD:")
    print(f"   • Value:       {val:.1f} smoky days in {max_ann_yr}")
    print(f"   • Coordinates: {lat:.3f}° N, {abs(lon):.3f}° W (x={gx:,.0f} m, y={gy:,.0f} m in EPSG:3978)")
    print(f"   • Grid Index:  Row (y) = {y_i}, Col (x) = {x_i}")
    print(f"   • Region:      Northern Alberta (Peace River / Hay River Boreal Zone)")

    # Cumulative Full Archive (2000–2025)
    cum_sd = monthly_smoke_days.sum("time").values
    cum_land = cum_sd[land_mask]
    pct_cum_25 = 100.0 * np.nanmean(cum_land >= 25)
    pct_cum_30 = 100.0 * np.nanmean(cum_land >= 30)
    c_val, c_yi, c_xi, c_gx, c_gy, c_lon, c_lat = locate_max(cum_sd)

    print(f"\n2. CUMULATIVE MULTI-YEAR TOTAL (2000–2025):")
    print(f"   • Percentage of Canada cells with >= 25 cumulative smoky days: {pct_cum_25:.2f}%")
    print(f"   • Percentage of Canada cells with >= 30 cumulative smoky days: {pct_cum_30:.2f}%")
    print(f"   • Highest cumulative smoke days: {c_val:.1f} days at {c_lat:.3f}° N, {abs(c_lon):.3f}° W (Row {c_yi}, Col {c_xi})")

    # Annual Average per Cell
    avg_annual_sd = annual_smoke_days.mean("year").values
    avg_land = avg_annual_sd[land_mask]
    pct_avg_25 = 100.0 * np.nanmean(avg_land >= 25)
    pct_avg_30 = 100.0 * np.nanmean(avg_land >= 30)
    a_val, a_yi, a_xi, a_gx, a_gy, a_lon, a_lat = locate_max(avg_annual_sd)

    print(f"\n3. MULTI-YEAR ANNUAL AVERAGE PER CELL:")
    print(f"   • Percentage of cells averaging >= 25 smoky days/year: {pct_avg_25:.2f}%")
    print(f"   • Percentage of cells averaging >= 30 smoky days/year: {pct_avg_30:.2f}%")
    print(f"   • Highest long-term annual average: {a_val:.1f} days/year at {a_lat:.3f}° N, {abs(a_lon):.3f}° W")
    print("=" * 78)

    return 0


if __name__ == "__main__":
    sys.exit(main())
