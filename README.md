# MODIS MAIAC smoke AOD over Canada

Monthly, 25 km aerosol optical depth and smoke-condition statistics over Canada from MODIS
MAIAC (`MCD19A2.061`), 2000-02 → 2025-12.

> **Attribution & Provenance:**  
> Adapted from **Andrew Dessler's** MODIS MAIAC analysis for CONUS (https://github.com/aedessler/smoke).  
> Adapted by **Mitchell Dickau** with **Gemini** to analyze smoke days, aerosol optical depth, and wildfire smoke frequency across **Canada**.

---

## Layout

```
maiac/        the pipeline — raw NASA granules on a GCE spot VM
              (see maiac/README.md for everything operational)
figures/      annual and monthly smoke frequency maps
```

Commands throughout assume the repository root as the working directory
(`python3 maiac/plot_month.py data/maiac/monthly/…`).

**Not in this repository.** The paths below are referenced throughout because that is where the
code writes and reads them, but they are outputs and working notes rather than source, so they
are not archived here. Running the pipeline recreates `data/`; the rest is local:

```
data/maiac_smoke_25km_monthly.nc    the combined archive — the file you want
data/maiac/monthly/                 311 per-month NetCDF checkpoints (2000-02 -> 2025-12)
data/maiac/manifest.csv             per-month run records
```

## Two implementations, one science question

The **current** product is `maiac/`: it downloads raw `MCD19A2.061` HDF granules from NASA over
authenticated HTTPS and aggregates them on a Google Compute Engine spot VM. No Earth Engine.
Full archive run of record — 311/311 months, 0 failures — is
documented in [`maiac/README.md`](maiac/README.md).

## The archive file

`data/maiac_smoke_25km_monthly.nc` — dims `time × y × x` = 311 × 188 × 219, EPSG:3978 Canada Atlas
Lambert Equal Area, 25 km, monthly 2000-02 → 2025-12. Built by `maiac/concat_archive.py` from the monthly
checkpoints, which re-checks that every month shares one grid before combining.

All six variables are built from **four accumulators** per 25 km cell per month, which is why
they are so tightly related:

```
one observation      →  QA screen: land, best quality (AOD_QA bits 8-11 == 0), AOD >= 0
one native pixel-day →  A = mean AOD over that day's valid obs (Terra+Aqua averaged)
   (1 km, 1 day)        B = 1 if any valid obs
                        C = mean AOD over obs where the smoke model was selected
                        D = 1 if any smoke obs
one 25 km cell/month →  sum A, B, C, D over all pixels and all days
```

The file stores `ΣB` and `ΣD` directly plus three ratios. Medians and NaN fractions below are
measured over the archive as shipped:

`crs` is not data. It is a CF grid-mapping stub carrying the projection (Canada Atlas Lambert
Equal Area, standard parallels 49.0/77.0, central meridian −95°, latitude of origin 49.0°, EPSG:3978);
every other variable points at it via `grid_mapping`. `units = 1` throughout means dimensionless — AOD
genuinely has no units.

Two identities hold in the file (verified, not just intended):

```
smoke_aod_index          = smoke_pixel_day_fraction × mean_smoke_aod_055   (to float32 noise)
smoke_pixel_day_fraction = smoke_pixel_day_weight / valid_pixel_day_weight (exact)
```

The index factorizes cleanly into occurrence × intensity. That is the whole reason all three
smoke variables exist rather than one.

### Three things to know before using it

**1. `mean_smoke_aod_055` is not smoke-only AOD.** MAIAC does not retrieve separate smoke and
background optical depths. The `AOD_QA` aerosol-model bits say which aerosol model the
retrieval *selected*, nothing more, and `Optical_Depth_055` remains total-column AOD under that
model. Reporting it as the smoke contribution to total AOD would be wrong.

The two companion fields split the signal: `smoke_pixel_day_fraction` isolates **occurrence**;
`smoke_aod_index` combines occurrence and magnitude, since non-smoke valid days contribute zero
to the numerator but still count in the denominator. For "how smoky was this month," the index
is usually what you want — the conditional mean is conditioned on an event whose frequency is
itself the signal.

It is also **NaN exactly where the answer is "no smoke."** Cells where `ΣD = 0` leave the
conditional mean undefined. Average it over time or space and you silently drop every clean
case, biasing the result high. `smoke_aod_index` is **0** in those same cells, which is the truthful value.

(The ~58 % NaN floor shared by the other ratios is the Canada land mask plus non-retrieval cells:
17,140 of the 41,172 cells in the bounding grid rectangle fall inside Canada's land boundary.)

**2. Aggregation is a ratio of sums, not a mean of means.** Screen at native 1 km → collapse to
pixel-days → sum through the month → bin to 25 km → *then* form ratios. The same logic applies
across **time**, which is why the weights are in the file:

```python
annual = ds["mean_aod_055"] * ds["valid_pixel_day_weight"]   # weight before averaging
annual = annual.sum("time") / ds["valid_pixel_day_weight"].sum("time")   # right
annual = ds["mean_aod_055"].mean("time")                                 # biased
```

The bias is worst exactly where it matters most — winter and the cloudy Pacific Northwest,
where valid-day counts swing hardest month to month.

**3. Read every ratio next to its weight.** A cell with 3 valid pixel-days and one with 23,000
both produce a `mean_aod_055`, and only one of them means anything:

```python
ds = ds.where(ds["valid_pixel_day_weight"] >= 5)   # reasonable first pass; vary it
```

Coverage runs 28.6–41.5 % of the grid rectangle; the low end is winter, where snow sets the QA
bits away from "land" and those pixels are correctly excluded rather than retrieved through
snow. `valid_pixel_day_weight` also **steps up in mid-2002** when Aqua joins Terra, roughly
doubling observations per pixel-day — a real feature of the record, but trends spanning that
boundary need care.


## Requirements

Reading the archive locally needs only `xarray`, `netCDF4`, `numpy`, `matplotlib` (and
`rioxarray` + `requests` for `maiac/plot_month.py`). The pipeline itself needs `earthaccess`,
`google-cloud-storage` and an HDF4-capable GDAL, but only on the VM —
`maiac/startup_script.sh` installs them there. Nothing in `maiac/` needs to run locally, and
its 44 tests run without credentials, network or GDAL.
