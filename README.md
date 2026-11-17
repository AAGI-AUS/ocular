# ocular

`ocular` retrieves Landsat and Sentinel-2 imagery and calculates spectral
indices for crop monitoring and agricultural research. Supply a point, date
range and analysis extent to produce a raster or time series, with optional
sample-free field boundary delineation. Reusable pipelines apply the same
processing across sites.

**Sample-free field boundary delineation in ocular** has not yet been
validated against independent field boundary reference data. It produces
candidate masks from spectral patterns around a supplied point. Inspect each
result against suitable imagery or reference boundaries before using it for
field-level analysis.

## Installation

Requires R 4.1 or later. Install from GitHub:

```r
if( !requireNamespace("pak", quietly = TRUE) ){
  install.packages("pak")
}
pak::pak("AAGI-AUS/ocular")
```

## Try it

This example retrieves Landsat imagery for a 1,000 m square around a published
example point, delineates a candidate field, and returns NDVI values. Retrieval
uses the Microsoft Planetary Computer STAC and requires internet access.

```r
library(ocular)

rs <- get_rs(
  longitude = -96.80342,
  latitude = 33.26425,
  start_date = "2005-03-15",
  end_date = "2005-09-05",
  x_metres = 1000,
  index_name = "NDVI",
  source = "landsat-5"
)

field <- boundary_delineation(rs)
vi_raster <- as_raster(field)
vi_series <- as_time_series(field)

terra::plot(vi_raster)
head(vi_series)
```

The raster shows median NDVI across retrieved scenes; the time series contains
`date` and `value`, with mean NDVI within the candidate field for each available
day. Inspect the delineation before using it in an analysis. An empty mask gives
an all-`NA` raster and no time-series rows.

The default cloud-cover filter selects scenes; per-pixel masking is optional.
To extract values within an existing field polygon, use the output functions'
`mask` argument instead of delineating a field.

## Guides

- [Getting started](vignettes/ocular-introduction.Rmd): inspect and save outputs,
  choose imagery and masks, and interpret missing results.
- [Field boundary delineation and reusable pipelines](vignettes/ocular-build-pipeline.Rmd):
  compare settings, refine a mask, and process multiple sites.
- [Experimental features](vignettes/ocular-experimental-features.Rmd): field
  boundary priors and Landsat-MODIS fusion, including their validation limits.

Use `help(package = "ocular")` for the function reference. Installations built
with vignettes also support `browseVignettes("ocular")`.
See [NEWS](NEWS.md) for changes in 0.1.2.
