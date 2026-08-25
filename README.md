# ocular

`ocular` provides tools for obtaining high quality spectral indices from
reusable remote sensing pipelines.

## Installation

```R
o <- options() # store original options

options(pkg.build_vignettes = TRUE)

if (!require("pak")) {
  install.packages("pak")
}

pak::pak("AAGI-AUS/ocular")
options(o) # reset options
```

Given a point and a date range, `ocular` provides spectral indices commonly
used in crop monitoring and agricultural research, enables batched retrieval
of vegetative indices (VI) values for multi-site datasets, and offers reusable
pipelines for bundling pre-processing, sample-free field boundary delineation,
and other applied steps.

## Retrieval and processing

Landsat and Sentinel-2 satellite imagery are retrieved from the Microsoft
Planetary Computer Data Catalog and stored as SpatioTemporal Asset Catalog
(STAC) scenes in an `ocular` object. This allows the quality and temporal
coverage of the retrieved data to be optimised before it is transformed to a
data frame or raster.

Cloud masking and cloud cover filtering are set during retrieval, while other
processing steps and output stages are configured and applied in reusable
pipelines. Sample-free field boundary delineation in `ocular` consists of four
modular functions. In custom pipelines, these functions can be reordered or
reconfigured. For a default sequence, use `boundary_delineation()`.

Experimental features include:

- Landsat-MODIS data fusion, which can increase temporal coverage
  (leave-one-out validation against held-out Landsat scenes is available to
  assess the viability of the fused estimates).
- Fields of The World (FTW) global field boundary data incorporated as an
  optional prior to guide field boundary delineation, provided they are
  supplied in a compatible GeoParquet file.

## Sample-free field boundary delineation

`ocular` field boundary delineation uses a multi-temporal spectral feature
stack. It minimises the need for labelled training data by using the supplied
point as a reference (seed). The four modular functions are:

- `segment_area()`: **initial search to classify a plausible field area** using
  a constrained breadth-first search flood fill.
- `trace_perimeter()`: **perimeter tracing** for refining irregular or noisy
  field boundaries.
- `segment_interior()`: **within-field segmentation** to remove areas whose VI
  values fall outside the user-specified bounds or differ from the seed
  signature for the targeted crop production zone.
- `split_area()`: **field splitting** applicable where a targeted zone is
  surrounded by other cropping areas, or multiple adjacent zones are required
  in a single mask.

Placing these functions at different stages produces different results, in
addition to modifying the parameter settings. `ocular` pipelines offer a
convenient way to test and apply post-retrieval steps.

Global field boundary data from Fields of The World (FTW) can be incorporated
as an optional prior to guide delineation, provided they are supplied in a
compatible field boundary GeoParquet file. Whether this improves accuracy
requires independent validation. As such, it is considered an
**experimental** feature.

## Landsat-MODIS data fusion

Landsat-MODIS data fusion uses daily MCD43A4 NBAR data to estimate additional
VI values not covered by Landsat. Leave-one-out validation is available
against held-out Landsat scenes. This provides reportable metrics for assessing
the viability of fused estimates for your workflow. However, these metrics do
not establish general accuracy or replace independent validation.
Landsat-MODIS data fusion is therefore also marked as **experimental**.
