# ocular 0.1.2

This update fixes retrieval, delineation and fusion errors. Some values, field
masks and included dates change from 0.1.1.

## Retrieval and outputs

- Remove the Sentinel-2 digital-number offset for processing baseline 04.00
  and later before scaling bands and calculating indices.
- For duplicate Sentinel-2 products from the same acquisition and tile, select
  the highest available processing baseline by default. Use
  `rs_params(s2_processing = "lowest_baseline")` to select the lowest available
  baseline. Unreadable metadata is reported, and ties are deterministic.
- Include the whole `end_date` in retrievals. Use MCD43A4 product dates and
  retain only dates within the requested period.
- Align detection scenes on a common grid when their geometry differs.
- Support single-scene raster output. Direct raster and time-series retrievals
  avoid unnecessary detection assets and feature preparation.
- Use the same analysis area for undelineated rasters and time series. Keep
  monthly and yearly rows in calendar order when fused scenes are present.
- Correct fusion band stacking and exclude estimates or validation holdouts
  without a sufficiently close MCD43A4 observation at the required date.

## Field boundary delineation

- Defaults not explicitly selected in `rs_params()` no longer replace stored
  or calibrated settings. Forward explicit settings when splitting or interior
  segmentation first needs to run area growth.
- Prevent windows with insufficient measurable support from forcing other
  calibration windows to expand. Extend the final fill window to the full
  grid and support grids one row or one column wide.
- Rerun delineation on the current mask with fresh stage state, retaining
  feature data and settings. All stages remain active.
- Retain eligible calibrated signatures when learned signatures are missing.
  This prevents later interior errors after calibrated fills and can change
  merging and subsequent boundary refinement.
- For named area raster stacks, clean each area, union the cleaned masks and
  apply the combined mask to the original value layers.
- Retain a usable FTW prior when an unrelated polygon cannot be tested,
  provided the unusable polygon cannot contain the supplied point.
- Avoid repeating an identical failed calibration refresh within one
  invocation. Later calls and changed calibration inputs still retry.

## Pipelines, batches and caches

- Reusable pipelines use the supplied input, capture positional settings
  correctly and run when ocular is loaded without being attached.
- Retry transient STAC failures in batch rounds and return completed results
  when interrupted. No-scene results and pipeline failures are final.
- Add an `outcomes` attribute to batch results, with one row per input site.
  It distinguishes success, empty success, failure and interruption, and records
  returned date counts, retrieval attempts and cache use. Cached progress
  messages also report date counts.
- Give equivalent numeric settings the same cache key. Include retrieval
  settings, retry policy and package version in the key. Batch results
  cached under an earlier package version are retrieved and processed again.
- Detect changes to local FTW source files in cache keys. Do not cache results
  with no source year when a date column was requested.

## Documentation and diagnostics

- Clarify errors for invalid settings, missing feature windows, unusable FTW
  references and unsupported fusion-validation input. Report dropped
  acquisitions and preserve the cause of failed catalogue queries.
- Revise the reference documentation, README and vignettes, and make code
  layout consistent. Reuse repeated delineation and fusion calculations.

FTW priors and Landsat-MODIS fusion remain experimental.

# ocular 0.1.1

- Correct raster assembly when band lists or scene lists have names. Preserve
  band identifiers during stacking and merge scenes across acquisition dates.

# ocular 0.1.0

## Parameter handling

- Explicit settings supplied through `rs_params()`, including values equal to
  the defaults, take precedence over calibration baselines.
- An omitted `interior_sensitivity` or `split_sensitivity` inherits the stored
  setting; an explicit `NULL` disables the corresponding stage.
- Removed the pre-publication categorical calibration extent and optional
  circular growth constraint. Self-calibration now expands a complete,
  seed-centred square until its precision condition is met or the feature grid
  prevents further symmetric expansion. Initial area growth has no hidden
  radial constraint.

## Runtime behaviour

- Field splitting now defaults to two consecutive measurable divergent
  line-of-sight profiles and evaluates each axis provisionally. An axis is
  rolled back if it would remove any cell in the finite, 8-connected
  seed-signature component containing the supplied point. Prepared ordinary
  synthetic regressions cover this software invariant, not field identity or
  delineation accuracy.
- Successful field boundary delineation stages are quiet by default rather
  than printing internal search diagnostics. `segment_area()` warns once when
  no candidate field is identified and returns an empty mask.
- FTW queries warn when an unavailable local `icu` extension prevents a source
  year from being read; the field boundary can still provide undated
  calibration support.
- Removed unused tuning and visual-cleanup controls, obsolete aliases and
  convenience wrappers, and an unreachable MOD13Q1 retrieval route before the
  first public release.

## Documentation

- Revised the package and API reference documentation to distinguish initial
  area growth, perimeter refinement, within-field segmentation, and field
  splitting, and to state the validation limits of experimental FTW prior and
  Landsat-MODIS fusion features.
- Reorganised the vignettes around three user tasks: a first retrieval and
  output workflow, field boundary delineation with reusable multi-site
  pipelines, and experimental features with their validation requirements.

# ocular 0.0.93

## Rename

- Renamed the package from `oculaR` to `ocular`.
- New objects use `ocular` and `ocular_pipeline` as their primary S3 classes.

# oculaR 0.0.92

## User-facing changes

- Added explicit `source =` selection to `get_rs()`. Specify `"landsat-8"` or
  `"landsat-5"` when a Landsat workflow is required, including Landsat-MODIS
  fusion and fusion validation.
- Standardised the `search_windows`, `vi_sensitivity`, and `vi_threshold`
  parameter names.
- Added `is_fused` to identify synthetic fused scenes. Fusion diagnostics
  report `anchor_mae` and `anchor_nmae` separately from held-out fusion error.
- Added optional field boundary priors from compatible GeoParquet files through
  `add_ftw_prior()` and `diagnose_against_ftw()`.
- Added resumable, cached multi-site processing with `batch_rs()`.

## Reliability and validation

- Strengthened validation of dates, geometry, sources, masks, aggregation
  functions, and delineation parameters; invalid inputs now fail early with
  clearer messages.
- Improved calibration, cache-key construction, and cache-clearing safeguards.
  Custom batch pipelines require an explicit cache key before their results are
  cached.
- Added offline tests across the public API and updated package documentation.
- Live STAC and remote field boundary GeoParquet queries, as well as scientific
  accuracy against independent reference data, remain separate validation
  tasks.

# oculaR 0.0.91

## Field-boundary reference support

- Added `diagnose_against_ftw()` for pixel-level agreement statistics (IoU,
  F1, precision, and recall) between a delineation and a reference polygon.
- Added the ability to use an attached field polygon as calibration support.
  Agreement with the same polygon is an in-sample diagnostic, not an
  independent accuracy assessment.

# oculaR 0.0.90

- Initial feature development: STAC retrieval, sample-free field boundary
  delineation, optional coarse-fine fusion, and leave-one-out fusion
  evaluation.
