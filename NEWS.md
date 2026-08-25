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
