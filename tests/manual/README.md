# Maintainer release checks

These are maintainer-run checks. They are not sourced by `testthat`,
`R CMD check`, or the vignettes. The live script installs nothing, does not
install DuckDB extensions, and makes no query unless
`OCULAR_RUN_LIVE` is set to `true` (case-insensitive, with surrounding
whitespace ignored). Values such as `1` and `yes` do not enable it.

None of the commands in this file was run when these instructions were
prepared. Run them from a fresh terminal, in order, and inspect the stated
outputs before continuing.

The complete offline sequence in sections 1--3 and 6 is automated by:

```sh
/bin/sh ocular/tests/manual/verify-offline-release.sh
```

Run that command from the workspace root. It regenerates documentation, runs
the focused and full test suites, builds and checks one checksum-tracked
tarball, verifies an isolated installation, and creates one new split design
diagnostic. It never enables the live cases, `--as-cran`, repository checks,
or publication. The individual commands remain below for inspection and
targeted reruns.

## 1. Establish a workspace-local check environment

Start at the workspace root: the directory that contains `AGENTS.md` and the
`ocular/` package directory. Confirm both files before setting any paths:

```sh
pwd
test -f AGENTS.md || exit 1
test -f ocular/DESCRIPTION || exit 1
test ! -L ocular || exit 1

workspace_root="$(pwd -P)"
package_root="$workspace_root/ocular"

ensure_workspace_dir() (
  target="$1"
  root="$2"
  if [ -L "$target" ]; then
    echo "Refusing symlinked directory: $target" >&2
    exit 1
  fi
  if [ -e "$target" ] && [ ! -d "$target" ]; then
    echo "Refusing non-directory path: $target" >&2
    exit 1
  fi
  if [ ! -d "$target" ]; then
    mkdir "$target" || exit 1
  fi
  resolved="$(cd "$target" && pwd -P)" || exit 1
  case "$resolved" in
    "$root"|"$root"/*) ;;
    *) echo "Directory resolves outside package: $target" >&2; exit 1 ;;
  esac
)
```

Create and verify each directory one level at a time. This prevents an
existing symlink from redirecting temporary files, caches, check output, or an
installed library outside the package:

```sh
ensure_workspace_dir "$package_root/.check-artifacts" "$package_root" || exit 1
ensure_workspace_dir "$package_root/.check-artifacts/offline" "$package_root" || exit 1
ensure_workspace_dir "$package_root/.check-artifacts/offline/tmp" "$package_root" || exit 1
ensure_workspace_dir "$package_root/.check-artifacts/offline/cache" "$package_root" || exit 1

export TMPDIR="$package_root/.check-artifacts/offline/tmp"
export R_USER_CACHE_DIR="$package_root/.check-artifacts/offline/cache"
```

The supplied R profile disables configured package repositories. These
selected checks are designed to be offline: the automated tests use mocks and
the remote vignette chunks are not evaluated. The profile is not an operating
system network sandbox; site startup code or a direct HTTP/GDAL call could
still use a network connection. Disable network access at the operating-system
level as well if a hard offline guarantee is required.

Check the required R packages and Pandoc without installing anything:

```sh
R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  Rscript -e 'needed <- c("devtools", "pkgload", "rstac", "sf", "terra", "testthat", "duckdb", "DBI", "knitr", "rmarkdown"); missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) { message("Missing: ", paste(missing, collapse = ", ")); quit(status = 1) }' || exit 1

R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  Rscript -e 'if (!rmarkdown::pandoc_available()) stop("Pandoc is unavailable") else cat(rmarkdown::pandoc_version(), "\n")' || exit 1
```

Stop if either command fails. Install or update dependencies only in a
separate, deliberately network-enabled step under your own control. A final
`--as-cran` check that builds the PDF manual also needs a working LaTeX
toolchain.

## 2. Document, test, build, and check

### Regenerate roxygen output only when required

This step changes `ocular/man/`, `ocular/NAMESPACE`, and potentially
`RoxygenNote`. Run it only after roxygen source in `ocular/R/` is final:

```sh
R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  Rscript -e 'devtools::document("ocular")' || exit 1
```

Inspect every generated change. A vignette-only edit does not require roxygen
regeneration.

### Run source tests

Run the focused split regression file first:

```sh
R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  Rscript -e 'devtools::test("ocular", filter = "split-guards", stop_on_failure = TRUE)' || exit 1
```

Then run the complete source suite:

```sh
R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  Rscript -e 'devtools::test("ocular", stop_on_failure = TRUE)' || exit 1
```

An optional spelling pass, when `spelling` is already installed, is:

```sh
R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  Rscript -e 'spelling::spell_check_package("ocular")' || exit 1
```

### Build the source package

Derive the release artifact name from `DESCRIPTION`:

```sh
ocular_version="$(R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" Rscript -e 'cat(read.dcf("ocular/DESCRIPTION")[1, "Version"])')" || exit 1
ocular_tarball="ocular_${ocular_version}.tar.gz"
checked_tarball="$workspace_root/$ocular_tarball"
if [ -e "$checked_tarball" ] || [ -L "$checked_tarball" ]; then
  echo "Archive the previous $ocular_tarball inside the workspace before building." >&2
  exit 1
fi

R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  R CMD build ocular || exit 1

test -f "$checked_tarball" || exit 1
test ! -L "$checked_tarball" || exit 1
tar -tf "$checked_tarball" || exit 1

export OCULAR_CHECKED_TARBALL="$checked_tarball"
checked_tarball_md5="$(R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" Rscript -e 'path <- Sys.getenv("OCULAR_CHECKED_TARBALL"); md5 <- unname(tools::md5sum(path)); if (length(md5) != 1L || is.na(md5) || !nzchar(md5)) stop("Could not calculate tarball checksum"); cat(md5)')" || exit 1
test -n "$checked_tarball_md5" || exit 1
export OCULAR_CHECKED_TARBALL_MD5="$checked_tarball_md5"
```

Confirm that the archive contains all three Rmd vignettes and their installed
vignette outputs. Confirm that it excludes `.Rhistory`, `.Rproj.user`,
`.DS_Store`, `.check-artifacts`, `.test-artifacts`, source-tree vignette HTML
previews, `tests/manual`, and prior check directories.

### Check the built package

Use a new output parent so a prior check is not replaced:

```sh
check_parent="$package_root/.check-artifacts/offline/check-${ocular_version}-$(date +%Y%m%d-%H%M%S)-$$"
test ! -e "$check_parent" && test ! -L "$check_parent" || exit 1
ensure_workspace_dir "$check_parent" "$package_root" || exit 1

R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  _R_CHECK_CRAN_INCOMING_REMOTE_=false \
  R CMD check --output="$check_parent" --no-manual "$checked_tarball" || exit 1
```

Require zero errors, warnings, and notes. Inspect the generated
`ocular.Rcheck` directory, including its test log, examples, installed help,
and `ocular/doc/` vignette index and HTML files. Checking the tarball verifies
the release object rather than whichever source files happen to be loaded in
the working directory.

## 3. Install the checked tarball into an isolated library

Keep the validation library and R temporary directory inside the package:

```sh
live_root="$package_root/.check-artifacts/live-release-validation"
live_session="$live_root/install-${ocular_version}-$(date +%Y%m%d-%H%M%S)-$$"
live_library="$live_session/library"
live_tmp="$live_session/tmp"
live_cache="$live_session/cache"
test ! -e "$live_session" && test ! -L "$live_session" || exit 1
ensure_workspace_dir "$live_root" "$package_root" || exit 1
ensure_workspace_dir "$live_session" "$package_root" || exit 1
ensure_workspace_dir "$live_library" "$package_root" || exit 1
ensure_workspace_dir "$live_tmp" "$package_root" || exit 1
ensure_workspace_dir "$live_cache" "$package_root" || exit 1

preinstall_md5="$(R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" Rscript -e 'path <- Sys.getenv("OCULAR_CHECKED_TARBALL"); md5 <- unname(tools::md5sum(path)); if (length(md5) != 1L || is.na(md5)) stop("Could not verify tarball checksum"); cat(md5)')" || exit 1
test "$preinstall_md5" = "$OCULAR_CHECKED_TARBALL_MD5" || exit 1

R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  R CMD INSTALL --library="$live_library" "$checked_tarball" || exit 1

export OCULAR_VALIDATION_LIBRARY="$live_library"
export OCULAR_INSTALL_RECEIPT="$live_session/install-receipt.csv"
export TMPDIR="$live_tmp"
export R_USER_CACHE_DIR="$live_cache"
cd "$package_root" || exit 1

test ! -e "$OCULAR_INSTALL_RECEIPT" && \
  test ! -L "$OCULAR_INSTALL_RECEIPT" || exit 1
R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  OCULAR_WRITE_INSTALL_RECEIPT=true \
  Rscript tests/manual/write-install-receipt.R || exit 1

if [ -n "${R_LIBS:-}" ]; then
  export R_LIBS="$OCULAR_VALIDATION_LIBRARY:$R_LIBS"
else
  export R_LIBS="$OCULAR_VALIDATION_LIBRARY"
fi

R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  Rscript -e 'cat(find.package("ocular"), "\n"); cat(as.character(packageVersion("ocular")), "\n")' || exit 1
```

The printed package path must be inside `OCULAR_VALIDATION_LIBRARY` and the
version must match the local `DESCRIPTION`. The scripts enforce both
conditions, carry the post-build tarball checksum through installation, and
verify the workspace-local install receipt before running. The timestamped
install directory prevents an older same-version ocular installation from
being reused silently. `R_LIBS` places that installation first without hiding
the ordinary user and site libraries, so dependencies can still resolve from
the other active R library paths. Those paths, the receipt, and both checksums
are recorded in `install-identity.csv`. Each run creates a new timestamped
directory below
`.check-artifacts/live-release-validation/`, including its own ocular cache
and terra temporary directory. It does not clear an existing run.

The named cases are:

- `intro`: the introductory Landsat 5/7 retrieval, delineation, raster, and
  time series;
- `mcd43a4`: a short direct probe of the current internal MCD43A4 asset
  surface;
- `fusion`: the experimental-vignette Landsat 8/9 query, `add_modis()`,
  field boundary delineation, and leave-one-out agreement assessment;
- `ftw-local`: a configured local field boundary GeoParquet; and
- `ftw-remote`: a configured HTTPS or S3 field boundary GeoParquet.

With no case names, all five are requested. An empty or unknown case selection,
or missing configuration for a requested case, fails rather than silently
passing.

## 4. Run and inspect the introductory, MCD43A4, and fusion cases

The coordinates, dates, and settings are retained from the current vignettes.
These commands deliberately enable network access:

```sh
R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  OCULAR_RUN_LIVE=true \
  Rscript tests/manual/live-release-validation.R intro || exit 1

R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  OCULAR_RUN_LIVE=true \
  Rscript tests/manual/live-release-validation.R mcd43a4 || exit 1

R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
  OCULAR_RUN_LIVE=true \
  Rscript tests/manual/live-release-validation.R fusion || exit 1
```

Any warning gives that case `WARN` status and causes the run to finish with a
failure after saving its artefacts. This makes partial scene, fusion, or
held-out failures visible instead of treating them as a clean release result.
`setup-manifest.csv` records package/library setup separately;
`install-identity.csv` records the checked tarball, checksum, verified install
receipt, installed ocular path, and active library paths; `manifest.csv`
records the requested cases.

Inspect at minimum:

- `intro-delineation.png` and `intro-field-mask.tif`: the seed-centred field
  and its spatial plausibility;
- `intro-time-series.png` and `intro-series.csv`: dates, missingness, range,
  and continuity for the selected site;
- `mcd43a4-band-summary.csv`: live asset names, finite-cell counts, missing
  cells, and scaled minima/maxima;
- `fusion-time-series.png` and `fusion-series.csv`: observed Landsat dates
  versus estimates after the vignette's delineation step; and
- `fusion-validation.csv`: support counts and query-specific held-out metrics.

The MCD43A4 probe establishes that the configured asset names can be queried
and records values after the package's current range and scale handling. It
cannot independently establish that raw fill value `32767`, valid range
`0:32766`, or scale `0.0001` is scientifically correct; no authoritative
product specification is supplied with the repository. Product-specific
MCD43A4 quality layers are not currently applied.

A successful live query is not evidence of general fusion accuracy. The saved
metrics describe leave-one-out agreement for that query.

## 5. Run local and remote field boundary GeoParquet cases

No approved GeoParquet fixture or remote source is supplied with the
repository. Supply the schema, CRS, and a point known to fall inside one
field; do not infer them from column names.

The local file must resolve inside the package workspace. The remote value must
begin with `https://` or `s3://`; local paths, `file://` values, and plain HTTP
are rejected by the remote case.

Set the non-credential shared configuration. Treat site coordinates and local
paths according to your own data-sensitivity requirements:

```sh
export OCULAR_FTW_LOCAL_SOURCE="$PWD/path/to/local-fields.parquet"
export OCULAR_FTW_LONGITUDE="<known longitude inside one field>"
export OCULAR_FTW_LATITUDE="<known latitude inside one field>"
export OCULAR_FTW_START_DATE="<YYYY-MM-DD>"
export OCULAR_FTW_END_DATE="<YYYY-MM-DD>"
export OCULAR_FTW_X_METRES="1000"
export OCULAR_FTW_IMAGERY_SOURCE="landsat-8"
export OCULAR_FTW_INDEX="NDVI"
export OCULAR_FTW_CRS="4326"
export OCULAR_FTW_GEOMETRY_COL="geometry"
export OCULAR_FTW_DATETIME_COL="determination_datetime"
```

Do not put a signed URL in a command, tracked file, or shell-history entry.
Inject it with a secret manager, or read it without terminal echo:

```zsh
printf 'Remote GeoParquet URL: '
IFS= read -r -s OCULAR_FTW_REMOTE_SOURCE
printf '\n'
export OCULAR_FTW_REMOTE_SOURCE
```

The silent `read -s` form above is for zsh or Bash.

Use `OCULAR_FTW_DATETIME_COL=NULL` when the source has no date column. Set
`OCULAR_FTW_EXPECT_YEAR=true` when the check must fail unless a source year can
be extracted.

If the local and remote sources need different settings, override a shared
value with a source-specific name such as `OCULAR_FTW_LOCAL_CRS`,
`OCULAR_FTW_REMOTE_GEOMETRY_COL`, `OCULAR_FTW_LOCAL_LONGITUDE`, or
`OCULAR_FTW_REMOTE_START_DATE`. Supported suffixes are:

```text
SOURCE LONGITUDE LATITUDE START_DATE END_DATE X_METRES INDEX IMAGERY_SOURCE
CRS GEOMETRY_COL DATETIME_COL EXPECT_YEAR
```

Run the cases explicitly. The conditional preserves a non-zero R exit status
while still removing the remote source from the shell environment:

```sh
if ! R_PROFILE_USER="$package_root/tests/manual/offline-check.Rprofile" \
    OCULAR_RUN_LIVE=true \
    Rscript tests/manual/live-release-validation.R ftw-local ftw-remote
then
  unset OCULAR_FTW_REMOTE_SOURCE
  exit 1
fi

unset OCULAR_FTW_REMOTE_SOURCE
```

The `duckdb`, `DBI`, and `sf` packages must already be installed. DuckDB's
`spatial` extension must already be available; remote sources also need the
local `httpfs` extension, and dated-source validation needs `icu`. Neither the
package nor the script installs or downloads an extension.

Inspect each `ftw-*-delineation.png`, `ftw-*-series.csv`, diagnostic CSV,
`setup-manifest.csv`, `manifest.csv`, and `run.log`. The script requires the
selected polygon to contain the seed, provide calibration support, produce a
non-empty field, and produce a usable time series. It imposes no threshold for
IoU or F1. Agreement against a polygon used as a prior is in-sample and is not
an independent accuracy assessment. Configured source values are redacted from
the script's managed log and manifest.

## 6. Offline split design comparison

`split-design-comparison.R` is an opt-in, non-network visual comparison for
the selected 0.1.0 source contract, pending author-run verification. Its lobe
counts, orientations, and 3 x 3 minority size come from the historical prose;
the 40 x 40 grid, coordinates, feature values, and some settings are explicit
reconstruction assumptions because the historical test source is unavailable.
The prepared ordinary tests encode the deterministic protection and rollback
contract. Exact retained-cell counts and the reconstructed figures still
require visual author review and are not claims about field identity or
delineation accuracy.

The successful pre-correction run is preserved at
`.check-artifacts/split-design-comparison/20260825-211428-33249/`. A
post-correction run must use the launch command below and write to a new
timestamped directory; do not replace the earlier evidence.

Run it from the current `ocular/` source root. This method diagnostic loads the
current source with `pkgload::load_all()`; it does not require a built tarball,
an installed copy of `ocular`, or the release-install receipt from section 3.
The `pkgload` package and ocular's ordinary dependencies must already be
installed. The script reports a missing dependency but never installs or
downloads one.

From an R or RStudio session whose working directory is the `ocular/` source
root:

```r
source("tests/manual/split-design-comparison.R")
```

Sourcing the file is the opt-in action. The launcher starts a
fresh, offline R worker because an existing R session's `tempdir()` cannot be
relocated. The worker receives package-local temporary and cache directories;
terra also receives a package-local per-run temporary directory. The parent R
session's working directory, environment variables, repository option, and
global environment are not changed. Each run is written to a new directory
below `.check-artifacts/split-design-comparison/`. Inspect:

- `reconstructed-fixtures.csv` for the exact geometry and feature values used;
- `historical-provenance.csv` for the workspace sources, cited line ranges,
  claims used, availability, and checksums;
- `execution-identity.csv`, `source-tree-context.csv`, and
  `effective-parameters.csv` for the loaded source path and version,
  package-local worker paths, mutable source-tree context, and complete
  parameter values and explicit-ownership flags;
- `suffix-stage-schedule.csv` for the calls, overrides, and selected saved
  states in the reconstructed post-segment sequence;
- `split-results.csv`, `split-masks.rds`, and `figures/` for the direct
  gate-1/gate-2 comparison;
- the `five_split_schedule` rows and figures as a repeated-split stress test,
  not as `boundary_delineation()` output;
- the `reconstructed_post_segment_boundary_suffix` rows and figures as the
  selected cumulative states from the current default stage order
  (`cleanup_boundary = 2L`) under the indicated split settings, starting from
  the manually constructed candidate mask, not as evidence that
  `segment_area()` would produce that mask.

Seed loss, complete mask removal, and disagreement with the July retained-cell
counts remain recorded diagnostic outcomes rather than structural launcher
failures. Both gate settings now exercise the source-tree overall-median rule,
protected point-connected component, and provisional per-axis rollback. They
do not reproduce or test the July near/far-half and removed-cell-majority
rules.

No post-correction command in this section was run while the source changes
and ordinary tests were prepared.

## 7. Coverage not supplied by this targeted harness

This script is not a claim that every public workflow has received a live
integration run. It does not currently include a verified Sentinel-2/SCL
masking query or a `batch_rs()` multi-site/cache query. Automated tests cover
those interfaces with local fixtures and mocks, but a live release record
needs author-approved locations, dates, and expected visual behaviour. Add
those cases only after such examples are supplied; do not invent coordinates
or parameter settings.

## 8. Final repository-facing checks

After `AAGI-AUS/ocular` is online:

1. Verify that the `URL` and `BugReports` fields recorded in `DESCRIPTION`
   resolve to the intended repository and issue tracker.
2. Repeat the document/test/build/check sequence for the final repository
   metadata.
3. Return to the workspace root, deliberately enable network access, and run
   the following. If this is a new terminal, repeat section 1 first.

   ```sh
   cd "$workspace_root" || exit 1
   ocular_version="$(Rscript -e 'cat(read.dcf("ocular/DESCRIPTION")[1, "Version"])')" || exit 1
   ocular_tarball="ocular_${ocular_version}.tar.gz"
   checked_tarball="$workspace_root/$ocular_tarball"
   test -f "$checked_tarball" && test ! -L "$checked_tarball" || exit 1
   as_cran_root="$package_root/.check-artifacts/as-cran"
   as_cran_parent="$as_cran_root/check-${ocular_version}-$(date +%Y%m%d-%H%M%S)-$$"
   as_cran_tmp="$as_cran_root/tmp"
   as_cran_cache="$as_cran_root/cache"
   test ! -e "$as_cran_parent" && test ! -L "$as_cran_parent" || exit 1
   ensure_workspace_dir "$as_cran_root" "$package_root" || exit 1
   ensure_workspace_dir "$as_cran_parent" "$package_root" || exit 1
   ensure_workspace_dir "$as_cran_tmp" "$package_root" || exit 1
   ensure_workspace_dir "$as_cran_cache" "$package_root" || exit 1
   export TMPDIR="$as_cran_tmp"
   export R_USER_CACHE_DIR="$as_cran_cache"

   R CMD check --output="$as_cran_parent" --as-cran "$checked_tarball" || exit 1
   Rscript -e 'if (!requireNamespace("urlchecker", quietly = TRUE)) stop("urlchecker is unavailable"); urlchecker::url_check("ocular")' || exit 1
   ```

4. Inspect the installed vignette index and all three rendered articles.
5. Verify `pak::pak("AAGI-AUS/ocular")` from a fresh, isolated R library.
6. Run Linux, macOS, and Windows checks for current R and the declared minimum
   R version.

The final `--as-cran` run is intentionally network-enabled and may contact
external services. Keep its output inside the workspace or request approval
before directing it elsewhere. Repository creation does not by itself validate
the remote-data paths; retain the live artefacts with the release record.
