#!/bin/sh

# Run ocular's complete offline release verification from a clean sequence.
#
# This script does not install dependencies, enable live checks, use a remote
# repository, or publish anything. All generated files are written below
# ocular/.check-artifacts/. Run it with:
#
#   /bin/sh ocular/tests/manual/verify-offline-release.sh

set -eu

umask 077

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
package_root=$(CDPATH= cd -- "$script_dir/../.." && pwd -P)
workspace_root=$(CDPATH= cd -- "$package_root/.." && pwd -P)
profile="$package_root/tests/manual/offline-check.Rprofile"

fail() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 1
}

ensure_directory() {
  target=$1
  root=$2

  [ ! -L "$target" ] || fail "Refusing symlinked directory: $target"
  if [ -e "$target" ] && [ ! -d "$target" ]; then
    fail "A non-directory blocks verification output: $target"
  fi
  if [ ! -d "$target" ]; then
    mkdir "$target" || fail "Could not create directory: $target"
  fi
  resolved=$(CDPATH= cd -- "$target" && pwd -P) ||
    fail "Could not resolve directory: $target"
  case "$resolved" in
    "$root"|"$root"/*) ;;
    *) fail "Directory resolves outside the allowed root: $target" ;;
  esac
}

run_logged() {
  step=$1
  shift
  log="$logs_root/$step.log"
  printf '\n==> %s\n' "$step"
  if "$@" >"$log" 2>&1; then
    cat "$log"
  else
    status=$?
    cat "$log" >&2
    fail "$step failed with status $status; see $log"
  fi
}

run_logged_in() {
  step=$1
  directory=$2
  shift 2
  log="$logs_root/$step.log"
  printf '\n==> %s\n' "$step"
  if (CDPATH= cd -- "$directory" && "$@") >"$log" 2>&1; then
    cat "$log"
  else
    status=$?
    cat "$log" >&2
    fail "$step failed with status $status; see $log"
  fi
}

[ "$PWD" = "$workspace_root" ] ||
  fail "Run this script from the workspace root: $workspace_root"
[ -f "$workspace_root/AGENTS.md" ] || fail "AGENTS.md is missing."
[ -f "$package_root/DESCRIPTION" ] || fail "ocular/DESCRIPTION is missing."
[ -f "$profile" ] || fail "The offline R profile is missing."
[ ! -L "$package_root" ] || fail "The ocular package root is a symlink."
command -v R >/dev/null 2>&1 || fail "R is unavailable."
command -v Rscript >/dev/null 2>&1 || fail "Rscript is unavailable."
command -v tar >/dev/null 2>&1 || fail "tar is unavailable."

package_name=$(sed -n 's/^Package:[[:space:]]*//p' "$package_root/DESCRIPTION")
[ "$package_name" = "ocular" ] || fail "The package is not ocular."
source_version=$(sed -n 's/^Version:[[:space:]]*//p' "$package_root/DESCRIPTION")
[ -n "$source_version" ] || fail "The package version is unavailable."

artifacts_root="$package_root/.check-artifacts"
verification_root="$artifacts_root/offline-release-verification"
run_id=$(date +%Y%m%d-%H%M%S)-$$
run_root="$verification_root/$run_id"
logs_root="$run_root/logs"
tmp_root="$run_root/tmp"
cache_root="$run_root/cache"
pre_document_root="$run_root/pre-documentation"
build_root="$run_root/build"
check_root="$run_root/check"
live_release_root="$artifacts_root/live-release-validation"
install_root="$live_release_root/install-$source_version-$run_id"
install_library="$install_root/library"
install_tmp="$install_root/tmp"
install_cache="$install_root/cache"
diagnostic_launch_tmp="$run_root/diagnostic-launch-tmp"
diagnostic_launch_cache="$run_root/diagnostic-launch-cache"

ensure_directory "$artifacts_root" "$package_root"
ensure_directory "$verification_root" "$package_root"
ensure_directory "$live_release_root" "$package_root"
[ ! -e "$run_root" ] && [ ! -L "$run_root" ] ||
  fail "Refusing to reuse verification run: $run_root"
ensure_directory "$run_root" "$package_root"
for directory in \
  "$logs_root" "$tmp_root" "$cache_root" "$pre_document_root" \
  "$build_root" "$check_root" "$install_root" "$install_library" \
  "$install_tmp" "$install_cache" "$diagnostic_launch_tmp" \
  "$diagnostic_launch_cache"
do
  ensure_directory "$directory" "$package_root"
done

export TMPDIR="$tmp_root"
export R_USER_CACHE_DIR="$cache_root"
export R_PROFILE_USER="$profile"
export OCULAR_RUN_LIVE=false

printf '%s\n' "$run_root" >"$run_root/run-root.txt"

run_logged preflight-dependencies env \
  R_PROFILE_USER="$profile" TMPDIR="$tmp_root" \
  R_USER_CACHE_DIR="$cache_root" OCULAR_RUN_LIVE=false \
  Rscript -e 'needed <- c("devtools", "pkgload", "rstac", "sf", "terra", "testthat", "duckdb", "DBI", "knitr", "rmarkdown"); missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) stop("Missing installed dependencies: ", paste(missing, collapse = ", "), call. = FALSE); cat("All required packages are already installed.\n")'

run_logged preflight-pandoc env \
  R_PROFILE_USER="$profile" TMPDIR="$tmp_root" \
  R_USER_CACHE_DIR="$cache_root" OCULAR_RUN_LIVE=false \
  Rscript -e 'if (!rmarkdown::pandoc_available()) stop("Pandoc is unavailable", call. = FALSE); cat("Pandoc ", as.character(rmarkdown::pandoc_version()), "\n", sep = "")'

ensure_directory "$pre_document_root/man" "$run_root"
cp "$package_root/DESCRIPTION" "$pre_document_root/DESCRIPTION"
cp "$package_root/NAMESPACE" "$pre_document_root/NAMESPACE"
cp "$package_root"/man/*.Rd "$pre_document_root/man/"

run_logged document env \
  R_PROFILE_USER="$profile" TMPDIR="$tmp_root" \
  R_USER_CACHE_DIR="$cache_root" OCULAR_RUN_LIVE=false \
  Rscript -e 'devtools::document("ocular")'

(
  CDPATH= cd -- "$pre_document_root"
  find man -type f -name '*.Rd' -print | sort
) >"$run_root/pre-documentation-files.txt"
(
  CDPATH= cd -- "$package_root"
  find man -type f -name '*.Rd' -print | sort
) >"$run_root/post-documentation-files.txt"
cmp -s "$run_root/pre-documentation-files.txt" \
  "$run_root/post-documentation-files.txt" ||
  fail "Documentation generation added or removed help topics."

: >"$run_root/generated-documentation-changes.txt"
for relative in DESCRIPTION NAMESPACE
do
  if ! cmp -s "$pre_document_root/$relative" "$package_root/$relative"; then
    printf '%s\n' "$relative" >>"$run_root/generated-documentation-changes.txt"
  fi
done
while IFS= read -r relative
do
  if ! cmp -s "$pre_document_root/$relative" "$package_root/$relative"; then
    printf '%s\n' "$relative" >>"$run_root/generated-documentation-changes.txt"
  fi
done <"$run_root/post-documentation-files.txt"

{
  diff -u "$pre_document_root/DESCRIPTION" "$package_root/DESCRIPTION" || true
  diff -u "$pre_document_root/NAMESPACE" "$package_root/NAMESPACE" || true
  diff -ru "$pre_document_root/man" "$package_root/man" || true
} >"$run_root/generated-documentation.diff"

printf '%s\n' 'man/rs_params.Rd' 'man/split_area.Rd' \
  >"$run_root/expected-documentation-changes.txt"
sort "$run_root/generated-documentation-changes.txt" \
  >"$run_root/generated-documentation-changes.sorted.txt"
sort "$run_root/expected-documentation-changes.txt" \
  >"$run_root/expected-documentation-changes.sorted.txt"
comm -23 "$run_root/generated-documentation-changes.sorted.txt" \
  "$run_root/expected-documentation-changes.sorted.txt" \
  >"$run_root/unexpected-documentation-changes.txt"
if [ -s "$run_root/unexpected-documentation-changes.txt" ]; then
    printf 'Observed generated changes:\n' >&2
    cat "$run_root/generated-documentation-changes.txt" >&2
    fail "Generated documentation changes differ from the reviewed expectation."
fi

run_logged focused-split-tests env \
  R_PROFILE_USER="$profile" TMPDIR="$tmp_root" \
  R_USER_CACHE_DIR="$cache_root" OCULAR_RUN_LIVE=false \
  Rscript -e 'devtools::test("ocular", filter = "split-guards", reporter = "summary", stop_on_failure = TRUE)'

run_logged full-tests env \
  R_PROFILE_USER="$profile" TMPDIR="$tmp_root" \
  R_USER_CACHE_DIR="$cache_root" OCULAR_RUN_LIVE=false \
  Rscript -e 'devtools::test("ocular", reporter = "summary", stop_on_failure = TRUE)'

ocular_version=$(env R_PROFILE_USER="$profile" TMPDIR="$tmp_root" \
  R_USER_CACHE_DIR="$cache_root" OCULAR_RUN_LIVE=false \
  Rscript -e 'cat(read.dcf("ocular/DESCRIPTION")[1L, "Version"])')
[ -n "$ocular_version" ] || fail "Could not read the ocular version."
[ "$ocular_version" = "$source_version" ] ||
  fail "The package version changed during verification."
tarball="$build_root/ocular_${ocular_version}.tar.gz"
[ ! -e "$tarball" ] && [ ! -L "$tarball" ] ||
  fail "Refusing to replace an existing tarball: $tarball"

run_logged_in build "$build_root" env \
  R_PROFILE_USER="$profile" TMPDIR="$tmp_root" \
  R_USER_CACHE_DIR="$cache_root" OCULAR_RUN_LIVE=false \
  R CMD build "$package_root"
[ -f "$tarball" ] && [ ! -L "$tarball" ] ||
  fail "R CMD build did not create the expected tarball: $tarball"

tar -tf "$tarball" >"$run_root/tarball-contents.txt"
for required in \
  ocular/tests/testthat/test-split-guards.R \
  ocular/man/split_area.Rd \
  ocular/man/rs_params.Rd \
  ocular/vignettes/ocular-introduction.Rmd \
  ocular/vignettes/ocular-build-pipeline.Rmd \
  ocular/vignettes/ocular-experimental-features.Rmd \
  ocular/inst/doc/ocular-introduction.html \
  ocular/inst/doc/ocular-build-pipeline.html \
  ocular/inst/doc/ocular-experimental-features.html
do
  grep -Fqx "$required" "$run_root/tarball-contents.txt" ||
    fail "The source tarball is missing $required"
done

if grep -E '/([.]DS_Store|[.]Rhistory|[.]Rproj[.]user|[.]check-artifacts|[.]test-artifacts)(/|$)|/tests/manual(/|$)|/vignettes/[^/]+[.]html$|/ocular[.]Rproj$|/LICENSE[.]md$' \
  "$run_root/tarball-contents.txt" >"$run_root/disallowed-tarball-entries.txt"
then
  cat "$run_root/disallowed-tarball-entries.txt" >&2
  fail "The source tarball contains excluded development artefacts."
fi

tarball_md5=$(env R_PROFILE_USER="$profile" TMPDIR="$tmp_root" \
  R_USER_CACHE_DIR="$cache_root" OCULAR_RUN_LIVE=false \
  OCULAR_CHECKED_TARBALL="$tarball" \
  Rscript -e 'value <- unname(tools::md5sum(Sys.getenv("OCULAR_CHECKED_TARBALL"))); if (length(value) != 1L || is.na(value) || !nzchar(value)) stop("Could not calculate tarball checksum", call. = FALSE); cat(value)')
[ -n "$tarball_md5" ] || fail "The tarball checksum is empty."
printf '%s  %s\n' "$tarball_md5" "$tarball" >"$run_root/tarball-md5.txt"

run_logged check env \
  R_PROFILE_USER="$profile" TMPDIR="$tmp_root" \
  R_USER_CACHE_DIR="$cache_root" OCULAR_RUN_LIVE=false \
  _R_CHECK_CRAN_INCOMING_REMOTE_=false \
  R CMD check --output="$check_root" --no-manual "$tarball"

check_log="$check_root/ocular.Rcheck/00check.log"
[ -f "$check_log" ] || fail "R CMD check did not write $check_log"
grep -Fq 'Status: OK' "$check_log" || {
  cat "$check_log" >&2
  fail "R CMD check did not finish with Status: OK."
}

preinstall_md5=$(env R_PROFILE_USER="$profile" TMPDIR="$install_tmp" \
  R_USER_CACHE_DIR="$install_cache" OCULAR_RUN_LIVE=false \
  OCULAR_CHECKED_TARBALL="$tarball" \
  Rscript -e 'cat(unname(tools::md5sum(Sys.getenv("OCULAR_CHECKED_TARBALL"))))')
[ "$preinstall_md5" = "$tarball_md5" ] ||
  fail "The checked tarball changed before installation."

run_logged install env \
  R_PROFILE_USER="$profile" TMPDIR="$install_tmp" \
  R_USER_CACHE_DIR="$install_cache" OCULAR_RUN_LIVE=false \
  R CMD INSTALL --library="$install_library" "$tarball"

install_receipt="$install_root/install-receipt.csv"
[ ! -e "$install_receipt" ] && [ ! -L "$install_receipt" ] ||
  fail "Refusing to replace install receipt: $install_receipt"
run_logged_in install-receipt "$package_root" env \
  R_PROFILE_USER="$profile" TMPDIR="$install_tmp" \
  R_USER_CACHE_DIR="$install_cache" OCULAR_RUN_LIVE=false \
  OCULAR_WRITE_INSTALL_RECEIPT=true \
  OCULAR_CHECKED_TARBALL="$tarball" \
  OCULAR_CHECKED_TARBALL_MD5="$tarball_md5" \
  OCULAR_VALIDATION_LIBRARY="$install_library" \
  OCULAR_INSTALL_RECEIPT="$install_receipt" \
  Rscript tests/manual/write-install-receipt.R

run_logged verify-install env \
  R_PROFILE_USER="$profile" TMPDIR="$install_tmp" \
  R_USER_CACHE_DIR="$install_cache" OCULAR_RUN_LIVE=false \
  OCULAR_VALIDATION_LIBRARY="$install_library" \
  OCULAR_EXPECTED_VERSION="$ocular_version" \
  R_LIBS="$install_library${R_LIBS:+:$R_LIBS}" \
  Rscript -e 'library_path <- normalizePath(Sys.getenv("OCULAR_VALIDATION_LIBRARY"), winslash = "/", mustWork = TRUE); installed_path <- normalizePath(find.package("ocular"), winslash = "/", mustWork = TRUE); expected_path <- normalizePath(file.path(library_path, "ocular"), winslash = "/", mustWork = TRUE); version <- as.character(packageVersion("ocular")); stopifnot(identical(installed_path, expected_path), identical(version, Sys.getenv("OCULAR_EXPECTED_VERSION"))); cat("Installed path: ", installed_path, "\nVersion: ", version, "\n", sep = "")'

comparison_root="$artifacts_root/split-design-comparison"
ensure_directory "$comparison_root" "$package_root"
find "$comparison_root" -mindepth 1 -maxdepth 1 -type d -name '20*' -print \
  | sort >"$run_root/diagnostic-runs-before.txt"

run_logged_in split-diagnostic "$package_root" env \
  R_PROFILE_USER="$profile" TMPDIR="$diagnostic_launch_tmp" \
  R_USER_CACHE_DIR="$diagnostic_launch_cache" OCULAR_RUN_LIVE=false \
  Rscript -e 'source("tests/manual/split-design-comparison.R")'

find "$comparison_root" -mindepth 1 -maxdepth 1 -type d -name '20*' -print \
  | sort >"$run_root/diagnostic-runs-after.txt"
comm -13 "$run_root/diagnostic-runs-before.txt" \
  "$run_root/diagnostic-runs-after.txt" >"$run_root/new-diagnostic-runs.txt"
new_diagnostic_count=$(wc -l <"$run_root/new-diagnostic-runs.txt" | tr -d ' ')
[ "$new_diagnostic_count" = "1" ] || {
  cat "$run_root/new-diagnostic-runs.txt" >&2
  fail "Expected exactly one new split diagnostic run; found $new_diagnostic_count."
}
diagnostic_root=$(sed -n '1p' "$run_root/new-diagnostic-runs.txt")
case "$diagnostic_root" in
  "$comparison_root"/*) ;;
  *) fail "The new diagnostic resolves outside its package-local root." ;;
esac
[ ! -L "$diagnostic_root" ] || fail "The new diagnostic is a symlink."

for required in \
  execution-identity.csv historical-provenance.csv \
  reconstructed-fixtures.csv effective-parameters.csv \
  suffix-stage-schedule.csv source-tree-context.csv split-results.csv \
  split-masks.rds manifest.csv interpretation.txt sessionInfo.txt
do
  [ -f "$diagnostic_root/$required" ] ||
    fail "The diagnostic is missing $required"
done
[ -d "$diagnostic_root/figures" ] || fail "The diagnostic has no figures directory."
figure_count=$(find "$diagnostic_root/figures" -type f -name '*.png' | wc -l | tr -d ' ')
[ "$figure_count" -gt 0 ] || fail "The diagnostic produced no PNG figures."
printf '%s\n' "$diagnostic_root" >"$run_root/diagnostic-root.txt"
printf '%s\n' "$figure_count" >"$run_root/diagnostic-figure-count.txt"

run_logged diagnostic-summary env \
  R_PROFILE_USER="$profile" TMPDIR="$tmp_root" \
  R_USER_CACHE_DIR="$cache_root" OCULAR_RUN_LIVE=false \
  OCULAR_DIAGNOSTIC_ROOT="$diagnostic_root" \
  Rscript -e 'path <- file.path(Sys.getenv("OCULAR_DIAGNOSTIC_ROOT"), "split-results.csv"); x <- utils::read.csv(path, stringsAsFactors = FALSE); cat("Rows: ", nrow(x), "\n", sep = ""); for (name in intersect(c("outside_input_alive", "seed_lost", "all_removed", "unchanged", "exact_seed_region"), names(x))) { values <- x[[name]]; if (is.logical(values)) { cat(name, ": TRUE=", sum(values %in% TRUE, na.rm = TRUE), ", FALSE=", sum(values %in% FALSE, na.rm = TRUE), ", NA=", sum(is.na(values)), "\n", sep = "") } else { cat(name, ": min=", min(values, na.rm = TRUE), ", max=", max(values, na.rm = TRUE), "\n", sep = "") } }'

{
  printf 'Offline ocular release verification completed successfully.\n'
  printf 'Run root: %s\n' "$run_root"
  printf 'Generated documentation changes: %s\n' \
    "$(tr '\n' ' ' <"$run_root/generated-documentation-changes.txt")"
  printf 'Focused test log: %s\n' "$logs_root/focused-split-tests.log"
  printf 'Full test log: %s\n' "$logs_root/full-tests.log"
  printf 'Checked tarball: %s\n' "$tarball"
  printf 'Tarball MD5: %s\n' "$tarball_md5"
  printf 'R CMD check log: %s\n' "$check_log"
  printf 'Install receipt: %s\n' "$install_receipt"
  printf 'Diagnostic root: %s\n' "$diagnostic_root"
  printf 'Diagnostic PNG count: %s\n' "$figure_count"
} | tee "$run_root/verification-summary.txt"
