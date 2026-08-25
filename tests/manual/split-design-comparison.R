#!/usr/bin/env Rscript

# Manual, offline comparison of the current split behaviour.
#
# The split fixtures in this file are reconstructions from the July prose
# record. The historical test source, grid dimensions, exact cell coordinates,
# and feature values are not available. Results are therefore evidence for
# author review, not release assertions. This script installs nothing, makes
# no network request, and writes generated artefacts only below
# ocular/.check-artifacts/.

local({

.enabled <- function(x) {
  identical(tolower(trimws(x)), "true")
}

.assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

.inside <- function(path, root, must_work = FALSE) {
  path <- normalizePath(path, winslash = "/", mustWork = must_work)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  identical(path, root) || startsWith(path, paste0(root, "/"))
}

.is_symlink <- function(path) {
  target <- Sys.readlink(path)
  length(target) == 1L && !is.na(target) && nzchar(target)
}

.is_scalar_text <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.ensure_directory <- function(path, root) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  candidate <- normalizePath(path, winslash = "/", mustWork = FALSE)
  .assert(.inside(candidate, root),
          "Refusing to create a diagnostic directory outside the package.")

  relative <- substring(candidate, nchar(root) + 2L)
  parts <- strsplit(relative, "/", fixed = TRUE)[[1L]]
  parts <- parts[nzchar(parts)]
  current <- root

  for (part in parts) {
    next_path <- file.path(current, part)
    .assert(!.is_symlink(next_path),
            paste0("Refusing symlinked diagnostic output: ", next_path))
    if (!dir.exists(next_path)) {
      .assert(!file.exists(next_path),
              paste0("A non-directory blocks diagnostic output: ", next_path))
      dir.create(next_path, showWarnings = FALSE)
    }
    .assert(dir.exists(next_path),
            paste0("Could not create diagnostic directory: ", next_path))
    resolved <- normalizePath(next_path, winslash = "/", mustWork = TRUE)
    .assert(.inside(resolved, root, must_work = TRUE),
            "A diagnostic output directory resolves outside the package.")
    current <- resolved
  }

  current
}

.launch_worker <- function(package_root, script_path) {
  artifacts_root <- .ensure_directory(
    file.path(package_root, ".check-artifacts"), package_root)
  comparison_root <- .ensure_directory(
    file.path(artifacts_root, "split-design-comparison"), package_root)
  session_temp <- .ensure_directory(
    file.path(comparison_root, "r-session-temp"), package_root)
  session_cache <- .ensure_directory(
    file.path(comparison_root, "r-session-cache"), package_root)

  rscript <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )
  .assert(isTRUE(file_test("-f", rscript)),
          paste0("Could not find Rscript at ", rscript, "."))

  status <- system2(
    rscript,
    args = c("--vanilla", shQuote(script_path)),
    stdout = "",
    stderr = "",
    env = c(
      paste0("TMPDIR=", shQuote(session_temp)),
      paste0("R_USER_CACHE_DIR=", shQuote(session_cache)),
      "OCULAR_SPLIT_DIAGNOSTIC_WORKER=true"
    ),
    wait = TRUE
  )
  .assert(
    is.numeric(status) && length(status) == 1L && !is.na(status) &&
      status == 0L,
    paste0(
      "The workspace-local split diagnostic process failed with status ",
      status, ". Review the preceding error."
    )
  )
  invisible(status)
}

.load_current_source <- function(package_root, source_version) {
  .assert(.is_scalar_text(source_version),
          "DESCRIPTION does not contain one valid package version.")
  .assert(
    requireNamespace("pkgload", quietly = TRUE),
    paste0(
      "The split diagnostic requires the already-installed pkgload package. ",
      "Nothing was installed or downloaded."
    )
  )

  tryCatch(
    pkgload::load_all(
      package_root,
      reset = TRUE,
      export_all = FALSE,
      helpers = FALSE,
      attach_testthat = FALSE,
      quiet = TRUE
    ),
    error = function(e) {
      stop(
        "Could not load the current ocular source: ",
        conditionMessage(e),
        "\nNothing was installed or downloaded.",
        call. = FALSE
      )
    }
  )

  namespace <- asNamespace("ocular")
  loaded_path <- unname(normalizePath(
    getNamespaceInfo(namespace, "path"), winslash = "/", mustWork = TRUE))
  loaded_name <- unname(as.character(getNamespaceName(namespace)))
  loaded_version <- unname(as.character(getNamespaceVersion(namespace)))
  .assert(
    .is_scalar_text(loaded_name) && .is_scalar_text(loaded_version),
    "The loaded namespace specification has no scalar package identity."
  )
  .assert(
    identical(loaded_path, package_root),
    paste0(
      "pkgload did not load ocular from the current package root. Loaded: ",
      loaded_path
    )
  )
  .assert(
    identical(loaded_name, "ocular"),
    paste0("pkgload loaded namespace ", loaded_name, " instead of ocular.")
  )
  .assert(
    identical(loaded_version, source_version),
    paste0(
      "Loaded ocular version ", loaded_version,
      " does not match source version ", source_version, "."
    )
  )

  list(
    namespace = namespace,
    path = loaded_path,
    version = loaded_version,
    method = "pkgload::load_all()"
  )
}

.slug <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "-", x))
  gsub("(^-+|-+$)", "", x)
}

.write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

.value_text <- function(x) {
  if (is.null(x)) return("NULL")
  paste(deparse(x, width.cutoff = 500L), collapse = " ")
}

.params_table <- function(params, scenario, fixture, configuration) {
  .assert(is.list(params) && !is.null(names(params)),
          paste0(configuration, ": effective parameters are not a named list."))
  user_set <- attr(params, "user_set", exact = TRUE)
  if (is.null(user_set)) user_set <- character()
  .assert(is.character(user_set) && !anyNA(user_set) &&
            all(nzchar(user_set)) && anyDuplicated(user_set) == 0L &&
            all(user_set %in% names(params)),
          paste0(configuration, ": parameter ownership is not character data."))
  data.frame(
    scenario = scenario,
    fixture = fixture,
    configuration = configuration,
    parameter = names(params),
    value = vapply(params, .value_text, character(1L)),
    explicitly_set = names(params) %in% user_set,
    stringsAsFactors = FALSE
  )
}

.save_png <- function(path, draw, width = 1600L, height = 1200L) {
  grDevices::png(path, width = width, height = height, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
  invisible(path)
}

.draw_split_mask <- function(mask, fixture, title, show_legend = TRUE) {
  code <- matrix(0L, nrow(mask), ncol(mask))
  code[fixture$seed_region & !mask] <- 1L
  code[fixture$seed_region & mask] <- 2L
  code[fixture$divergent_region & !mask] <- 3L
  code[fixture$divergent_region & mask] <- 4L
  code[!fixture$alive & mask] <- 5L

  colours <- c("#F7F7F7", "#E69F00", "#0072B2",
               "#D9D9D9", "#D55E00", "#CC79A7")
  display <- t(code[nrow(code):1L, , drop = FALSE])
  graphics::image(
    x = seq_len(ncol(code)), y = seq_len(nrow(code)), z = display,
    breaks = seq(-0.5, 5.5, by = 1), col = colours, asp = 1,
    axes = FALSE, xlab = "", ylab = "", main = title
  )
  graphics::box()
  graphics::points(
    fixture$seed[2L], nrow(code) - fixture$seed[1L] + 1L,
    pch = 4, lwd = 2, cex = 1.1
  )
  if (isTRUE(show_legend)) {
    graphics::legend(
      "bottom", inset = -0.28, xpd = NA, horiz = TRUE, bty = "n", cex = 0.72,
      fill = colours[c(2L, 3L, 4L, 5L, 6L)],
      legend = c("seed region removed", "seed region retained",
                 "divergent region removed", "divergent region retained",
                 "outside input alive")
    )
  }
  invisible()
}

.save_split_mask <- function(path, mask, fixture, title) {
  .save_png(path, function() {
    old <- graphics::par(mar = c(3, 3, 4, 1))
    on.exit(graphics::par(old), add = TRUE)
    .draw_split_mask(mask, fixture, title, show_legend = TRUE)
  })
}

.save_split_sequence <- function(path, masks, fixture, title) {
  panel_cols <- 3L
  panel_rows <- as.integer(ceiling(length(masks) / panel_cols))
  .save_png(path, function() {
    old <- graphics::par(
      mfrow = c(panel_rows, panel_cols), mar = c(1.5, 1.5, 3, 0.5),
      oma = c(0, 0, 3, 0)
    )
    on.exit(graphics::par(old), add = TRUE)
    for (name in names(masks))
      .draw_split_mask(masks[[name]], fixture, name, show_legend = FALSE)
    unused <- panel_rows * panel_cols - length(masks)
    if (unused > 0L) for (i in seq_len(unused)) graphics::plot.new()
    graphics::mtext(title, outer = TRUE, cex = 1.1)
  }, width = 1800L, height = max(1000L, panel_rows * 520L))
}

.rectangle <- function(n, rows, columns) {
  out <- matrix(FALSE, n, n)
  out[rows, columns] <- TRUE
  out
}

.make_fixture <- function(kind, n = 40L) {
  feature <- matrix(0.1, n, n)

  if (identical(kind, "side_by_side")) {
    seed_region <- .rectangle(n, 13:27, 10:19)
    divergent_region <- .rectangle(n, 13:27, 20:30)
    seed <- c(20L, 15L)
    expected <- c(seed = 150L, divergent = 165L, total = 315L)
  } else if (identical(kind, "stacked")) {
    seed_region <- .rectangle(n, 10:19, 13:27)
    divergent_region <- .rectangle(n, 20:30, 13:27)
    seed <- c(15L, 20L)
    expected <- c(seed = 150L, divergent = 165L, total = 315L)
  } else if (identical(kind, "minority_3x3")) {
    alive <- .rectangle(n, 13:27, 10:30)
    divergent_region <- .rectangle(n, 19:21, 19:21)
    seed_region <- alive & !divergent_region
    seed <- c(20L, 15L)
    expected <- c(seed = 306L, divergent = 9L, total = 315L)
  } else {
    stop("Unknown reconstructed fixture: ", kind, call. = FALSE)
  }

  alive <- seed_region | divergent_region
  feature[seed_region] <- 0.7
  feature[divergent_region] <- 0.2

  list(
    name = kind,
    n = n,
    seed = seed,
    seed_region = seed_region,
    divergent_region = divergent_region,
    alive = alive,
    feature = feature,
    seed_signature = 0.7,
    expected = expected
  )
}

.validate_fixture <- function(fixture) {
  n <- fixture$n
  .assert(identical(dim(fixture$alive), c(n, n)),
          paste0(fixture$name, ": invalid alive-mask dimensions."))
  .assert(!any(fixture$seed_region & fixture$divergent_region),
          paste0(fixture$name, ": reconstructed regions overlap."))
  .assert(identical(fixture$alive,
                    fixture$seed_region | fixture$divergent_region),
          paste0(fixture$name, ": alive mask is not the region union."))
  .assert(sum(fixture$seed_region) == fixture$expected[["seed"]],
          paste0(fixture$name, ": seed-region count is wrong."))
  .assert(sum(fixture$divergent_region) == fixture$expected[["divergent"]],
          paste0(fixture$name, ": divergent-region count is wrong."))
  .assert(sum(fixture$alive) == fixture$expected[["total"]],
          paste0(fixture$name, ": total reconstructed count is wrong."))
  .assert(fixture$seed_region[fixture$seed[1L], fixture$seed[2L]],
          paste0(fixture$name, ": seed is not in the seed region."))
  .assert(all(is.finite(fixture$feature[fixture$alive])),
          paste0(fixture$name, ": alive feature values are not finite."))
  invisible(fixture)
}

.fixture_cells <- function(fixture) {
  cells <- which(matrix(TRUE, fixture$n, fixture$n), arr.ind = TRUE)
  region <- rep("outside", nrow(cells))
  region[fixture$seed_region[cells]] <- "seed_region"
  region[fixture$divergent_region[cells]] <- "divergent_region"
  data.frame(
    fixture = fixture$name,
    row = cells[, "row"],
    column = cells[, "col"],
    region = region,
    alive = fixture$alive[cells],
    feature_value = fixture$feature[cells],
    seed = cells[, "row"] == fixture$seed[1L] &
      cells[, "col"] == fixture$seed[2L],
    stringsAsFactors = FALSE
  )
}

.validate_mask <- function(mask, fixture, label) {
  .assert(is.matrix(mask) && is.logical(mask),
          paste0(label, ": output is not a logical matrix."))
  .assert(identical(dim(mask), c(fixture$n, fixture$n)),
          paste0(label, ": output dimensions changed."))
  .assert(!anyNA(mask), paste0(label, ": output contains NA."))
  mask
}

.make_split_rs <- function(fixture, gate, multiple, linear,
                           new_ocular, validate_ocular) {
  rs <- new_ocular()
  rs$spec$index_name <- "EVI2"
  rs$geom$source <- "synthetic"
  rs$geom$pixel_size_m <- 10
  rs$geom$nr <- fixture$n
  rs$geom$nc <- fixture$n
  rs$geom$centre_rc <- fixture$seed
  rs$internals$feat_array <- array(
    fixture$feature, dim = c(fixture$n, fixture$n, 1L)
  )
  rs$internals$mu <- fixture$seed_signature
  rs$state$alive_mat <- fixture$alive
  rs$state$decided_mat <- matrix(FALSE, fixture$n, fixture$n)
  rs$params <- ocular::rs_params(
    split_sensitivity = 0.05,
    split_gate = as.integer(gate),
    split_linear = linear,
    multiple_areas = multiple
  )
  validate_ocular(rs)
  rs
}

.july_reference_total <- function(fixture, scenario, gate, multiple, linear) {
  if (!identical(scenario, "direct_split") || gate != 2L || isTRUE(linear))
    return(NA_integer_)
  if (identical(fixture$name, "side_by_side") && !isTRUE(multiple))
    return(150L)
  if (identical(fixture$name, "side_by_side") && isTRUE(multiple))
    return(165L)
  if (identical(fixture$name, "stacked") && !isTRUE(multiple))
    return(150L)
  if (identical(fixture$name, "minority_3x3") && !isTRUE(multiple))
    return(315L)
  NA_integer_
}

.measure_split <- function(mask, fixture, scenario, gate, multiple, linear,
                           pass, stage) {
  seed_alive <- mask[fixture$seed[1L], fixture$seed[2L]]
  seed_total <- sum(fixture$seed_region)
  divergent_total <- sum(fixture$divergent_region)
  seed_alive_n <- sum(mask & fixture$seed_region)
  divergent_alive_n <- sum(mask & fixture$divergent_region)

  data.frame(
    fixture = fixture$name,
    scenario = scenario,
    gate = as.integer(gate),
    multiple_areas = multiple,
    split_linear = linear,
    pass = as.integer(pass),
    stage = stage,
    input_alive = sum(fixture$alive),
    output_alive = sum(mask),
    seed_region_alive = seed_alive_n,
    divergent_region_alive = divergent_alive_n,
    outside_input_alive = sum(mask & !fixture$alive),
    seed_region_retained_fraction = seed_alive_n / seed_total,
    divergent_region_removed_fraction =
      1 - divergent_alive_n / divergent_total,
    seed_alive = seed_alive,
    seed_lost = !seed_alive,
    unchanged = identical(mask, fixture$alive),
    exact_seed_region = identical(mask, fixture$seed_region),
    all_removed = !any(mask),
    july_prose_reference_total = .july_reference_total(
      fixture, scenario, gate, multiple, linear
    ),
    stringsAsFactors = FALSE
  )
}

.run_boundary_suffix <- function(rs, fixture) {
  states <- list(input_post_segment = rs$state$alive_mat)

  out <- ocular::trace_perimeter(rs)
  states$after_post_segment_trace <- .validate_mask(
    out$state$alive_mat, fixture, "post-segment trace"
  )

  out <- ocular::split_area(out)
  states$after_split_1 <- .validate_mask(
    out$state$alive_mat, fixture, "boundary suffix split 1"
  )
  out <- ocular::trace_perimeter(out)
  out <- ocular::segment_interior(out)
  out <- ocular::trace_perimeter(out)

  out <- ocular::split_area(out)
  states$after_split_2 <- .validate_mask(
    out$state$alive_mat, fixture, "boundary suffix split 2"
  )
  out <- ocular::trace_perimeter(out)

  out <- ocular::split_area(out)
  states$after_split_3 <- .validate_mask(
    out$state$alive_mat, fixture, "boundary suffix split 3"
  )
  out <- ocular::trace_perimeter(out)

  out <- ocular::split_area(out)
  states$after_split_4 <- .validate_mask(
    out$state$alive_mat, fixture, "boundary suffix split 4"
  )
  out <- ocular::trace_perimeter(
    out, params = list(perimeter_margins = c(0.3, 0.9))
  )

  out <- ocular::split_area(out)
  states$after_split_5 <- .validate_mask(
    out$state$alive_mat, fixture, "boundary suffix split 5"
  )
  out <- ocular::trace_perimeter(
    out, params = list(perimeter_margins = c(0.5, 0.5, 1))
  )
  states$final_after_closing_trace <- .validate_mask(
    out$state$alive_mat, fixture, "boundary suffix final trace"
  )

  list(object = out, states = states)
}

.source_tree_context <- function(package_root) {
  r_files <- sort(list.files(
    file.path(package_root, "R"), pattern = "[.]R$", full.names = FALSE))
  relative <- c(
    "DESCRIPTION", "NAMESPACE", file.path("R", r_files),
    "tests/manual/split-design-comparison.R"
  )
  paths <- file.path(package_root, relative)
  present <- vapply(paths, function(path) isTRUE(file_test("-f", path)),
                    logical(1L))
  if (any(present)) {
    resolved <- normalizePath(paths[present], winslash = "/", mustWork = TRUE)
    .assert(all(vapply(
      resolved, .inside, logical(1L), root = package_root, must_work = TRUE
    )), "A source-tree context file resolves outside the package.")
  }
  hashes <- rep(NA_character_, length(paths))
  if (any(present))
    hashes[present] <- unname(tools::md5sum(paths[present]))
  data.frame(path = relative, present = present, md5 = hashes,
             stringsAsFactors = FALSE)
}

.historical_provenance <- function(workspace_root) {
  design_notes <- file.path(
    "vignette and manuscript", "oculaR", "ocular-docs-review-August-2026",
    "DESIGN-NOTES.md"
  )
  july_manuscript <- file.path(
    "vignette and manuscript", "oculaR", "ocular-docs-review-August-2026",
    "all-manuscript", "manuscript-july", "manuscript_walker-rebuild.md"
  )
  july_package <- file.path(
    "vignette and manuscript", "oculaR", "ocular-docs-review-August-2026",
    "all-manuscript", "manuscript-july", "ocular-package.Rmd"
  )
  out <- data.frame(
    path = c(
      design_notes,
      july_manuscript, july_manuscript,
      july_package
    ),
    lines = c("138-164", "117-128", "142-145", "178-192"),
    recorded_evidence = c(
      paste0("near/far-half split rule, per-axis acceptance, destructive ",
             "315-cell ablation, retained totals, and gate-2 decision"),
      paste0("guarded split semantics, gate-2 default, cut modes, and ",
             "per-axis acceptance"),
      paste0("40 x 40 is stated for DISC-EXACT; SPLIT-2LOBE counts, ",
             "orientations, minority size, retained totals, and ablation ",
             "are stated separately without a split-grid dimension"),
      paste0("SPLIT-2LOBE counts, orientations, minority size, retained ",
             "totals, and destructive ablation")
    ),
    use_in_diagnostic = c(
      "historical split comparison target",
      "historical split comparison target",
      "fixture facts and explicit 40 x 40 reconstruction limit",
      "corroborating historical fixture account"
    ),
    stringsAsFactors = FALSE
  )
  full_paths <- file.path(workspace_root, out$path)
  out$present <- vapply(
    full_paths, function(path) isTRUE(file_test("-f", path)), logical(1L))
  if (any(out$present)) {
    resolved <- normalizePath(
      full_paths[out$present], winslash = "/", mustWork = TRUE)
    .assert(all(vapply(
      resolved, .inside, logical(1L), root = workspace_root, must_work = TRUE
    )), "A historical provenance source resolves outside the workspace.")
  }
  out$md5 <- NA_character_
  if (any(out$present))
    out$md5[out$present] <- unname(tools::md5sum(full_paths[out$present]))
  .assert(all(out$present),
          paste0("One or more historical provenance sources are unavailable; ",
                 "run this diagnostic from the complete ocular workspace."))
  .assert(!anyNA(out$md5) && all(nzchar(out$md5)),
          "A historical provenance checksum could not be calculated.")
  out
}

.suffix_stage_schedule <- function() {
  data.frame(
    order = seq_len(14L),
    call = c(
      "manually reconstructed post-segment input",
      "trace_perimeter()", "split_area()", "trace_perimeter()",
      "segment_interior()", "trace_perimeter()", "split_area()",
      "trace_perimeter()", "split_area()", "trace_perimeter()",
      "split_area()", "trace_perimeter()", "split_area()",
      "trace_perimeter()"
    ),
    parameter_override = c(
      NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
      "perimeter_margins=c(0.3, 0.9)", NA,
      "perimeter_margins=c(0.5, 0.5, 1)"
    ),
    saved_state = c(
      "input_post_segment", "after_post_segment_trace", "after_split_1",
      NA, NA, NA, "after_split_2", NA, "after_split_3", NA,
      "after_split_4", NA, "after_split_5",
      "final_after_closing_trace"
    ),
    context = paste0(
      "selected cumulative states from the current default stage order ",
      "(cleanup_boundary=2L), beginning with a reconstructed mask"
    ),
    stringsAsFactors = FALSE
  )
}

.run <- function() {
  package_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  description_path <- file.path(package_root, "DESCRIPTION")
  .assert(file.exists(description_path),
          "Run this script from the ocular package root.")
  description <- read.dcf(description_path)
  .assert(identical(unname(description[1L, "Package"]), "ocular"),
          "The working directory is not the ocular package root.")
  workspace_root <- dirname(package_root)
  script_path <- normalizePath(
    file.path(package_root, "tests", "manual", "split-design-comparison.R"),
    winslash = "/", mustWork = TRUE)
  .assert(.inside(script_path, package_root, must_work = TRUE) &&
            !.is_symlink(script_path),
          "The split diagnostic script must be a regular package-local file.")

  is_worker <- .enabled(
    Sys.getenv("OCULAR_SPLIT_DIAGNOSTIC_WORKER", unset = ""))
  if (!is_worker)
    return(.launch_worker(package_root, script_path))

  options(repos = character())
  .assert(
    .inside(tempdir(), package_root, must_work = TRUE),
    paste0(
      "The workspace-local diagnostic process did not receive a package-local ",
      "temporary directory."
    )
  )
  session_cache_raw <- Sys.getenv("R_USER_CACHE_DIR", unset = "")
  .assert(nzchar(session_cache_raw),
          "The workspace-local diagnostic process has no cache directory.")
  session_cache <- normalizePath(
    session_cache_raw, winslash = "/", mustWork = TRUE)
  .assert(.inside(session_cache, package_root, must_work = TRUE),
          "The diagnostic cache directory is outside the ocular package.")

  run_id <- paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-", Sys.getpid())
  artifacts_root <- .ensure_directory(
    file.path(package_root, ".check-artifacts"), package_root)
  comparison_root <- .ensure_directory(
    file.path(artifacts_root, "split-design-comparison"), package_root)
  output_candidate <- file.path(comparison_root, run_id)
  .assert(!file.exists(output_candidate) && !dir.exists(output_candidate) &&
            !.is_symlink(output_candidate),
          paste0("Refusing to reuse an existing diagnostic run: ",
                 output_candidate))
  output_root <- .ensure_directory(output_candidate, package_root)
  cache_root <- .ensure_directory(file.path(output_root, "cache"), package_root)
  figures_root <- .ensure_directory(
    file.path(output_root, "figures"), package_root
  )
  terra_temp_root <- .ensure_directory(
    file.path(output_root, "terra-temp"), package_root)
  Sys.setenv(R_USER_CACHE_DIR = cache_root)

  source_version <- unname(description[1L, "Version"])
  load_info <- .load_current_source(package_root, source_version)
  terra::terraOptions(tempdir = terra_temp_root)

  .write_csv(data.frame(
    item = c(
      "execution_mode", "source_package_version", "loaded_package_version",
      "package_root", "loaded_namespace_path", "load_method",
      "worker_tempdir", "worker_startup_cache", "run_cache", "terra_tempdir",
      "active_library_paths"
    ),
    value = c(
      "current_source_tree", source_version, load_info$version,
      package_root, load_info$path, load_info$method,
      tempdir(), session_cache, cache_root, terra_temp_root,
      paste(.libPaths(), collapse = "; ")
    ),
    stringsAsFactors = FALSE
  ), file.path(output_root, "execution-identity.csv"))
  .write_csv(.historical_provenance(workspace_root),
             file.path(output_root, "historical-provenance.csv"))
  .write_csv(.suffix_stage_schedule(),
             file.path(output_root, "suffix-stage-schedule.csv"))

  namespace <- load_info$namespace
  new_ocular <- get(".newOcular", envir = namespace, inherits = FALSE)
  validate_ocular <- get(".validateOcular", envir = namespace,
                         inherits = FALSE)

  fixtures <- list(
    side_by_side = .make_fixture("side_by_side"),
    stacked = .make_fixture("stacked"),
    minority_3x3 = .make_fixture("minority_3x3")
  )
  invisible(lapply(fixtures, .validate_fixture))
  fixture_table <- do.call(rbind, lapply(fixtures, .fixture_cells))
  .write_csv(fixture_table, file.path(output_root, "reconstructed-fixtures.csv"))

  split_rows <- list()
  split_masks <- list()
  parameter_rows <- list()

  for (fixture in fixtures) {
    for (gate in c(1L, 2L)) {
      for (multiple in c(FALSE, TRUE)) {
        for (linear in c(FALSE, TRUE)) {
          label <- paste(
            "direct", fixture$name, paste0("gate", gate),
            paste0("multiple", tolower(multiple)),
            paste0("linear", tolower(linear)), sep = "-"
          )
          rs <- .make_split_rs(
            fixture, gate, multiple, linear, new_ocular, validate_ocular
          )
          parameter_rows[[length(parameter_rows) + 1L]] <- .params_table(
            rs$params, "direct_split", fixture$name, label
          )
          out <- ocular::split_area(rs)
          mask <- .validate_mask(out$state$alive_mat, fixture, label)
          row <- .measure_split(
            mask, fixture, "direct_split", gate, multiple, linear,
            pass = 1L, stage = "after_split"
          )
          .assert(row$outside_input_alive == 0L,
                  paste0(label, ": split_area revived cells outside its input."))
          split_rows[[length(split_rows) + 1L]] <- row
          split_masks[[label]] <- mask
          .save_split_mask(
            file.path(figures_root, paste0(.slug(label), ".png")),
            mask, fixture, label
          )
        }
      }
    }
  }

  for (fixture in fixtures[c("side_by_side", "stacked")]) {
    for (gate in c(1L, 2L)) {
      for (multiple in c(FALSE, TRUE)) {
        for (linear in c(FALSE, TRUE)) {
          label <- paste(
            "five-split-schedule", fixture$name, paste0("gate", gate),
            paste0("multiple", tolower(multiple)),
            paste0("linear", tolower(linear)), sep = "-"
          )
          out <- .make_split_rs(
            fixture, gate, multiple, linear, new_ocular, validate_ocular
          )
          parameter_rows[[length(parameter_rows) + 1L]] <- .params_table(
            out$params, "five_split_schedule", fixture$name, label
          )
          states <- list(input = out$state$alive_mat)
          for (pass in seq_len(5L)) {
            out <- ocular::split_area(out)
            mask <- .validate_mask(
              out$state$alive_mat, fixture,
              paste0(label, " pass ", pass)
            )
            .assert(!any(mask & !fixture$alive),
                    paste0(label, ": split pass revived outside-input cells."))
            stage <- paste0("after_split_", pass)
            states[[stage]] <- mask
            split_rows[[length(split_rows) + 1L]] <- .measure_split(
              mask, fixture, "five_split_schedule", gate, multiple, linear,
              pass = pass, stage = stage
            )
            split_masks[[paste(label, stage, sep = "-")]] <- mask
          }
          .save_split_sequence(
            file.path(figures_root, paste0(.slug(label), ".png")),
            states, fixture, label
          )
        }
      }
    }
  }

  for (fixture in fixtures[c("side_by_side", "stacked")]) {
    for (gate in c(1L, 2L)) {
      for (multiple in c(FALSE, TRUE)) {
        linear <- FALSE
        label <- paste(
          "reconstructed-post-segment-boundary-suffix", fixture$name,
          paste0("gate", gate), paste0("multiple", tolower(multiple)),
          sep = "-"
        )
        rs <- .make_split_rs(
          fixture, gate, multiple, linear, new_ocular, validate_ocular
        )
        parameter_rows[[length(parameter_rows) + 1L]] <- .params_table(
          rs$params, "reconstructed_post_segment_boundary_suffix",
          fixture$name, label
        )
        suffix <- .run_boundary_suffix(rs, fixture)
        for (stage in names(suffix$states)) {
          mask <- suffix$states[[stage]]
          split_no <- suppressWarnings(as.integer(sub("^after_split_", "", stage)))
          if (!startsWith(stage, "after_split_")) split_no <- NA_integer_
          split_rows[[length(split_rows) + 1L]] <- .measure_split(
            mask, fixture, "reconstructed_post_segment_boundary_suffix",
            gate, multiple, linear, pass = split_no, stage = stage
          )
          split_masks[[paste(label, stage, sep = "-")]] <- mask
        }
        .save_split_sequence(
          file.path(figures_root, paste0(.slug(label), ".png")),
          suffix$states, fixture, label
        )
      }
    }
  }

  split_table <- do.call(rbind, split_rows)
  .write_csv(split_table, file.path(output_root, "split-results.csv"))
  saveRDS(split_masks, file.path(output_root, "split-masks.rds"))

  manifest <- data.frame(
    item = c(
      "source_package_version", "loaded_package_version",
      "split_fixture_status", "split_fixture_grid",
      "split_feature_values", "split_sensitivity"
    ),
    value = c(
      source_version,
      load_info$version,
      paste0("lobe counts, orientations, and minority size from July prose; ",
             "grid and remaining fixture details reconstructed"),
      "40 x 40",
      "background=0.1; seed region=0.7; divergent region=0.2",
      "0.05"
    ),
    stringsAsFactors = FALSE
  )
  .write_csv(manifest, file.path(output_root, "manifest.csv"))
  .write_csv(do.call(rbind, parameter_rows),
             file.path(output_root, "effective-parameters.csv"))
  .write_csv(.source_tree_context(package_root),
             file.path(output_root, "source-tree-context.csv"))
  writeLines(capture.output(utils::sessionInfo()),
             file.path(output_root, "sessionInfo.txt"))
  writeLines(c(
    "Interpretation limits",
    "",
    paste0(
      "The 40 x 40 grid is an explicit reconstruction assumption. The lobe ",
      "counts, orientations, and minority size come from the historical ",
      "prose, but the grid dimensions, exact coordinates, feature values, ",
      "and some settings were not recovered and require author review."
    ),
    paste0(
      "The direct_split rows isolate one public split_area() call. The ",
      "five_split_schedule rows apply split_area() five times as an ",
      "accumulation stress test; they are not boundary_delineation()."
    ),
    paste0(
      "The reconstructed_post_segment_boundary_suffix rows mirror the five ",
      "split calls and interleaved cleanup stages after segment_area() in the ",
      "current default stage order (cleanup_boundary = 2L), under the ",
      "indicated split settings. They save selected cumulative states rather ",
      "than every interleaved stage, and do not establish that segment_area() ",
      "would produce the manually joined input mask."
    ),
    paste0(
      "Seed loss, complete removal, and disagreement with the July retained-",
      "cell totals are diagnostic outcomes and were deliberately not treated ",
      "as structural failures."
    ),
    paste0(
      "Both gate settings use the source-tree overall-median rule, protected ",
      "point-connected component, and provisional per-axis rollback. They do ",
      "not reproduce the July near/far-half and removed-cell-majority rules. ",
      "The preserved pre-correction run is 20260825-211428-33249."
    )
  ), file.path(output_root, "interpretation.txt"))

  cat("Diagnostic artefacts written to:\n", output_root, "\n", sep = "")
  invisible(output_root)
}

.run()
})
