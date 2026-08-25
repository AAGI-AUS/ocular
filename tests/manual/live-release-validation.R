#!/usr/bin/env Rscript

# Manual, opt-in release validation for live data services.
#
# This script is not part of the automated test suite. It performs network
# queries only when OCULAR_RUN_LIVE=true and never installs packages or DuckDB
# extensions. Run it from the ocular package root; see tests/manual/README.md.

.is_true <- function(x) {
  tolower(trimws(x)) %in% c("1", "true", "yes")
}

.live_enabled <- function(x) {
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

.ensure_directory <- function(path, root) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  candidate <- normalizePath(path, winslash = "/", mustWork = FALSE)
  .assert(.inside(candidate, root),
          "Refusing to create a validation directory outside the package.")
  relative <- substring(candidate, nchar(root) + 2L)
  parts <- strsplit(relative, "/", fixed = TRUE)[[1L]]
  parts <- parts[nzchar(parts)]
  current <- root
  for (part in parts) {
    next_path <- file.path(current, part)
    .assert(!nzchar(Sys.readlink(next_path)),
            paste0("Refusing symlinked validation output: ", next_path))
    if (!dir.exists(next_path)) {
      .assert(!file.exists(next_path),
              paste0("A non-directory blocks validation output: ", next_path))
      dir.create(next_path, showWarnings = FALSE)
    }
    .assert(dir.exists(next_path),
            paste0("Could not create validation directory: ", next_path))
    resolved <- normalizePath(next_path, winslash = "/", mustWork = TRUE)
    .assert(.inside(resolved, root, must_work = TRUE),
            "A validation output directory resolves outside the package.")
    current <- resolved
  }
  current
}

.env_value <- function(primary, fallback = NULL, default = NULL,
                       required = FALSE) {
  value <- Sys.getenv(primary, unset = "")
  if (!nzchar(value) && !is.null(fallback))
    value <- Sys.getenv(fallback, unset = "")
  if (!nzchar(value) && !is.null(default)) value <- default
  if (required && !nzchar(value))
    stop("Missing required environment variable `", primary, "`.",
         call. = FALSE)
  value
}

.env_number <- function(primary, fallback = NULL, default = NULL,
                        required = FALSE) {
  value <- .env_value(primary, fallback, default, required)
  out <- suppressWarnings(as.numeric(value))
  if (length(out) != 1L || !is.finite(out))
    stop("Environment variable `", primary, "` must be a finite number.",
         call. = FALSE)
  out
}

.write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
  invisible(path)
}

.checked_tarball_identity <- function(workspace_root, source_version) {
  raw_path <- Sys.getenv("OCULAR_CHECKED_TARBALL", unset = "")
  .assert(nzchar(raw_path),
          paste0("OCULAR_CHECKED_TARBALL is not set; follow the build and ",
                 "install sequence in tests/manual/README.md."))
  .assert(!nzchar(Sys.readlink(raw_path)),
          "OCULAR_CHECKED_TARBALL must not be a symbolic link.")
  path <- normalizePath(raw_path, winslash = "/", mustWork = TRUE)
  .assert(isTRUE(file_test("-f", path)),
          "OCULAR_CHECKED_TARBALL must identify a regular file.")
  .assert(.inside(path, workspace_root, must_work = TRUE),
          "The checked tarball must be inside the ocular workspace.")
  expected_name <- paste0("ocular_", source_version, ".tar.gz")
  .assert(identical(basename(path), expected_name),
          paste0("The checked tarball must be named ", expected_name, "."))
  md5 <- unname(tools::md5sum(path))
  .assert(length(md5) == 1L && !is.na(md5) && nzchar(md5),
          "Could not calculate the checked tarball checksum.")
  expected_md5 <- Sys.getenv("OCULAR_CHECKED_TARBALL_MD5", unset = "")
  .assert(nzchar(expected_md5) && identical(md5, expected_md5),
          paste0("The checked tarball checksum does not match the value ",
                 "recorded immediately after the build."))
  list(path = path, md5 = md5)
}

.verify_install_receipt <- function(package_root, checked_tarball,
                                    installed_version, installed_path,
                                    validation_library) {
  raw_path <- Sys.getenv("OCULAR_INSTALL_RECEIPT", unset = "")
  .assert(nzchar(raw_path),
          paste0("OCULAR_INSTALL_RECEIPT is not set; create the receipt ",
                 "immediately after installation as shown in ",
                 "tests/manual/README.md."))
  .assert(!nzchar(Sys.readlink(raw_path)),
          "OCULAR_INSTALL_RECEIPT must not be a symbolic link.")
  path <- normalizePath(raw_path, winslash = "/", mustWork = TRUE)
  .assert(isTRUE(file_test("-f", path)) &&
            .inside(path, package_root, must_work = TRUE),
          "The install receipt must be a regular file inside the package.")
  .assert(identical(dirname(path), dirname(validation_library)),
          paste0("The install receipt must be in the timestamped directory ",
                 "that contains OCULAR_VALIDATION_LIBRARY."))
  receipt <- utils::read.csv(path, stringsAsFactors = FALSE,
                             check.names = FALSE)
  .assert(identical(names(receipt), c("item", "value")) &&
            nrow(receipt) > 0L && !anyNA(receipt$item) &&
            anyDuplicated(receipt$item) == 0L,
          "The install receipt has an invalid structure.")
  receipt_value <- function(item) {
    hit <- which(receipt$item == item)
    .assert(length(hit) == 1L,
            paste0("The install receipt lacks a unique `", item, "` entry."))
    receipt$value[[hit]]
  }
  recorded_tarball <- normalizePath(
    receipt_value("checked_tarball"), winslash = "/", mustWork = TRUE)
  recorded_installed <- normalizePath(
    receipt_value("installed_package_path"), winslash = "/", mustWork = TRUE)
  recorded_library <- normalizePath(
    receipt_value("validation_library"), winslash = "/", mustWork = TRUE)
  .assert(identical(recorded_tarball, checked_tarball$path) &&
            identical(receipt_value("checked_tarball_md5"),
                      checked_tarball$md5) &&
            identical(receipt_value("installed_package_version"),
                      installed_version) &&
            identical(recorded_installed, installed_path) &&
            identical(recorded_library, validation_library),
          paste0("The installed ocular package, validation library, and ",
                 "checked tarball do not match the install receipt."))
  receipt_md5 <- unname(tools::md5sum(path))
  .assert(length(receipt_md5) == 1L && !is.na(receipt_md5) &&
            nzchar(receipt_md5),
          "Could not calculate the install receipt checksum.")
  list(path = path, md5 = receipt_md5)
}

.save_png <- function(path, draw, width = 1600L, height = 1000L) {
  grDevices::png(path, width = width, height = height, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw()
  invisible(path)
}

.finite_values <- function(raster) {
  values <- as.numeric(terra::values(raster, mat = FALSE))
  values[is.finite(values)]
}

.field_mask_raster <- function(field) {
  template <- terra::rast(field$internals$feature_stack[[1L]])
  terra::values(template) <- NA_real_
  cells <- which(field$state$alive_mat, arr.ind = TRUE)
  if (nrow(cells) > 0L) {
    cell_numbers <- terra::cellFromRowCol(
      template, cells[, "row"], cells[, "col"])
    template[cell_numbers] <- 1
  }
  names(template) <- "delineated_field"
  template
}

.series_value_name <- function(series) {
  candidates <- c("value", "union", grep("^area_", names(series), value = TRUE))
  candidates <- candidates[candidates %in% names(series)]
  .assert(length(candidates) > 0L,
          "The time series has no recognised value column.")
  candidates[[1L]]
}

.check_series <- function(series, label) {
  .assert(is.data.frame(series) && "date" %in% names(series) &&
            nrow(series) > 0L,
          paste0(label, " is not a non-empty date series."))
  .assert(!anyNA(as.Date(series$date)),
          paste0(label, " contains an invalid date."))
  value_name <- .series_value_name(series)
  values <- series[[value_name]]
  .assert(is.numeric(values) && any(is.finite(values)),
          paste0(label, " has no finite values."))
  invisible(value_name)
}

.plot_series <- function(series, path, title, fused = FALSE) {
  value_name <- .series_value_name(series)
  dates <- as.Date(series$date)
  values <- series[[value_name]]
  .assert(any(is.finite(values)), "The time series has no finite values.")
  .save_png(path, function() {
    if (isTRUE(fused) && "is_fused" %in% names(series)) {
      colours <- ifelse(series$is_fused, "#D55E00", "#0072B2")
      graphics::plot(dates, values, type = "n", xlab = "Date", ylab = value_name,
                     main = title)
      graphics::lines(dates, values, col = "grey55")
      graphics::points(dates, values, pch = 19, col = colours)
      graphics::legend("topleft", legend = c("Landsat observation", "Estimate"),
                       pch = 19, col = c("#0072B2", "#D55E00"), bty = "n")
    } else {
      graphics::plot(dates, values, type = "o", pch = 19,
                     xlab = "Date", ylab = value_name, main = title)
    }
  })
}

.check_metric_range <- function(x, label) {
  .assert(is.numeric(x) && length(x) == 1L && is.finite(x) &&
            x >= 0 && x <= 1,
          paste0(label, " must be one finite value in [0, 1]."))
}

.check_count <- function(x, label) {
  .assert(is.numeric(x) && length(x) == 1L && is.finite(x) &&
            x >= 0 && x == floor(x),
          paste0(label, " must be one non-negative integer count."))
}

.run <- function() {
  old_options <- options(warn = 2)
  on.exit(options(old_options), add = TRUE)
  package_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  description_path <- file.path(package_root, "DESCRIPTION")
  .assert(file.exists(description_path),
          "Run this script from the ocular package root (the directory containing DESCRIPTION).")
  description <- read.dcf(description_path)
  .assert(identical(unname(description[1L, "Package"]), "ocular"),
          "The working directory is not the ocular package root.")
  workspace_root <- dirname(package_root)
  source_version <- unname(description[1L, "Version"])
  session_temp <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  .assert(.inside(session_temp, package_root, must_work = TRUE),
          paste0("R's session tempdir is outside the ocular package. ",
                 "Set TMPDIR to an existing package-local directory before ",
                 "starting R."))

  allowed_cases <- c("intro", "mcd43a4", "fusion", "ftw-local", "ftw-remote")
  requested <- commandArgs(trailingOnly = TRUE)
  requested <- requested[requested != "--args"]
  if (length(requested) == 0L) {
    requested <- allowed_cases
  }
  requested <- unique(requested[nzchar(requested)])
  .assert(length(requested) > 0L, "No live validation cases were requested.")
  unknown <- setdiff(requested, allowed_cases)
  .assert(length(unknown) == 0L,
          paste0("Unknown validation case(s): ", paste(unknown, collapse = ", "),
                 ". Allowed cases: ", paste(allowed_cases, collapse = ", "), "."))

  run_id <- paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), "-", Sys.getpid())
  artifacts_root <- .ensure_directory(
    file.path(package_root, ".check-artifacts"), package_root)
  validation_root <- .ensure_directory(
    file.path(artifacts_root, "live-release-validation"), package_root)
  output_candidate <- file.path(validation_root, run_id)
  .assert(!file.exists(output_candidate) && !dir.exists(output_candidate) &&
            !nzchar(Sys.readlink(output_candidate)),
          paste0("Refusing to reuse an existing validation run: ",
                 output_candidate))
  output_root <- .ensure_directory(
    output_candidate, package_root)
  cache_root <- .ensure_directory(file.path(output_root, "cache"), package_root)
  terra_temp <- .ensure_directory(
    file.path(output_root, "terra-temp"), package_root)
  Sys.setenv(R_USER_CACHE_DIR = cache_root)

  configured_sources <- c(
    Sys.getenv("OCULAR_FTW_LOCAL_SOURCE", unset = ""),
    Sys.getenv("OCULAR_FTW_REMOTE_SOURCE", unset = ""))
  configured_sources <- configured_sources[nzchar(configured_sources)]
  redact <- function(message) {
    for (source in configured_sources)
      message <- gsub(source, "<configured field boundary source>", message,
                      fixed = TRUE)
    message <- gsub("(https://|s3://)[^[:space:]\"']+",
                    "<remote URL>", message, ignore.case = TRUE)
    message
  }

  log_connection <- file(file.path(output_root, "run.log"), open = "wt")
  on.exit(close(log_connection), add = TRUE)
  log_line <- function(...) {
    line <- redact(paste0(..., collapse = ""))
    cat(line, "\n", sep = "")
    writeLines(line, log_connection)
    flush(log_connection)
    invisible(line)
  }

  setup_messages <- character()
  setup_error <- NULL
  checked_tarball <- NULL
  install_receipt <- NULL
  installed_version <- NA_character_
  validation_library <- NA_character_
  installed_path <- NA_character_
  setup_ok <- tryCatch(
    withCallingHandlers({
      required_packages <- c("ocular", "sf", "terra")
      missing_packages <- required_packages[!vapply(
        required_packages, requireNamespace, logical(1L), quietly = TRUE)]
      .assert(length(missing_packages) == 0L,
              paste0("Install the required package(s) before running this script: ",
                     paste(missing_packages, collapse = ", "), "."))
      installed_version <- as.character(utils::packageVersion("ocular"))
      .assert(identical(installed_version, source_version),
              paste0("Installed ocular version ", installed_version,
                     " does not match source version ", source_version,
                     ". Reinstall the current source before live validation."))
      validation_library <- .env_value(
        "OCULAR_VALIDATION_LIBRARY", required = TRUE)
      validation_library <- normalizePath(
        validation_library, winslash = "/", mustWork = TRUE)
      .assert(dir.exists(validation_library) &&
                .inside(validation_library, package_root, must_work = TRUE),
              paste0("OCULAR_VALIDATION_LIBRARY must be an existing ",
                     "directory inside the ocular package."))
      installed_path <- normalizePath(
        find.package("ocular"), winslash = "/", mustWork = TRUE)
      .assert(identical(dirname(installed_path), validation_library),
              paste0("The loaded ocular package is not from ",
                     "OCULAR_VALIDATION_LIBRARY. Use the isolated library ",
                     "created from the checked tarball."))
      checked_tarball <- .checked_tarball_identity(
        workspace_root, source_version)
      install_receipt <- .verify_install_receipt(
        package_root, checked_tarball, installed_version, installed_path,
        validation_library)
      terra::terraOptions(tempdir = terra_temp)
      sf::sf_use_s2(TRUE)
      TRUE
    }, message = function(message) {
      setup_messages <<- c(
        setup_messages, redact(conditionMessage(message)))
      invokeRestart("muffleMessage")
    }),
    error = function(error) {
      setup_error <<- redact(conditionMessage(error))
      FALSE
    }
  )
  setup_status <- if (isTRUE(setup_ok)) "PASS" else "FAIL"
  setup_detail <- if (isTRUE(setup_ok)) "completed" else setup_error
  .write_csv(data.frame(
    stage = "setup", status = setup_status, detail = setup_detail,
    messages = paste(unique(setup_messages), collapse = " | "),
    stringsAsFactors = FALSE),
    file.path(output_root, "setup-manifest.csv"))
  if (isTRUE(setup_ok)) {
    .write_csv(data.frame(
      item = c(
        "source_package_version", "checked_tarball",
        "checked_tarball_md5", "install_receipt", "install_receipt_md5",
        "installed_package_version",
        "installed_package_path", "validation_library",
        "active_library_paths", "dependency_isolation_scope"
      ),
      value = c(
        source_version, checked_tarball$path, checked_tarball$md5,
        install_receipt$path, install_receipt$md5, installed_version,
        installed_path, validation_library,
        paste(.libPaths(), collapse = "; "),
        paste0("ocular is required to load from validation_library; ",
               "dependencies may resolve from other active library paths")
      ),
      stringsAsFactors = FALSE
    ), file.path(output_root, "install-identity.csv"))
  }
  if (length(setup_messages) > 0L)
    for (line in unique(setup_messages)) log_line("[setup] message: ", line)
  log_line("[setup] ", setup_status, ": ", setup_detail)
  if (!isTRUE(setup_ok))
    stop(setup_detail, call. = FALSE)

  results <- list()
  values <- new.env(parent = emptyenv())
  run_case <- function(name, code) {
    started <- Sys.time()
    case_warnings <- character()
    case_messages <- character()
    case_output <- character()
    case_error <- NULL
    value <- NULL
    log_line("", "[", name, "] starting")
    result <- tryCatch({
      case_output <- utils::capture.output({
        value <- tryCatch(
          withCallingHandlers(
            code(),
            warning = function(warning) {
              case_warnings <<- c(
                case_warnings, redact(conditionMessage(warning)))
              invokeRestart("muffleWarning")
            },
            message = function(message) {
              case_messages <<- c(
                case_messages, redact(conditionMessage(message)))
              invokeRestart("muffleMessage")
            }
          ),
          error = function(error) {
            case_error <<- error
            NULL
          }
        )
      }, type = "output")
      if (!is.null(case_error)) stop(case_error)
      assign(name, value, envir = values)
      if (length(case_warnings) > 0L) {
        list(status = "WARN",
             detail = paste0("completed with ", length(case_warnings),
                             " warning(s)"))
      } else {
        list(status = "PASS", detail = "completed")
      }
    }, error = function(error) {
      list(status = "FAIL", detail = redact(conditionMessage(error)))
    })
    if (length(case_output) > 0L)
      for (line in redact(case_output)) log_line("[", name, "] output: ", line)
    if (length(case_messages) > 0L)
      for (line in unique(case_messages)) log_line("[", name, "] message: ", line)
    if (length(case_warnings) > 0L)
      for (line in unique(case_warnings)) log_line("[", name, "] warning: ", line)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    log_line("[", name, "] ", result$status, ": ", result$detail)
    results[[length(results) + 1L]] <<- data.frame(
      case = name, status = result$status, seconds = round(elapsed, 2),
      detail = result$detail, warning_count = length(case_warnings),
      warnings = paste(unique(case_warnings), collapse = " | "),
      stringsAsFactors = FALSE)
    invisible(result$status)
  }

  intro_query <- list(
    longitude = -96.80342,
    latitude = 33.26425,
    start_date = "2005-03-15",
    end_date = "2005-09-05",
    x_metres = 1000,
    index_name = "NDVI",
    source = "landsat-5"
  )
  fusion_query <- list(
    longitude = 147.227906,
    latitude = -33.0675167,
    start_date = "2016-05-18",
    end_date = "2016-09-21",
    x_metres = 1000,
    index_name = "NDVI",
    source = "landsat-8"
  )

  if ("intro" %in% requested) run_case("intro", function() {
    rs <- do.call(ocular::get_rs, intro_query)
    .assert(ocular::is_rs(rs), "get_rs() did not return an ocular object.")
    .assert(identical(rs$geom$source, "landsat-5"),
            "The introductory query did not retain the Landsat 5/7 source.")
    .assert(length(rs$scenes) > 0L, "The introductory query returned no scenes.")
    scene_dates <- as.Date(vapply(rs$scenes, function(scene)
      as.character(scene$date), character(1L)))
    .assert(!anyNA(scene_dates), "One or more introductory scene dates are invalid.")
    .assert(inherits(rs$internals$feature_stack, "SpatRaster"),
            "The introductory query has no feature stack.")

    field <- ocular::boundary_delineation(rs)
    alive <- field$state$alive_mat
    .assert(is.matrix(alive) && is.logical(alive),
            "Delineation did not return a logical field mask.")
    .assert(identical(dim(alive), c(field$geom$nr, field$geom$nc)),
            "The delineated mask dimensions do not match the feature grid.")
    .assert(any(alive),
            "The representative introductory delineation returned an empty mask.")

    raster <- ocular::as_raster(field)
    series <- ocular::as_time_series(field)
    .assert(inherits(raster, "SpatRaster") && terra::nlyr(raster) > 0L,
            "as_raster() did not return a populated SpatRaster.")
    .assert(length(.finite_values(raster)) > 0L,
            "The introductory raster contains no finite values.")
    .assert(identical(attr(raster, "rs_meta")$index_name, "NDVI"),
            "The introductory raster metadata does not record NDVI.")
    .assert("value" %in% names(series),
            "The introductory series lacks its documented value column.")
    .check_series(series, "The introductory time series")

    mask <- .field_mask_raster(field)
    terra::writeRaster(raster, file.path(output_root, "intro-delineated-raster.tif"))
    terra::writeRaster(mask, file.path(output_root, "intro-field-mask.tif"))
    .write_csv(series, file.path(output_root, "intro-series.csv"))
    .save_png(file.path(output_root, "intro-delineation.png"), function() {
      old <- graphics::par(mfrow = c(1, 2))
      on.exit(graphics::par(old), add = TRUE)
      terra::plot(raster, main = "Delineated NDVI composite")
      terra::plot(mask, main = "Delineated field mask", legend = FALSE)
    })
    .plot_series(series, file.path(output_root, "intro-time-series.png"),
                 "Introductory NDVI time series")
    list(rs = rs, field = field, raster = raster, series = series)
  })

  if ("mcd43a4" %in% requested) run_case("mcd43a4", function() {
    bbox <- ocular::point_to_bbox(
      fusion_query$longitude, fusion_query$latitude,
      fusion_query$x_metres, fusion_query$x_metres)
    scenes <- ocular:::.fetchStac(
      bbox = bbox,
      start_date = "2016-05-18",
      end_date = "2016-05-24",
      index_name = "NDVI",
      max_cloud_cover = 100,
      scl_classes = NULL,
      source = "mcd43a4")
    .assert(length(scenes) > 0L, "The direct MCD43A4 probe returned no scenes.")
    expected_assets <- c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band1")
    .assert(all(vapply(scenes, function(scene)
      all(expected_assets %in% names(scene$bands)), logical(1L))),
      "A live MCD43A4 scene lacks an expected configured reflectance band.")
    .assert(all(vapply(scenes, function(scene)
      all(expected_assets %in% names(scene$item_assets)), logical(1L))),
      "The live MCD43A4 item assets do not expose the configured band names.")

    band_rows <- lapply(seq_along(scenes), function(i) {
      do.call(rbind, lapply(expected_assets, function(asset) {
        raster <- scenes[[i]]$bands[[asset]]
        values <- as.numeric(terra::values(raster, mat = FALSE))
        finite <- values[is.finite(values)]
        .assert(length(finite) > 0L,
                paste0("MCD43A4 ", asset, " has no finite cells."))
        .assert(all(finite >= 0 & finite <= 3.2766 + 1e-8),
                paste0("MCD43A4 ", asset,
                       " falls outside the package's configured scaled range."))
        data.frame(
          date = as.character(scenes[[i]]$date), asset = asset,
          n_cells = length(values), n_finite = length(finite),
          n_missing = sum(!is.finite(values)), min = min(finite), max = max(finite),
          stringsAsFactors = FALSE)
      }))
    })
    band_summary <- do.call(rbind, band_rows)
    index_values <- unlist(lapply(scenes, function(scene)
      .finite_values(scene$index)), use.names = FALSE)
    .assert(length(index_values) > 0L,
            "The direct MCD43A4 NDVI probe has no finite index values.")
    .assert(all(index_values >= -1 - 1e-8 & index_values <= 1 + 1e-8),
            "The direct MCD43A4 probe returned NDVI outside [-1, 1].")
    .write_csv(band_summary, file.path(output_root, "mcd43a4-band-summary.csv"))
    .write_csv(data.frame(
      date = vapply(scenes, function(scene) as.character(scene$date), character(1L)),
      assets = paste(expected_assets, collapse = ","), stringsAsFactors = FALSE),
      file.path(output_root, "mcd43a4-asset-surface.csv"))
    list(scenes = scenes, band_summary = band_summary)
  })

  if ("fusion" %in% requested) run_case("fusion", function() {
    rs <- do.call(ocular::get_rs, fusion_query)
    real_count <- length(rs$scenes)
    detection_before <- rs$internals$detection_scenes
    .assert(is.list(detection_before) && length(detection_before) > 0L,
            "The Landsat query has no detection scenes before fusion.")
    detection_dates_before <- unname(as.Date(vapply(
      detection_before, function(scene) as.character(scene$date), character(1L))))
    .assert(!anyNA(detection_dates_before),
            "A pre-fusion detection scene has an invalid date.")
    fused <- ocular::add_modis(rs)
    is_fused <- vapply(fused$scenes, function(scene) isTRUE(scene$fused), logical(1L))
    .assert(sum(!is_fused) == real_count,
            "Fusion changed the number of real Landsat scenes.")
    .assert(any(is_fused), "add_modis() did not append an estimated scene.")
    real_dates <- as.Date(vapply(fused$scenes[!is_fused], function(scene)
      as.character(scene$date), character(1L)))
    fused_dates <- as.Date(vapply(fused$scenes[is_fused], function(scene)
      as.character(scene$date), character(1L)))
    nearest_days <- vapply(fused_dates, function(date)
      min(abs(as.numeric(date - real_dates))), numeric(1L))
    .assert(all(nearest_days > 3),
            "A fused date lies within three days of a real Landsat scene.")
    .assert(all(vapply(fused$scenes[is_fused], function(scene)
      is.numeric(scene$anchor_days) && length(scene$anchor_days) == 1L &&
        scene$anchor_days > 3, logical(1L))),
      "A fused scene has an invalid anchor_days value.")
    real_template <- fused$scenes[[which(!is_fused)[1L]]]$index
    .assert(all(vapply(fused$scenes[is_fused], function(scene)
      terra::compareGeom(scene$index, real_template, stopOnError = FALSE),
      logical(1L))), "A fused scene does not use the Landsat grid.")
    detection_after <- fused$internals$detection_scenes
    .assert(is.list(detection_after) &&
              length(detection_after) == length(detection_before),
            "Fusion changed the number of Landsat detection scenes.")
    detection_dates_after <- unname(as.Date(vapply(
      detection_after, function(scene) as.character(scene$date), character(1L))))
    .assert(identical(detection_dates_after, detection_dates_before),
            "Fusion changed the Landsat detection-scene dates.")
    detection_is_fused <- vapply(detection_after,
                                 function(scene) isTRUE(scene$fused), logical(1L))
    .assert(!any(detection_is_fused),
            "A fused scene leaked into the delineation feature source.")

    field <- ocular::boundary_delineation(fused)
    .assert(is.matrix(field$state$alive_mat) && any(field$state$alive_mat),
            "The fusion-vignette delineation returned an empty mask.")
    series <- ocular::as_time_series(field)
    required_series <- c("date", "value", "is_fused", "anchor_days",
                         "anchor_mae", "anchor_nmae")
    .assert(all(required_series %in% names(series)),
            "The fused time series lacks one or more documented columns.")
    .check_series(series, "The fused time series")
    .assert(any(series$is_fused) && any(!series$is_fused),
            "The fused series does not contain both observations and estimates.")
    .assert(any(is.finite(series$value[series$is_fused])) &&
              any(is.finite(series$anchor_mae[series$is_fused])),
            "The fused series has no finite estimated values or anchor MAE.")
    .write_csv(series, file.path(output_root, "fusion-series.csv"))
    .plot_series(series, file.path(output_root, "fusion-time-series.png"),
                 "Observed and estimated NDVI", fused = TRUE)

    validation_query <- fusion_query
    names(validation_query)[names(validation_query) == "longitude"] <- "x"
    validation <- do.call(ocular::validate_data_fusion, validation_query)
    .assert(is.data.frame(validation) && nrow(validation) > 0L,
            "validate_data_fusion() returned no rows.")
    required_metrics <- c("date", "rmse", "mae", "bias", "r", "n_pixels")
    .assert(all(required_metrics %in% names(validation)),
            "Fusion validation lacks one or more documented metric columns.")
    .assert(!anyNA(as.Date(validation$date)),
            "Fusion validation contains an invalid date.")
    expected_holdout_dates <- unname(sort(as.Date(vapply(
      rs$scenes, function(scene) as.character(scene$date), character(1L))))
    validation_dates <- unname(sort(as.Date(validation$date)))
    .assert(identical(validation_dates, expected_holdout_dates),
            paste0("Fusion validation did not return one held-out row for ",
                   "every retrieved Landsat scene."))
    .assert(all(is.finite(validation$n_pixels) & validation$n_pixels >= 0 &
                  validation$n_pixels == floor(validation$n_pixels)),
            "Fusion validation n_pixels must be non-negative integer counts.")
    min_valid_n <- 4L
    supported <- validation$n_pixels >= min_valid_n
    .assert(any(supported),
            "No held-out scene meets the default four-pixel support threshold.")
    agreement <- c("rmse", "mae", "bias")
    .assert(all(vapply(validation[agreement], function(metric)
      all(is.finite(metric[supported])), logical(1L))),
      "A supported held-out scene has a missing agreement metric.")
    below_threshold_metrics <- c(agreement, "r")
    .assert(all(vapply(validation[below_threshold_metrics], function(metric)
      all(is.na(metric[!supported])), logical(1L))),
      "A below-threshold held-out scene has a non-missing agreement metric.")
    .assert(all(validation$rmse[supported] >= 0) &&
              all(validation$mae[supported] >= 0),
            "Fusion RMSE or MAE contains a negative supported value.")
    finite_r <- validation$r[is.finite(validation$r)]
    .assert(all(finite_r >= -1 - 1e-8 & finite_r <= 1 + 1e-8),
            "Fusion validation correlation falls outside [-1, 1].")
    .write_csv(validation, file.path(output_root, "fusion-validation.csv"))
    list(rs = fused, field = field, series = series, validation = validation)
  })

  ftw_config <- function(kind) {
    upper <- toupper(kind)
    prefix <- paste0("OCULAR_FTW_", upper, "_")
    source <- .env_value(paste0(prefix, "SOURCE"), required = TRUE)
    if (identical(kind, "local")) {
      source <- normalizePath(source, winslash = "/", mustWork = TRUE)
      .assert(.inside(source, workspace_root),
              "The configured local GeoParquet must be inside the package workspace.")
      .assert(file.exists(source) && !dir.exists(source),
              "The configured local GeoParquet must be a regular file.")
    } else {
      .assert(grepl("^(https://|s3://)", source, ignore.case = TRUE),
              "The remote GeoParquet source must use HTTPS or S3.")
    }
    datetime_col <- .env_value(
      paste0(prefix, "DATETIME_COL"), "OCULAR_FTW_DATETIME_COL",
      "determination_datetime")
    if (toupper(datetime_col) %in% c("NULL", "NONE")) datetime_col <- NULL
    crs <- .env_value(paste0(prefix, "CRS"), "OCULAR_FTW_CRS", "4326")
    if (grepl("^[0-9]+$", crs)) crs <- as.integer(crs)
    list(
      source = source,
      longitude = .env_number(paste0(prefix, "LONGITUDE"),
                              "OCULAR_FTW_LONGITUDE", required = TRUE),
      latitude = .env_number(paste0(prefix, "LATITUDE"),
                             "OCULAR_FTW_LATITUDE", required = TRUE),
      start_date = .env_value(paste0(prefix, "START_DATE"),
                              "OCULAR_FTW_START_DATE", required = TRUE),
      end_date = .env_value(paste0(prefix, "END_DATE"),
                            "OCULAR_FTW_END_DATE", required = TRUE),
      x_metres = .env_number(paste0(prefix, "X_METRES"),
                             "OCULAR_FTW_X_METRES", "1000"),
      index_name = .env_value(paste0(prefix, "INDEX"),
                              "OCULAR_FTW_INDEX", "NDVI"),
      imagery_source = .env_value(paste0(prefix, "IMAGERY_SOURCE"),
                                  "OCULAR_FTW_IMAGERY_SOURCE", "landsat-8"),
      crs = crs,
      geometry_col = .env_value(paste0(prefix, "GEOMETRY_COL"),
                                "OCULAR_FTW_GEOMETRY_COL", "geometry"),
      datetime_col = datetime_col,
      expect_year = .is_true(.env_value(paste0(prefix, "EXPECT_YEAR"),
                                        "OCULAR_FTW_EXPECT_YEAR", "false"))
    )
  }

  run_ftw <- function(kind) {
    config <- ftw_config(kind)
    optional_packages <- c("duckdb", "DBI")
    missing <- optional_packages[!vapply(
      optional_packages, requireNamespace, logical(1L), quietly = TRUE)]
    .assert(length(missing) == 0L,
            paste0("FTW validation needs installed package(s): ",
                   paste(missing, collapse = ", "), "."))
    rs <- ocular::get_rs(
      longitude = config$longitude,
      latitude = config$latitude,
      start_date = config$start_date,
      end_date = config$end_date,
      x_metres = config$x_metres,
      index_name = config$index_name,
      source = config$imagery_source)
    prior <- ocular::add_ftw_prior(
      rs, source = config$source, crs = config$crs, refresh = TRUE,
      geometry_col = config$geometry_col,
      datetime_col = config$datetime_col)
    attached <- prior$geom$ftw_prior
    .assert(is.list(attached) && inherits(attached$polygon, "sf") &&
              nrow(attached$polygon) == 1L &&
              !all(sf::st_is_empty(attached$polygon)) &&
              !is.na(sf::st_crs(attached$polygon)),
            "No usable one-field polygon was attached.")
    seed <- sf::st_sfc(sf::st_point(c(config$longitude, config$latitude)),
                       crs = 4326)
    seed <- sf::st_transform(seed, sf::st_crs(attached$polygon))
    .assert(any(sf::st_intersects(seed, attached$polygon, sparse = FALSE)),
            "The selected field boundary does not contain the configured seed.")
    if (isTRUE(config$expect_year))
      .assert(length(attached$source_year) == 1L &&
                is.numeric(attached$source_year) &&
                is.finite(attached$source_year),
              "The configured dated source did not yield a source year.")
    .assert(identical(prior$internals$calibration$source, "ftw"),
            "The attached polygon did not provide FTW calibration support.")

    field <- ocular::boundary_delineation(prior)
    .assert(is.matrix(field$state$alive_mat) && any(field$state$alive_mat),
            "Field-prior delineation returned an empty mask.")
    diagnostic <- ocular::diagnose_against_ftw(field)
    metric_names <- c("iou", "f1", "precision", "recall")
    count_names <- c("n_pred", "n_ref", "n_intersection", "n_union")
    .assert(all(c(metric_names, count_names) %in% names(diagnostic)),
            "The field boundary prior diagnostic lacks documented fields.")
    for (name in metric_names) .check_metric_range(diagnostic[[name]], name)
    for (name in count_names) .check_count(diagnostic[[name]], name)
    .assert(diagnostic$n_intersection <= diagnostic$n_pred &&
              diagnostic$n_intersection <= diagnostic$n_ref &&
              diagnostic$n_pred <= diagnostic$n_union &&
              diagnostic$n_ref <= diagnostic$n_union,
            "Field-prior diagnostic counts are inconsistent.")

    prefix <- paste0("ftw-", kind)
    diagnostic_row <- as.data.frame(diagnostic, stringsAsFactors = FALSE)
    diagnostic_row$source_kind <- kind
    source_year <- attached$source_year
    if (length(source_year) != 1L || !is.numeric(source_year) ||
        !is.finite(source_year))
      source_year <- NA_integer_
    diagnostic_row$source_year <- source_year
    diagnostic_row <- diagnostic_row[c("source_kind", "source_year",
                                       metric_names, count_names)]
    .write_csv(diagnostic_row,
               file.path(output_root, paste0(prefix, "-diagnostic.csv")))
    series <- ocular::as_time_series(field)
    .check_series(series, paste("The", kind, "field boundary prior time series"))
    .write_csv(series, file.path(output_root, paste0(prefix, "-series.csv")))
    mask <- .field_mask_raster(field)
    prior_vector <- terra::project(terra::vect(attached$polygon), terra::crs(mask))
    .save_png(file.path(output_root, paste0(prefix, "-delineation.png")), function() {
      terra::plot(mask, main = paste("Delineation and", kind, "field prior"),
                  legend = FALSE)
      terra::lines(prior_vector, col = "#D55E00", lwd = 2)
    })
    list(field = field, diagnostic = diagnostic_row)
  }

  if ("ftw-local" %in% requested)
    run_case("ftw-local", function() run_ftw("local"))
  if ("ftw-remote" %in% requested)
    run_case("ftw-remote", function() run_ftw("remote"))

  manifest <- do.call(rbind, results)
  .write_csv(manifest, file.path(output_root, "manifest.csv"))
  writeLines(c(
    paste0("ocular source version: ", source_version),
    paste0("run id: ", run_id),
    paste0("requested cases: ", paste(requested, collapse = ", ")),
    paste0("MCD43A4 quality layers are not applied by ocular ", source_version,
           "."),
    "Fusion metrics are query-specific leave-one-out agreement, not general accuracy.",
    "FTW diagnostics against an attached prior are in-sample, not independent accuracy."
  ), file.path(output_root, "interpretation.txt"))
  writeLines(capture.output(utils::sessionInfo()),
             file.path(output_root, "session-info.txt"))
  cat("\nArtifacts: ", output_root, "\n", sep = "")
  failed <- manifest$case[manifest$status != "PASS"]
  if (length(failed) > 0L)
    stop("Live validation failed: ", paste(failed, collapse = ", "),
         ". Inspect manifest.csv and run.log.", call. = FALSE)
  invisible(output_root)
}

if (.live_enabled(Sys.getenv("OCULAR_RUN_LIVE", unset = "false"))) {
  .run()
} else {
  message("Live validation was not run. Set OCULAR_RUN_LIVE=true explicitly; ",
          "see tests/manual/README.md.")
}
