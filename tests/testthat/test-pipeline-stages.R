test_that("segment_area preserves stored params with call-local strictness", {
  rs <- .make_test_rs()
  before <- rs$params
  out <- suppressMessages(segment_area(rs, strictness = "tight"))
  expect_true(is.matrix(out$state$alive_mat))
  expect_true(is.logical(out$state$alive_mat))
  expect_identical(dim(out$state$alive_mat), c(21L, 21L))
  expect_gt(sum(out$state$alive_mat), 0L)
  expect_identical(out$params, before)
  expect_identical(out$internals$calibration$strictness, "tight")
})

test_that("stage params-list priors and formal priors are equivalent", {
  rs <- .make_test_rs()
  formal <- suppressMessages(segment_area(rs, strictness = "tight"))
  listed <- suppressMessages(segment_area(
    rs, params = list(strictness = "tight")))
  expect_identical(formal$state$alive_mat, listed$state$alive_mat)
  expect_identical(formal$params, listed$params)
  expect_identical(listed$internals$calibration$strictness, "tight")
})

test_that("segment_area consumes a newly recalibrated baseline immediately", {
  rs <- .make_test_rs()
  rs$internals$calibration <- list(
    schema_version = ocular:::.CALIBRATION_SCHEMA_VERSION,
    strictness = "balanced",
    area_threshold = rs$params$area_threshold, index_name = "EVI2",
    baselines = list(area_sensitivity = 0.9))
  seen <- numeric()
  testthat::local_mocked_bindings(
    .calibrateScene = function(x) list(
      schema_version = ocular:::.CALIBRATION_SCHEMA_VERSION,
      source = "self",
      strictness = x$params$strictness,
      area_threshold = x$params$area_threshold, index_name = "EVI2",
      baselines = list(area_sensitivity = 0.05)),
    .searchAreaWithWidens = function(...){
      args <- list(...)
      seen <<- c(seen, args$area_sensitivity)
      matrix(FALSE, args$nr_full, args$nc_full)
    },
    .package = "ocular")
  expect_warning(
    segment_area(rs, strictness = "tight"),
    "no candidate field was identified",
    fixed = TRUE)
  expect_true(length(seen) > 0L)
  expect_equal(seen, rep(0.05, length(seen)))
})

test_that("public delineation stages are quiet on successful calls", {
  area <- expect_silent(segment_area(.make_test_rs()))
  perimeter <- expect_silent(trace_perimeter(
    area, perimeter_margins = c(0.5, 0.9), trace_iter = 2L))
  interior <- expect_silent(segment_interior(perimeter))
  split <- expect_silent(split_area(interior))
  field <- expect_silent(boundary_delineation(
    .make_test_rs(), cleanup_boundary = 0L))

  for (x in list(area, perimeter, interior, split, field))
    expect_true(is_rs(x))
})

test_that("segment_area warns when it returns an empty mask", {
  testthat::local_mocked_bindings(
    .searchAreaWithWidens = function(...){
      args <- list(...)
      matrix(FALSE, args$nr_full, args$nc_full)
    },
    .package = "ocular")

  expect_warning(
    out <- segment_area(.make_test_rs()),
    "no candidate field was identified; the returned mask is empty",
    fixed = TRUE)
  expect_false(any(out$state$alive_mat))
})

test_that("failed recalibration preserves a valid prior record", {
  rs <- .make_test_rs()
  before <- rs$internals$calibration
  testthat::local_mocked_bindings(
    .calibrateScene = function(...) list(error = "no support"),
    .package = "ocular")
  expect_warning(
    out <- ocular:::.applyScenePriors(
      rs, utils::modifyList(rs$params, list(strictness = "tight"))),
    "preserving")
  expect_identical(out$internals$calibration, before)
})

test_that("recalibration reads the full feature stack after a window subset", {
  rs <- .make_test_rs()
  rs$internals$feature_stack <- c(rs$internals$feature_stack,
                                  rs$internals$feature_stack[[1L]] + 0.04)
  names(rs$internals$feature_stack) <- c("w1", "w2", "w3")
  rs$internals$feat_array <- array(0.5, dim = c(rs$geom$nr, rs$geom$nc, 1L))
  seed <- ocular:::.seedNeighborhood(rs)
  expect_identical(ncol(seed), 3L)
})

test_that("public cleanup stages preserve dimensions on a synthetic field", {
  rs <- suppressMessages(segment_area(.make_test_rs()))
  traced <- suppressMessages(trace_perimeter(rs, perimeter_margins = c(0.5, 0.9),
                                              trace_iter = 2L))
  interior <- suppressMessages(segment_interior(traced))
  split <- suppressMessages(split_area(interior))
  for (x in list(traced, interior, split)) {
    expect_true(is_rs(x))
    expect_identical(dim(x$state$alive_mat), c(21L, 21L))
    expect_false(anyNA(x$state$alive_mat))
  }
})

test_that("explicit NULL disables optional stages while omission inherits", {
  rs <- suppressMessages(segment_area(.make_test_rs()))
  interior_calls <- 0L
  split_calls <- 0L
  testthat::local_mocked_bindings(
    .runDetectObjects = function(rs, params, log_msg){
      interior_calls <<- interior_calls + 1L
      rs
    },
    .splitAreaRule = function() list(
      postprocess = function(rs, params, log_msg) rs),
    .runStagedWalk = function(rs, rule, params, log_msg){
      split_calls <<- split_calls + 1L
      rs
    },
    .package = "ocular")

  suppressMessages(segment_interior(rs))
  suppressMessages(split_area(rs))
  expect_identical(interior_calls, 1L)
  expect_identical(split_calls, 1L)

  interior <- suppressMessages(segment_interior(
    rs, interior_sensitivity = NULL,
    params = list(interior_sensitivity = 0.3)))
  split <- suppressMessages(split_area(
    rs, split_sensitivity = NULL,
    params = list(split_sensitivity = 0.3)))
  expect_identical(interior_calls, 1L)
  expect_identical(split_calls, 1L)
  expect_identical(interior$state$alive_mat, rs$state$alive_mat)
  expect_identical(split$state$alive_mat, rs$state$alive_mat)
})

test_that("explicit NULL stage formals retain call ownership", {
  rs <- suppressMessages(segment_area(.make_test_rs()))
  seen <- list()
  testthat::local_mocked_bindings(
    .applyScenePriors = function(rs, effective_params){
      seen[[length(seen) + 1L]] <<-
        attr(effective_params, "user_set", exact = TRUE)
      rs
    },
    .package = "ocular")

  suppressMessages(segment_interior(rs, interior_sensitivity = NULL))
  suppressMessages(split_area(rs, split_sensitivity = NULL))

  expect_true("interior_sensitivity" %in% seen[[1L]])
  expect_true("split_sensitivity" %in% seen[[2L]])
})

test_that("boundary generic sensitivity preserves explicit NULL disables", {
  rs <- .make_test_rs()
  interior_disabled <- logical()
  split_disabled <- logical()
  pass <- function(x, ...) x
  capture_interior <- function(x, interior_sensitivity = 1, ...){
    interior_disabled <<- c(
      interior_disabled,
      !missing(interior_sensitivity) && is.null(interior_sensitivity))
    x
  }
  capture_split <- function(x, split_sensitivity = 1, ...){
    split_disabled <<- c(
      split_disabled,
      !missing(split_sensitivity) && is.null(split_sensitivity))
    x
  }
  testthat::local_mocked_bindings(
    trace_perimeter = pass,
    segment_area = pass,
    segment_interior = capture_interior,
    split_area = capture_split,
    .package = "ocular")

  boundary_delineation(
    rs, vi_sensitivity = 0.2, interior_sensitivity = NULL,
    split_sensitivity = NULL, cleanup_boundary = 0L)

  expect_true(length(interior_disabled) > 0L)
  expect_true(length(split_disabled) > 0L)
  expect_true(all(interior_disabled))
  expect_true(all(split_disabled))
})

test_that("standalone trace_perimeter treats zero as background", {
  r <- terra::rast(nrows = 9, ncols = 9, xmin = 0, xmax = 9,
                   ymin = 0, ymax = 9, crs = "EPSG:32750")
  vals <- matrix(0, 9, 9)
  vals[3:7, 3:7] <- 1
  terra::values(r) <- as.vector(t(vals))
  out <- trace_perimeter(r, perimeter_margins = c(0, 1), trace_iter = 1L)
  out_vals <- terra::values(out, mat = FALSE)
  expect_true(any(is.na(out_vals)))
  expect_lte(sum(!is.na(out_vals)), 25L)
  expect_error(trace_perimeter(r, perimeter_margins = c(0, 1),
                               trace_iter = 1.5), "positive integer")
})

test_that("build_pipeline captures a reusable idiomatic stage chain", {
  bucket <- "monthly"
  pipeline <- build_pipeline(
    as_time_series(time_aggregate = bucket, aggregate_function = "mean"))
  bucket <- "yearly"
  expect_s3_class(pipeline, "ocular_pipeline")
  expect_identical(class(pipeline), c("ocular_pipeline", "function"))
  out <- pipeline(.make_test_rs(field = FALSE))
  expect_equal(nrow(out), 1L)
  expect_identical(out$date, "2021-01")
  expect_error(pipeline(list()), "ocular object")
  expect_error(build_pipeline(stats::median()), "chain of stage")
})

test_that("a captured pipeline retains an explicit NULL stage disable", {
  rs <- suppressMessages(segment_area(.make_test_rs()))
  interior_calls <- 0L
  testthat::local_mocked_bindings(
    .runDetectObjects = function(rs, params, log_msg){
      interior_calls <<- interior_calls + 1L
      rs
    },
    .package = "ocular")
  pipeline <- build_pipeline(
    segment_interior(interior_sensitivity = NULL))

  pipeline(rs)

  expect_identical(interior_calls, 0L)
})

test_that("the real delineation pipeline works across synthetic grid sizes", {
  for (n in c(9L, 15L, 21L, 31L)) {
    rs <- .make_test_rs(n = n)
    ## An empty result is covered by a dedicated warning test above.
    field <- suppressWarnings(
      boundary_delineation(rs, cleanup_boundary = 0L))
    expect_true(is_rs(field), info = paste("grid", n))
    expect_identical(dim(field$state$alive_mat), c(n, n),
                     info = paste("grid", n))
    expect_false(anyNA(field$state$alive_mat), info = paste("grid", n))

    raster <- as_raster(field)
    daily <- as_time_series(field)
    expect_s4_class(raster, "SpatRaster")
    expect_s3_class(daily, "data.frame")
    if (sum(field$state$alive_mat) == 0L) {
      expect_true(all(is.na(terra::values(raster))), info = paste("grid", n))
      expect_equal(nrow(daily), 0L, info = paste("grid", n))
    } else {
      expect_equal(nrow(daily), length(rs$scenes), info = paste("grid", n))
    }
  }
})

test_that("boundary_delineation rejects ambiguous or unknown arguments early", {
  expect_error(boundary_delineation(.make_test_rs(), cleanup_boundary = 1.5),
               "integer")
  expect_error(boundary_delineation(.make_test_rs(), unknown_control = 1),
               "unknown parameter")
  expect_error(boundary_delineation(117.8, unknown_control = 1),
               "unknown argument")
})

test_that("closing boundary passes retain the whole-pipeline params overlay", {
  rs <- .make_test_rs()
  trace_params <- list()
  pass <- function(x, ...) x
  trace <- function(x, params = NULL, ...){
    trace_params[[length(trace_params) + 1L]] <<- params
    x
  }
  testthat::local_mocked_bindings(
    trace_perimeter = trace, segment_area = pass, split_area = pass,
    segment_interior = pass, .package = "ocular")
  boundary_delineation(rs, params = list(multiple_areas = TRUE),
                       cleanup_boundary = 0L)
  expect_true(isTRUE(trace_params[[length(trace_params) - 1L]]$multiple_areas))
  expect_true(isTRUE(trace_params[[length(trace_params)]]$multiple_areas))
})
