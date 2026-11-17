.make_failed_refresh_rs <- function(){

  rs <- .make_test_rs()
  centre <- terra::xyFromCell(rs$internals$feature_stack, 221L)
  polygon <- terra::as.polygons(
    terra::ext(centre[1L] - 1, centre[1L] + 1,
               centre[2L] - 1, centre[2L] + 1),
    crs = terra::crs(rs$internals$feature_stack))
  rs$geom$ftw_prior <- list(polygon = polygon, source_year = 2021L)
  rs$internals$calibration$schema_version <- 0L
  return(rs)
}

.capture_refresh <- function(expr){

  warnings <- character()
  value <- withCallingHandlers(expr,
    warning = function(w){
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m){
      invokeRestart("muffleMessage")
    })
  return(list(value = value, warnings = warnings))
}

test_that("one boundary invocation attempts an identical failed refresh once", {
  rs <- .make_failed_refresh_rs()
  original <- ocular:::.calibrateScene
  attempts <- 0L
  testthat::local_mocked_bindings(
    .calibrateScene = function(rs){
      attempts <<- attempts + 1L
      return(original(rs))
    },
    .package = "ocular")

  out <- .capture_refresh(boundary_delineation(rs))

  expect_identical(attempts, 1L)
  expect_identical(out$warnings, paste0(
    "Recalibration for this call failed. The prior calibration was preserved: ",
    "no support (feature data unavailable)"))
  expect_identical(out$value$internals$calibration, rs$internals$calibration)
  expect_identical(out$value$params, rs$params)
  expect_equal(terra::values(out$value$internals$feature_stack),
               terra::values(rs$internals$feature_stack))
  expect_null(attr(out$value, ".calibration_context", exact = TRUE))
  expect_null(attr(rs, ".calibration_context", exact = TRUE))
  expect_true(any(out$value$state$alive_mat))

  again <- .capture_refresh(boundary_delineation(out$value))
  expect_identical(attempts, 2L)
  expect_identical(again$warnings, out$warnings)
})

test_that("nested stages share failures and standalone calls retry", {
  rs <- .make_failed_refresh_rs()
  original <- ocular:::.calibrateScene
  attempts <- 0L
  testthat::local_mocked_bindings(
    .calibrateScene = function(rs){
      attempts <<- attempts + 1L
      return(original(rs))
    },
    .package = "ocular")

  for( stage in list(split_area, segment_interior) ){
    previous <- attempts
    out <- .capture_refresh(stage(rs))
    expect_identical(attempts - previous, 1L)
    expect_length(out$warnings, 1L)
    expect_null(attr(out$value, ".calibration_context", exact = TRUE))
  }

  pipeline <- build_pipeline(
    boundary_delineation() |> split_area() |> segment_interior())
  for( iteration in seq_len(2L) ){
    previous <- attempts
    out <- .capture_refresh(pipeline(rs))
    expect_identical(attempts - previous, 1L)
    expect_length(out$warnings, 1L)
    expect_null(attr(out$value, ".calibration_context", exact = TRUE))
  }

  previous <- attempts
  .capture_refresh(segment_area(rs))
  .capture_refresh(segment_area(rs))
  expect_identical(attempts - previous, 2L)
})

test_that("changed calibration inputs retry within the same invocation", {
  rs <- .make_failed_refresh_rs()
  invocation <- ocular:::.beginCalibrationInvocation(rs)
  rs <- invocation$rs
  attempts <- 0L
  testthat::local_mocked_bindings(
    .calibrateScene = function(rs){
      attempts <<- attempts + 1L
      return(list(error = "unavailable test support"))
    },
    .package = "ocular")
  apply <- function(x, params = x$params){
    return(.capture_refresh(ocular:::.applyScenePriors(x, params)))
  }
  expect_length(apply(rs)$warnings, 1L)
  expect_length(apply(rs)$warnings, 0L)

  changes <- list(
    strictness = function(x){
      x$params$strictness <- "tight"
      return(x)
    },
    threshold = function(x){
      x$params$area_threshold <- c(0.2, 0.9)
      return(x)
    },
    index = function(x){
      x$params$search_index <- "NDVI"
      return(x)
    },
    feature_values = function(x){
      terra::values(x$internals$feature_stack)[1L, 1L] <- 0.2
      return(x)
    },
    feature_extent = function(x){
      terra::ext(x$internals$feature_stack) <-
        terra::ext(x$internals$feature_stack) + 1
      return(x)
    },
    feature_crs = function(x){
      terra::crs(x$internals$feature_stack) <- "EPSG:32749"
      return(x)
    },
    centre = function(x){
      x$geom$centre_rc <- c(10L, 11L)
      return(x)
    },
    point = function(x){
      x$spec$point$longitude <- 117.9
      return(x)
    },
    polygon = function(x){
      x$geom$ftw_prior$polygon <- terra::shift(x$geom$ftw_prior$polygon, dx = 10)
      return(x)
    },
    calibration = function(x){
      x$internals$calibration$baselines$area_sensitivity <- 0.123
      return(x)
    })
  for( change in changes ){
    altered <- change(rs)
    previous <- attempts
    expect_length(apply(altered)$warnings, 1L)
    expect_length(apply(altered)$warnings, 0L)
    expect_identical(attempts - previous, 1L)
  }
  expect_length(apply(rs)$warnings, 0L)
  expect_identical(attempts, 1L + length(changes))

  irrelevant <- rs
  irrelevant$state$alive_mat <- matrix(TRUE, 21L, 21L)
  irrelevant$internals$feat_array <- array(0, c(21L, 21L, 2L))
  irrelevant$geom$ftw_prior$source_year <- 2020L
  irrelevant$params$perimeter_margins <- c(0.1, 0.9)
  expect_length(apply(irrelevant)$warnings, 0L)
  expect_identical(attempts, 1L + length(changes))
})

test_that("caller stages can change features and retry a failed refresh", {
  rs <- .make_failed_refresh_rs()
  original <- ocular:::.calibrateScene
  original_trace <- ocular::trace_perimeter
  attempts <- 0L
  testthat::local_mocked_bindings(
    .calibrateScene = function(rs){
      attempts <<- attempts + 1L
      return(original(rs))
    },
    .package = "ocular")
  trace_perimeter <- function(x, ...){
    terra::values(x$internals$feature_stack)[1L, 1L] <- 0.2
    return(original_trace(x, ...))
  }
  pipeline <- build_pipeline(segment_area() |> trace_perimeter())

  out <- .capture_refresh(pipeline(rs))

  expect_identical(attempts, 2L)
  expect_length(out$warnings, 2L)
  expect_null(attr(out$value, ".calibration_context", exact = TRUE))
})

test_that("errors discard failure context even for a captured intermediate object", {
  rs <- .make_failed_refresh_rs()
  captured <- NULL
  trace_perimeter <- function(x, ...){
    captured <<- x
    stop("caller stage failed")
  }
  pipeline <- build_pipeline(segment_area() |> trace_perimeter())
  expect_error(.capture_refresh(pipeline(rs)), "caller stage failed")
  context <- attr(captured, ".calibration_context", exact = TRUE)
  expect_false(context$active)
  expect_null(context$failures)

  out <- .capture_refresh(split_area(captured))
  expect_length(out$warnings, 1L)
  expect_null(attr(out$value, ".calibration_context", exact = TRUE))
})

test_that("changed support can succeed without replacing the prior on failure", {
  rs <- .make_failed_refresh_rs()
  invocation <- ocular:::.beginCalibrationInvocation(rs)
  rs <- invocation$rs
  failed <- .capture_refresh(ocular:::.applyScenePriors(rs, rs$params))
  expect_length(failed$warnings, 1L)
  expect_identical(failed$value$internals$calibration, rs$internals$calibration)
  rs$geom$ftw_prior <- NULL
  refreshed <- .capture_refresh(ocular:::.applyScenePriors(rs, rs$params))
  expect_length(refreshed$warnings, 0L)
  expect_identical(refreshed$value$internals$calibration$schema_version,
                   ocular:::.CALIBRATION_SCHEMA_VERSION)
  expect_identical(refreshed$value$params, rs$params)
})

test_that("thrown calibration errors are scoped and unreadable snapshots still retry", {
  rs <- .make_failed_refresh_rs()
  invocation <- ocular:::.beginCalibrationInvocation(rs)
  rs <- invocation$rs
  attempts <- 0L
  testthat::local_mocked_bindings(
    .calibrateScene = function(rs){
      attempts <<- attempts + 1L
      stop("test calibration error")
    },
    .package = "ocular")
  for( iteration in seq_len(2L) ){
    .capture_refresh(ocular:::.applyScenePriors(rs, rs$params))
  }
  expect_identical(attempts, 1L)
  testthat::local_mocked_bindings(
    .calibrationInputSnapshot = function(...){
      stop("unreadable feature state")
    },
    .package = "ocular")
  for( iteration in seq_len(2L) ){
    out <- .capture_refresh(ocular:::.applyScenePriors(rs, rs$params))
    expect_identical(out$warnings, paste0(
      "Recalibration for this call failed. The prior calibration was preserved: ",
      "test calibration error"))
  }
  expect_identical(attempts, 3L)
})

test_that("current calibration incurs no input snapshot", {
  rs <- .make_test_rs()
  invocation <- ocular:::.beginCalibrationInvocation(rs)
  testthat::local_mocked_bindings(
    .calibrationInputSnapshot = function(...){
      stop("snapshot should not run")
    },
    .package = "ocular")
  out <- ocular:::.applyScenePriors(invocation$rs, rs$params)
  expect_identical(out$internals$calibration, rs$internals$calibration)
})


test_that("interrupts and warnings promoted to errors close the invocation", {
  rs <- .make_failed_refresh_rs()
  captured <- NULL
  testthat::local_mocked_bindings(
    .calibrateScene = function(rs){
      captured <<- rs
      signalCondition(structure(list(message = "test interrupt", call = NULL),
                                class = c("interrupt", "condition")))
    },
    .package = "ocular")
  interrupted <- tryCatch(boundary_delineation(rs),
                           interrupt = function(e){
                             return("interrupted")
                           })
  expect_identical(interrupted, "interrupted")
  context <- attr(captured, ".calibration_context", exact = TRUE)
  expect_false(context$active)
  expect_null(context$failures)

  testthat::local_mocked_bindings(
    .calibrateScene = function(rs){
      captured <<- rs
      return(list(error = "unavailable test support"))
    },
    .package = "ocular")
  previous_warn <- getOption("warn")
  run_with_warning_errors <- function(){

    previous_options <- options(warn = 2L)
    on.exit(options(previous_options), add = TRUE)
    return(boundary_delineation(rs))
  }
  expect_error(run_with_warning_errors(), "Recalibration for this call failed")
  expect_identical(getOption("warn"), previous_warn)
  context <- attr(captured, ".calibration_context", exact = TRUE)
  expect_false(context$active)
  expect_null(context$failures)
})
