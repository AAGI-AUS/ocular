.deferred_test_scenes <- function(dates = c("2021-01-05", "2021-01-20")){

  rs <- .make_test_rs()
  scenes <- lapply(seq_along(dates), function(i){
    scene <- rs$scenes[[1L]]
    scene$date <- as.Date(dates[[i]])
    scene$index <- scene$index + (i - 1L) * 0.02
    scene
  })
  return(scenes)
}

test_that("retrieval retains observations without building detection features", {
  scenes <- .deferred_test_scenes()
  n_fetch <- 0L
  testthat::local_mocked_bindings(
    .fetchStac = function(...){
      n_fetch <<- n_fetch + 1L
      scenes
    }, .package = "ocular")

  rs <- get_rs(117.8, -32, "2021-01-01", "2021-02-01",
               source = "landsat-8", x_metres = 210)
  expect_null(rs$internals$feature_stack)
  expect_null(rs$internals$calibration)
  expect_null(rs$state$decided_mat)
  expect_identical(rs$scenes, scenes)
  expect_invisible(ocular:::.validateOcular(rs))
  expect_equal(nrow(as_time_series(rs)), 2L)
  expect_s4_class(as_raster(rs), "SpatRaster")

  field <- boundary_delineation(rs)
  expect_s4_class(field$internals$feature_stack, "SpatRaster")
  expect_true(is.list(field$internals$calibration))
  expect_true(any(field$state$alive_mat))
  expect_null(rs$internals$feature_stack)
  expect_identical(n_fetch, 1L)
})

test_that("empty detection windows do not block extraction", {
  scenes <- .deferred_test_scenes(c("2021-01-05", "2021-01-21", "2021-05-14"))
  n_fetch <- 0L
  testthat::local_mocked_bindings(
    .fetchStac = function(...){
      n_fetch <<- n_fetch + 1L
      scenes
    }, .package = "ocular")
  rs <- get_rs(117.8, -32, "2021-01-01", "2021-06-01",
               source = "landsat-8", x_metres = 210)
  expect_equal(nrow(as_time_series(rs)), 3L)
  expect_s4_class(as_raster(rs), "SpatRaster")
  for( stage in list(segment_area, segment_interior, split_area,
                     boundary_delineation) ){
    expect_error(stage(rs), "no scenes in window 2")
  }
  pipeline <- build_pipeline(boundary_delineation() |> as_time_series())
  expect_error(pipeline(rs), "no scenes in window 2")
  expect_null(trace_perimeter(rs)$internals$feature_stack)
  expect_identical(n_fetch, 1L)
})

test_that("alternate detection bands are fetched once and retained for later use", {
  scenes <- .deferred_test_scenes()
  scenes <- lapply(scenes, function(scene){
    scene$bands <- list(nir08 = scene$index, swir22 = scene$index * 0 + 0.2)
    scene$index <- ocular:::.computeIndexFromBands(
      scene$bands, ocular:::landsat_index_list[["NBR"]])
    scene$item_assets <- list(red = "https://fixture.invalid/red.tif")
    scene
  })
  n_fetch <- 0L
  n_assets <- 0L
  testthat::local_mocked_bindings(
    .fetchStac = function(...){
      n_fetch <<- n_fetch + 1L
      scenes
    },
    .fetchAssets = function(...){
      n_assets <<- n_assets + 1L
      list(red = scenes[[1L]]$index * 0 + 0.1)
    }, .package = "ocular")
  rs <- get_rs(117.8, -32, "2021-01-01", "2021-02-01",
               source = "landsat-8", index_name = "NBR", x_metres = 210,
               search_start_date = "2021-01-10", search_end_date = "2021-02-01")
  expect_identical(n_assets, 1L)
  expect_null(rs$internals$feature_stack)
  expect_length(rs$scenes, 2L)
  expect_length(rs$internals$detection_scenes, 1L)
  expect_identical(rs$internals$detection_scenes[[1L]]$date, as.Date("2021-01-20"))
  expect_true("red" %in% names(rs$internals$detection_scenes[[1L]]$bands))
  field <- boundary_delineation(rs)
  expect_true(any(field$state$alive_mat))
  expect_identical(n_fetch, 1L)
  expect_identical(n_assets, 1L)
})

test_that("unavailable detection data stay unavailable until a new retrieval", {
  scenes <- .deferred_test_scenes()
  scenes <- lapply(scenes, function(scene){
    scene$bands <- list(nir08 = scene$index, swir22 = scene$index)
    scene$item_assets <- list(red = "https://fixture.invalid/red.tif")
    scene
  })
  n_assets <- 0L
  testthat::local_mocked_bindings(
    .fetchStac = function(...) scenes,
    .fetchAssets = function(...){
      n_assets <<- n_assets + 1L
      stop("synthetic unavailable detection asset")
    }, .package = "ocular")
  warnings <- character()
  rs <- withCallingHandlers(
    get_rs(117.8, -32, "2021-01-01", "2021-02-01",
            source = "landsat-8", index_name = "NBR", x_metres = 210),
    warning = function(w){
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_length(warnings, 2L)
  expect_true(all(grepl("Could not compose EVI2", warnings, fixed = TRUE)))
  expect_identical(rs$internals$detection_scenes, list())
  expect_equal(nrow(as_time_series(rs)), 2L)
  expect_error(boundary_delineation(rs), "no scenes available")
  expect_identical(n_assets, 2L)
})

test_that("deferred prior calibration preserves self support on both failure paths", {
  scenes <- .deferred_test_scenes()
  testthat::local_mocked_bindings(.fetchStac = function(...) scenes,
                                  .package = "ocular")
  rs <- get_rs(117.8, -32, "2021-01-01", "2021-02-01",
               source = "landsat-8", x_metres = 210)
  prepared <- ocular:::.fetchFeatureStack(rs)
  field <- sf::st_as_sf(terra::project(
    terra::as.polygons(terra::ext(prepared$internals$feature_stack),
                       crs = terra::crs(prepared$internals$feature_stack)),
    "EPSG:4326"))
  field$source_year <- 2020L
  calibrate <- ocular:::.calibrateScene
  failure_mode <- "throw"
  testthat::local_mocked_bindings(
    .ftwCacheRead = function(...) NULL,
    .ftwQuery = function(...) field,
    .ftwCacheWrite = function(...) invisible(NULL),
    .ftwCacheEvictLRU = function(...) invisible(0L),
    .calibrateScene = function(x){
      if( is.null(x$geom$ftw_prior) ){
        return(calibrate(x))
      }
      if( identical(failure_mode, "throw") ){
        stop("synthetic prior calibration failure")
      }
      return(list(source = "ftw", error = "synthetic unavailable prior support"))
    }, .package = "ocular")

  attached <- suppressMessages(add_ftw_prior(rs, source = "local.parquet"))
  expect_null(attached$internals$calibration)
  for( mode in c("throw", "record") ){
    failure_mode <- mode
    expect_warning(out <- segment_area(attached),
                    "FTW calibration was unavailable.*preserved")
    expect_identical(out$internals$calibration, prepared$internals$calibration)
    expect_identical(out$geom$ftw_prior, attached$geom$ftw_prior)
    expect_true(any(out$state$alive_mat))
  }
})

test_that("replacing a deferred prior preserves its calibration when replacement fails", {
  scenes <- .deferred_test_scenes()
  testthat::local_mocked_bindings(.fetchStac = function(...) scenes,
                                  .package = "ocular")
  rs <- get_rs(117.8, -32, "2021-01-01", "2021-02-01",
               source = "landsat-8", x_metres = 210)
  prepared <- ocular:::.fetchFeatureStack(rs)
  field <- sf::st_as_sf(terra::project(
    terra::as.polygons(terra::ext(prepared$internals$feature_stack),
                       crs = terra::crs(prepared$internals$feature_stack)),
    "EPSG:4326"))
  field$source_year <- 2020L
  replacement <- field
  replacement$source_year <- 2019L
  query_result <- field
  calibrate <- ocular:::.calibrateScene
  testthat::local_mocked_bindings(
    .ftwCacheRead = function(...) NULL,
    .ftwQuery = function(...) query_result,
    .ftwCacheWrite = function(...) invisible(NULL),
    .ftwCacheEvictLRU = function(...) invisible(0L),
    .calibrateScene = function(x){
      if( identical(x$geom$ftw_prior$source_year, 2019L) ){
        return(list(source = "ftw", error = "synthetic unavailable replacement support"))
      }
      return(calibrate(x))
    }, .package = "ocular")

  attached <- suppressMessages(add_ftw_prior(rs, source = "first.parquet"))
  expected <- suppressMessages(add_ftw_prior(prepared, source = "first.parquet"))
  expect_identical(expected$internals$calibration$source, "ftw")
  query_result <- replacement
  expect_warning(
    out <- suppressMessages(add_ftw_prior(attached, source = "second.parquet")),
    "FTW calibration was unavailable.*preserved")
  expect_identical(out$internals$calibration, expected$internals$calibration)
  expect_identical(out$geom$ftw_prior$source_year, 2019L)
  expect_identical(segment_area(out)$internals$calibration,
                   expected$internals$calibration)
})
