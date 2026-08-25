test_that("additive fusion applies the coarse delta and physical bounds", {
  make_r <- function(values) {
    r <- terra::rast(nrows = 1, ncols = length(values), xmin = 0,
                     xmax = length(values), ymin = 0, ymax = 1,
                     crs = "EPSG:32750")
    terra::setValues(r, values)
  }
  fine <- list(nir08 = make_r(c(0.4, 0.95)),
               red = make_r(c(0.2, 0.1)))
  c1 <- list(make_r(c(0.4, 0.95)), make_r(c(0.2, 0.1)))
  c2 <- list(make_r(c(0.5, 1.2)), make_r(c(0.25, -0.2)))
  delta <- ocular:::.fuseDelta(fine, c1, c2)
  expect_equal(as.numeric(terra::values(delta[[1L]])[, 1L]), c(0.5, 1),
               tolerance = 1e-7)
  expect_equal(as.numeric(terra::values(delta[[2L]])[, 1L]), c(0.25, 0),
               tolerance = 1e-7)

  idx_fun <- function(r) (r$nir08 - r$red) / (r$nir08 + r$red)
  fused <- ocular:::.fuseLandsatMODIS(fine, c1, c2, idx_fun)
  expect_true(inherits(fused$index, "SpatRaster"))
  expect_true(all(is.finite(terra::values(fused$index))))
})

test_that("agreement metrics have the documented signs", {
  m <- ocular:::.agreementMetrics(c(1, 2, 3), c(2, 3, 4))
  expect_equal(m$rmse, 1)
  expect_equal(m$mae, 1)
  expect_equal(m$bias, 1)
  expect_equal(m$r, 1)
  expect_true(is.na(ocular:::.agreementMetrics(1, 1)$r))
})

test_that("fusion validation rejects invalid counts and non-Landsat inputs", {
  expect_error(validate_data_fusion(1, min_valid_n = 1.5),
               "positive integer")
  x <- data.frame(date = "2021-01-01", value = 0.5)
  attr(x, "rs_meta") <- list(
    longitude = 117.8, latitude = -32, start_date = "2021-01-01",
    end_date = "2021-02-01", index_name = "EVI2",
    source = "sentinel-2", x_metres = 100, y_metres = 100,
    max_cloud_cover = 50, scl_classes = NULL, params = rs_params())
  expect_error(validate_data_fusion(x), "must originate from Landsat")

  y <- x[, "value", drop = FALSE]
  attr(y, "rs_meta") <- utils::modifyList(attr(x, "rs_meta"),
                                           list(source = "landsat-8"))
  expect_error(validate_data_fusion(y), "must contain `date`")
})

test_that("fusion validation succeeds in standalone, data, and raster modes", {
  rs <- .make_test_rs(field = FALSE)
  loo <- lapply(rs$scenes, function(sc){
    list(date = sc$date, predicted_rast = sc$index + 0.1,
         target_rast = sc$index)
  })
  testthat::local_mocked_bindings(
    .validateBuildRs = function(...) rs,
    .validateCompute = function(...) loo,
    .package = "ocular")

  standalone <- validate_data_fusion(
    117.8, latitude = -32, start_date = "2021-01-01",
    end_date = "2021-02-01", x_metres = 100)
  expect_equal(nrow(standalone), 2L)
  expect_equal(standalone$mae, c(0.1, 0.1), tolerance = 1e-7)

  series <- as_time_series(rs)
  series_out <- validate_data_fusion(series)
  expect_true(all(c("rmse", "mae", "bias", "r", "n_pixels") %in%
                    names(series_out)))
  expect_equal(series_out$mae, c(0.1, 0.1), tolerance = 1e-7)

  raster <- as_raster(rs)
  raster_out <- validate_data_fusion(raster)
  expect_s4_class(raster_out, "SpatRaster")
  expect_identical(names(raster_out),
                   c("rmse", "mae", "bias", "r", "n_scenes"))
})

test_that("unsupported fusion indices derive from the actual coarse lookup", {
  unsupported <- ocular:::.fusionUnsupportedIndices()
  expect_setequal(unsupported,
                  setdiff(names(ocular:::landsat_index_list),
                          names(ocular:::mcd43a4_index_list)))
})

test_that("per-date coarse retrieval failures are observable", {
  testthat::local_mocked_bindings(
    .fetchStac = function(...) stop("synthetic per-date failure"),
    .package = "ocular")
  expect_warning(
    out <- ocular:::.fetchCoarseNearDates(
      bbox = sf::st_bbox(c(xmin = 117, ymin = -33, xmax = 118, ymax = -32),
                         crs = 4326),
      fine_dates = as.Date(c("2021-01-01", "2021-01-09")),
      index_name = "EVI2", max_cloud_cover = 50),
    "2 fine-date window.*synthetic per-date failure")
  expect_length(out, 0L)
})

test_that("add_modis freezes fine-only detection scenes and is idempotent", {
  rs <- .make_test_rs()
  expect_error(add_modis(rs, typo = TRUE), "does not accept")
  rs$internals$detection_scenes <- NULL
  testthat::local_mocked_bindings(
    .fetchStac = function(...) NULL, .package = "ocular")
  out <- suppressMessages(add_modis(rs))
  expect_length(out$internals$detection_scenes, length(rs$scenes))

  fused <- rs$scenes[[1L]]
  fused$fused <- TRUE
  rs$scenes <- c(rs$scenes, list(fused))
  rs$internals$detection_scenes <- NULL
  out2 <- suppressMessages(add_modis(rs))
  expect_length(out2$scenes, length(rs$scenes))
  expect_true(all(!vapply(out2$internals$detection_scenes,
                          function(sc) isTRUE(sc$fused), logical(1L))))
})

test_that("add_modis reports coarse retrieval failures", {
  rs <- .make_test_rs()
  testthat::local_mocked_bindings(
    .fetchStac = function(...) stop("offline provider failure"),
    .package = "ocular")
  expect_warning(
    out <- suppressMessages(add_modis(rs)),
    "MCD43A4 retrieval failed: offline provider failure")
  expect_identical(out$scenes, rs$scenes)
})

test_that("add_modis appends a synthetic scene from an offline coarse provider", {
  rs <- .make_test_rs()
  coarse <- function(date, offset){
    r <- rs$scenes[[1L]]$index + offset
    list(date = as.Date(date),
         bands = list(Nadir_Reflectance_Band2 = r,
                      Nadir_Reflectance_Band1 = r / 2),
         index = r, item_assets = list())
  }
  coarse_scenes <- list(coarse("2021-01-05", 0),
                         coarse("2021-01-10", 0.01))
  testthat::local_mocked_bindings(
    .fetchStac = function(...) coarse_scenes,
    .fuseLandsatMODIS = function(fine_bands, ...){
      list(index = fine_bands[[1L]], bands = fine_bands,
           residual = fine_bands[[1L]] * 0,
           landsat_t1 = fine_bands[[1L]])
    },
    .package = "ocular")
  out <- suppressMessages(add_modis(rs))
  expect_equal(sum(vapply(out$scenes, function(sc) isTRUE(sc$fused),
                          logical(1L))), 1L)
  expect_length(out$internals$detection_scenes, 2L)
})

test_that("fusion failures are surfaced instead of silently discarded", {
  rs <- .make_test_rs()
  coarse <- function(date){
    r <- rs$scenes[[1L]]$index
    list(date = as.Date(date),
         bands = list(Nadir_Reflectance_Band2 = r,
                      Nadir_Reflectance_Band1 = r / 2),
         index = r, item_assets = list())
  }
  coarse_scenes <- list(coarse("2021-01-05"), coarse("2021-01-10"),
                         coarse("2021-01-20"))
  testthat::local_mocked_bindings(
    .fetchStac = function(...) coarse_scenes,
    .fuseLandsatMODIS = function(...) stop("synthetic fusion failure"),
    .package = "ocular")
  expect_warning(
    out <- suppressMessages(add_modis(rs)),
    "fusion failed.*synthetic fusion failure")
  expect_equal(length(out$scenes), length(rs$scenes))

  testthat::local_mocked_bindings(
    .fetchCoarseNearDates = function(...) coarse_scenes,
    .fuseLandsatMODIS = function(...) stop("synthetic validation failure"),
    .package = "ocular")
  expect_warning(
    loo <- ocular:::.validateCompute(rs, rs$params),
    "fusion failed.*synthetic validation failure")
  expect_null(loo)
})
