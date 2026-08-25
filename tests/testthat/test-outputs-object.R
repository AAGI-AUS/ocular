test_that("ocular objects validate and print without dumping raster contents", {
  rs <- .make_test_rs()
  expect_invisible(ocular:::.validateOcular(rs))
  printed <- capture.output(print(rs))
  expect_true(any(grepl("<ocular", printed, fixed = TRUE)))
  expect_true(any(grepl("sample-free field boundary delineation", printed,
                        fixed = TRUE)))
  expect_true(any(grepl("2 scenes", printed, fixed = TRUE)))
  expect_invisible(summary(rs))

  empty <- ocular:::.newOcular()
  empty$scenes <- list()
  expect_invisible(summary(empty))
  broken <- rs
  broken$state$alive_mat <- matrix(TRUE, 2, 2)
  expect_error(ocular:::.validateOcular(broken), "alive_mat")

  bad_date <- rs
  bad_date$scenes[[1L]]$date <- 1
  expect_error(ocular:::.validateOcular(bad_date), "invalid scene field")

  bad_params <- rs
  bad_params$params$not_an_ocular_parameter <- TRUE
  expect_error(ocular:::.validateOcular(bad_params), "unknown parameter")
})

test_that("new objects use the ocular class", {
  current <- ocular:::.newOcular()
  expect_identical(class(current), "ocular")
  expect_true(is_rs(current))
})

test_that("as_raster composites synthetic scenes and validates controls", {
  rs <- .make_test_rs(field = FALSE)
  out <- as_raster(rs, composite_function = "mean")
  expect_s4_class(out, "SpatRaster")
  expect_equal(terra::nlyr(out), 1L)
  expect_equal(as.numeric(terra::global(out, "mean", na.rm = TRUE)[[1L]]),
               0.11, tolerance = 1e-7)
  meta <- attr(out, "rs_meta")
  expect_identical(meta$params, rs$params)
  expect_equal(meta$max_cloud_cover, 50)
  expect_error(as_raster(rs, composite_function = "mode"),
               "should be one of")
  expect_error(as_raster(rs, mask = 1), "mask")

  no_crs <- sf::st_sfc(sf::st_point(c(0, 0)))
  expect_error(as_raster(rs, mask = no_crs), "coordinate reference")
})

test_that("an empty delineation never falls back to the whole AOI", {
  rs <- .make_test_rs(field = FALSE)
  rs$state$alive_mat <- matrix(FALSE, rs$geom$nr, rs$geom$nc)
  raster <- as_raster(rs)
  expect_true(all(is.na(terra::values(raster))))
  series <- as_time_series(rs)
  expect_equal(nrow(series), 0L)
  expect_false(is.null(attr(series, "rs_meta")))
})

test_that("all-missing sums remain missing rather than becoming zero", {
  rs <- .make_test_rs(field = FALSE)
  rs$scenes <- lapply(rs$scenes, function(sc){
    sc$index <- sc$index * NA_real_
    sc
  })
  raster <- as_raster(rs, composite_function = "sum")
  expect_true(all(is.na(terra::values(raster))))
  series <- as_time_series(rs, aggregate_function = "sum")
  expect_equal(nrow(series), 0L)
})

test_that("as_time_series supports explicit buckets and scalar aggregation", {
  rs <- .make_test_rs(field = FALSE)
  daily <- as_time_series(rs)
  expect_s3_class(daily, "data.frame")
  expect_equal(nrow(daily), 2L)
  expect_equal(daily$value, c(0.1, 0.12), tolerance = 1e-7)

  monthly <- as_time_series(rs, time_aggregate = "monthly",
                            aggregate_function = "mean")
  expect_equal(nrow(monthly), 1L)
  expect_equal(monthly$value, 0.11, tolerance = 1e-7)
  expect_error(as_time_series(rs, time_aggregate = "month"),
               "daily, monthly, or yearly")
  expect_error(as_time_series(rs, time_aggregate = "year"),
               "daily, monthly, or yearly")
  expect_error(as_time_series(rs, time_aggregate = "weekly"),
               "daily, monthly, or yearly")
  expect_error(as_time_series(rs, aggregate_function = "not_a_function"),
               "unknown aggregate_function")
  expect_error(as_time_series(rs, aggregate_function = quantile),
               "return one numeric")
  expect_error(as_time_series(rs, mask = list()), "mask")
})

test_that("fusion diagnostics use explicit anchor-mismatch names", {
  rs <- .make_test_rs(field = FALSE, fused = TRUE)
  out <- as_time_series(rs)
  expect_true(all(c("is_fused", "anchor_mae", "anchor_nmae") %in% names(out)))
  expect_false(any(c("MAPE", "MAE", "NMAE") %in% names(out)))
  expect_true(is.finite(out$anchor_nmae[out$is_fused][[1L]]))
})

test_that("band harmonisation checks full geometry, not x resolution alone", {
  a <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 20,
                   ymin = 0, ymax = 20, crs = "EPSG:32750")
  b <- terra::rast(nrows = 2, ncols = 2, xmin = 5, xmax = 25,
                   ymin = 0, ymax = 20, crs = "EPSG:32750")
  terra::values(a) <- 1:4
  terra::values(b) <- 5:8
  got <- ocular:::.harmoniseBands(list(a = a, b = b))
  expect_true(terra::compareGeom(got[[1L]], got[[2L]], stopOnError = FALSE))

  zero <- terra::setValues(a, c(0, 1, 1, 1))
  idx <- ocular:::.computeIndexFromBands(
    list(nir08 = zero, red = zero),
    list(assets = c("nir08", "red"),
         fun = function(r) (r$nir08 - r$red) / (r$nir08 + r$red)))
  expect_true(is.na(terra::values(idx)[[1L]]))
})
