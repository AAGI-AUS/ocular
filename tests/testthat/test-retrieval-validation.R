test_that("source resolution supports explicit sensors without fetching", {
  expect_identical(ocular:::.resolveSource("2021-01-01", "EVI2",
                                            source = "landsat-8")$source,
                   "landsat-8")
  expect_identical(ocular:::.resolveSource("2021-01-01", "EVI2",
                                            source = "sentinel-2")$source,
                   "sentinel-2")
  expect_identical(suppressMessages(
    ocular:::.resolveSource("2014-01-01", "EVI2", source = "auto")$source),
    "landsat-8")
  expect_error(ocular:::.resolveSource("2010-01-01", "EVI2",
                                        source = "sentinel-2"),
               "unavailable")
  expect_error(ocular:::.resolveSource("2010-01-01", "NDRE",
                                        source = "landsat-5"),
               "Unknown index")
})

test_that("get_rs rejects invalid inputs before any STAC request", {
  base <- list(longitude = 117.8, latitude = -32,
               start_date = "2021-01-01", end_date = "2021-02-01",
               x_metres = 100)
  call_get <- function(...) do.call(get_rs, utils::modifyList(base, list(...)))
  expect_error(call_get(longitude = Inf), "longitude")
  expect_error(call_get(latitude = 85), "UTM-supported")
  expect_error(call_get(source = "unknown"), "should be one of")
  expect_error(call_get(start_date = "bad"), "could not be parsed")
  expect_error(call_get(end_date = "2021-01-01"), "too short")
  expect_error(call_get(x_metres = 0), "positive")
  expect_error(call_get(max_cloud_cover = NA_real_), "max_cloud_cover")
  expect_error(call_get(scl_classes = c(1, 1)), "unique")
  expect_error(call_get(source = "sentinel-2", scl_classes = 12),
               "must be in")
  expect_error(call_get(source = "landsat-8", scl_classes = 16),
               "must be in")
  expect_error(call_get(index_name = "not-an-index"), "Unknown index_name")
  expect_error(call_get(params = list(search_start_date = "2020-12-01")),
               "must fall within")
})

test_that("calibration records the detection index, not the output index", {
  rs <- .make_test_rs()
  rs$spec$index_name <- "NDVI"
  rs$params$search_index <- "EVI2"
  cal <- ocular:::.calibrateScene(rs)
  expect_identical(cal$index_name, "EVI2")
  expect_equal(cal$c_noise, ocular:::.cNoise("EVI2"))
})

test_that("get_rs assembles a valid object from an offline scene provider", {
  fake_fetch <- function(bbox, start_date, end_date, index_name,
                         max_cloud_cover, scl_classes = NULL, source = NULL) {
    template <- terra::rast(nrows = 20, ncols = 20,
                            xmin = bbox[["xmin"]], xmax = bbox[["xmax"]],
                            ymin = bbox[["ymin"]], ymax = bbox[["ymax"]],
                            crs = "EPSG:4326")
    r1 <- terra::setValues(template, rep(0.5, terra::ncell(template)))
    r2 <- terra::setValues(template, rep(0.55, terra::ncell(template)))
    make_scene <- function(date, r) {
      list(date = as.Date(date), bands = list(nir08 = r, red = r),
           index = r, item_assets = list())
    }
    list(make_scene("2021-01-05", r1), make_scene("2021-01-20", r2))
  }
  testthat::local_mocked_bindings(.fetchStac = fake_fetch,
                                  .package = "ocular")
  explicit <- rs_params(area_sensitivity = 0.2)
  rs <- suppressMessages(get_rs(117.8, -32, "2021-01-01", "2021-02-01",
                                source = "landsat-8", x_metres = 300,
                                params = explicit, search_index = "EVI2"))
  expect_true(is_rs(rs))
  expect_identical(rs$geom$source, "landsat-8")
  expect_equal(length(rs$scenes), 2L)
  expect_s4_class(rs$internals$feature_stack, "SpatRaster")
  expect_setequal(attr(rs$params, "user_set", exact = TRUE),
                  c("area_sensitivity", "search_index"))
  expect_invisible(ocular:::.validateOcular(rs))

  partial <- suppressMessages(get_rs(
    117.8, -32, "2021-01-01", "2021-02-01",
    source = "landsat-8", x_metres = 300,
    params = list(strictness = "tight")))
  expect_identical(attr(partial$params, "user_set", exact = TRUE),
                   "strictness")
})

test_that("same-date adjacent fine tiles are mosaicked once", {
  make_tile <- function(xmin, xmax, value, href) {
    r <- terra::rast(nrows = 2, ncols = 2, xmin = xmin, xmax = xmax,
                     ymin = 0, ymax = 2, crs = "EPSG:4326")
    r <- terra::setValues(r, value)
    list(date = as.Date("2021-01-05"),
         bands = list(nir08 = r, red = r / 2), index = r,
         item_assets = list(nir08 = href, red = paste0(href, "-red")),
         pixel_mask = NULL)
  }
  scenes <- list(make_tile(0, 2, 1, "tile-a"),
                 make_tile(2, 4, 2, "tile-b"))
  got <- ocular:::.mosaicFineScenesByDate(
    scenes, ocular:::landsat_index_list[["NDVI"]])
  expect_length(got, 1L)
  expect_equal(c(terra::xmin(got[[1L]]$index),
                 terra::xmax(got[[1L]]$index)), c(0, 4))
  expect_equal(got[[1L]]$item_assets$nir08, c("tile-a", "tile-b"))
})

test_that("alternate detection assets inherit the stored pixel mask", {
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2,
                   ymin = 0, ymax = 2, crs = "EPSG:4326")
  r <- terra::setValues(r, 1:4)
  pm <- terra::setValues(r, c(TRUE, FALSE, FALSE, FALSE))
  scene <- list(date = as.Date("2021-01-01"), bands = list(), index = r,
                item_assets = list(extra = "unused"), pixel_mask = pm)
  testthat::local_mocked_bindings(
    .fetchAssets = function(...) list(extra = r), .package = "ocular")
  got <- ocular:::.composeSceneFromBands(
    scene, list(assets = "extra", fun = function(x) x[[1L]]),
    sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2, ymax = 2), crs = 4326),
    "sentinel-2")
  expect_true(is.na(terra::values(got$index, mat = FALSE)[[1L]]))
  expect_equal(sum(is.finite(terra::values(got$index, mat = FALSE))), 3L)
})

test_that("feature composites preserve all-missing cells as NA", {
  rs <- .make_test_rs()
  rs$internals$feature_stack <- NULL
  rs$internals$calibration <- NULL
  rs$scenes <- lapply(rs$scenes, function(sc){
    values <- terra::values(sc$index, mat = FALSE)
    values[[1L]] <- NA_real_
    sc$index <- terra::setValues(sc$index, values)
    sc
  })
  out <- suppressMessages(ocular:::.fetchFeatureStack(rs))
  values <- terra::values(out$internals$feature_stack, mat = FALSE)
  expect_false(any(is.infinite(values)))
})

test_that("sinusoidal search padding remains within world bounds", {
  seen <- NULL
  testthat::local_mocked_bindings(
    .fetchStacItems = function(..., bbox_vec){
      seen <<- bbox_vec
      list(features = list())
    }, .package = "ocular")
  bbox <- sf::st_bbox(c(xmin = 179.97, ymin = 89.95,
                        xmax = 179.99, ymax = 89.99), crs = 4326)
  expect_null(ocular:::.fetchStac(
    bbox, "2021-01-01", "2021-01-03", "NDVI", 100,
    source = "mcd43a4"))
  expect_gte(seen[[1L]], -180)
  expect_gte(seen[[2L]], -90)
  expect_lte(seen[[3L]], 180)
  expect_lte(seen[[4L]], 90)
})
