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

test_that("get_rs names search_index when the source lacks it", {
  testthat::local_mocked_bindings(
    .fetchStac = function(...) stop("STAC fetch reached"),
    .package = "ocular")
  expect_error(
    get_rs(117.8, -32, "2021-01-01", "2021-03-01", source = "landsat-8",
           x_metres = 300, search_index = "NDRE"),
    "search_index 'NDRE' is not known for source 'landsat-8'.",
    fixed = TRUE)
})

test_that("get_rs names the source for an unsupported output index", {
  n_fetch <- 0L
  testthat::local_mocked_bindings(
    .fetchStac = function(...){
      n_fetch <<- n_fetch + 1L
      stop("STAC fetch reached")
    }, .package = "ocular")
  cases <- data.frame(
    source = c("landsat-8", "landsat-5", "sentinel-2", "auto", "auto", "auto"),
    date = c("2021-01-01", "2010-01-01", "2021-01-01",
             "2021-01-01", "2014-01-01", "2010-01-01"),
    resolved = c("landsat-8", "landsat-5", "sentinel-2",
                 "sentinel-2", "landsat-8", "landsat-5"))
  for( i in seq_len(nrow(cases)) ){
    expect_error(
      suppressMessages(get_rs(
        117.8, -32, cases$date[i], "2021-03-01", source = cases$source[i],
        x_metres = 300, index_name = "NOPE")),
      paste0("Unknown index_name 'NOPE' for source '", cases$resolved[i], "'."),
      fixed = TRUE)
  }
  expect_error(
    get_rs(117.8, -32, "2021-01-01", "2021-03-01", source = "landsat-8",
           x_metres = 300, index_name = "NDRE", search_index = "NDRE"),
    "Unknown index_name 'NDRE' for source 'landsat-8'.", fixed = TRUE)
  expect_identical(n_fetch, 0L)
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
                         max_cloud_cover, scl_classes = NULL, source = NULL,
                         s2_processing = "highest_baseline") {
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

test_that("standalone output calls skip the feature stack however written", {
  ## The middle composite window has no scene, so building the feature
  ## stack would stop.
  fake_fetch <- function(bbox, start_date, end_date, index_name,
                         max_cloud_cover, scl_classes = NULL, source = NULL,
                         s2_processing = "highest_baseline") {
    template <- terra::rast(nrows = 10, ncols = 10,
                            xmin = bbox[["xmin"]], xmax = bbox[["xmax"]],
                            ymin = bbox[["ymin"]], ymax = bbox[["ymax"]],
                            crs = "EPSG:4326")
    r <- terra::setValues(template, rep(0.5, terra::ncell(template)))
    lapply(c("2021-01-05", "2021-01-21", "2021-05-14"), function(date) {
      list(date = as.Date(date), bands = list(nir08 = r, red = r),
           index = r, item_assets = list())
    })
  }
  testthat::local_mocked_bindings(.fetchStac = fake_fetch,
                                  .package = "ocular")
  args <- list(117.8, latitude = -32, start_date = "2021-01-01",
               end_date = "2021-06-01", source = "landsat-8",
               x_metres = 300)
  bare <- suppressMessages(do.call("as_time_series", args))
  expect_equal(nrow(bare), 3L)
  expect_identical(
    suppressMessages(ocular::as_time_series(
      117.8, latitude = -32, start_date = "2021-01-01",
      end_date = "2021-06-01", source = "landsat-8", x_metres = 300)),
    bare)
  expect_identical(suppressMessages(do.call(as_time_series, args)), bare)
  expect_identical(
    suppressMessages(lapply(117.8, as_time_series, latitude = -32,
                            start_date = "2021-01-01",
                            end_date = "2021-06-01",
                            source = "landsat-8", x_metres = 300))[[1L]],
    bare)
  expect_s4_class(suppressMessages(do.call(as_raster, args)), "SpatRaster")
  expect_error(
    suppressMessages(get_rs(117.8, -32, "2021-01-01", "2021-06-01",
                            source = "landsat-8", x_metres = 300)),
    "no scenes in window")
})

test_that("standalone output calls read no detection assets", {
  ## NBR scenes lack the red band that the EVI2 detection index needs.
  fake_fetch <- function(bbox, start_date, end_date, index_name,
                         max_cloud_cover, scl_classes = NULL, source = NULL,
                         s2_processing = "highest_baseline") {
    template <- terra::rast(nrows = 10, ncols = 10,
                            xmin = bbox[["xmin"]], xmax = bbox[["xmax"]],
                            ymin = bbox[["ymin"]], ymax = bbox[["ymax"]],
                            crs = "EPSG:4326")
    r <- terra::setValues(template, rep(0.5, terra::ncell(template)))
    lapply(c("2021-01-05", "2021-01-20"), function(date) {
      list(date = as.Date(date), bands = list(nir08 = r, swir22 = r),
           index = r, item_assets = list(red = "https://example.invalid/red.tif"))
    })
  }
  n_reads <- 0L
  failing_read <- function(...) {
    n_reads <<- n_reads + 1L
    stop("asset unreadable")
  }
  testthat::local_mocked_bindings(.fetchStac = fake_fetch,
                                  .fetchAssets = failing_read,
                                  .package = "ocular")
  args <- list(117.8, latitude = -32, start_date = "2021-01-01",
               end_date = "2021-02-01", source = "landsat-8",
               x_metres = 300, index_name = "NBR")
  ts <- suppressMessages(do.call(as_time_series, args))
  expect_equal(nrow(ts), 2L)
  expect_s4_class(suppressMessages(do.call(as_raster, args)), "SpatRaster")
  expect_identical(n_reads, 0L)
  expect_error(
    suppressWarnings(suppressMessages(do.call(get_rs, args))),
    "No usable detection scenes remain")
  expect_identical(n_reads, 2L)
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

test_that("feature composites keep detection scenes in another CRS or extent", {
  lon <- 119.995
  lat <- -32
  bbox <- point_to_bbox(lon, lat, 600, 600)
  aoi <- sf::st_as_sfc(bbox)
  scene_in <- function(date, crs, value, east_cut = 2000){
    xy <- terra::geom(terra::project(
      terra::vect(cbind(lon, lat), crs = "EPSG:4326"), crs))[1L, c("x", "y")]
    x0 <- floor((xy[["x"]] - 2000) / 10) * 10
    y0 <- floor((xy[["y"]] - 2000) / 10) * 10
    tile <- terra::rast(xmin = x0, xmax = floor((xy[["x"]] + east_cut) / 10) * 10,
                        ymin = y0, ymax = y0 + 4000, resolution = 10, crs = crs)
    terra::values(tile) <- value
    e <- terra::ext(sf::st_bbox(sf::st_transform(aoi, terra::crs(tile))))
    list(date = as.Date(date), bands = list(),
         index = terra::crop(tile, e, snap = "out"),
         item_assets = list(), pixel_mask = NULL)
  }
  stack_of <- function(scenes){
    rs <- ocular:::.newOcular()
    rs$spec$point$longitude <- lon
    rs$spec$point$latitude <- lat
    rs$spec$start_date <- "2021-08-01"
    rs$spec$end_date <- "2021-08-15"
    rs$spec$x_metres <- 600
    rs$spec$y_metres <- 600
    rs$spec$index_name <- "NDVI"
    rs$geom$bbox <- bbox
    rs$geom$source <- "sentinel-2"
    rs$params <- rs_params(search_index = "EVI2")
    rs$internals$detection_scenes <- scenes
    out <- suppressWarnings(ocular:::.fetchFeatureStack(rs))
    return(out$internals$feature_stack)
  }

  ## The window maximum includes the later date in the neighbouring UTM zone.
  fs <- stack_of(list(scene_in("2021-08-03", "EPSG:32750", 0.4),
                      scene_in("2021-08-08", "EPSG:32751", 0.6)))
  expect_equal(max(terra::values(fs), na.rm = TRUE), 0.6, tolerance = 1e-6)
  expect_false(anyNA(terra::values(fs)))

  ## A first scene cut by its tile edge does not crop the later full scene.
  full <- stack_of(list(scene_in("2021-08-03", "EPSG:32750", 0.4),
                        scene_in("2021-08-08", "EPSG:32750", 0.6)))
  cut <- stack_of(list(scene_in("2021-08-03", "EPSG:32750", 0.4, east_cut = 100),
                       scene_in("2021-08-08", "EPSG:32750", 0.6)))
  expect_identical(dim(cut), dim(full))
  expect_equal(as.vector(terra::ext(cut)), as.vector(terra::ext(full)))
})

test_that("asset read errors are reported once without the URL query", {
  original_rast <- terra::rast
  dir <- .test_artifact_dir("read-errors")
  bbox <- point_to_bbox(117.8, -32, 300, 300)
  e <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(bbox), "EPSG:32750"))
  fixture <- original_rast(xmin = e[["xmin"]] - 100, xmax = e[["xmax"]] + 100,
                           ymin = e[["ymin"]] - 100, ymax = e[["ymax"]] + 100,
                           res = 10, crs = "EPSG:32750")
  fixture <- terra::setValues(fixture, rep(3000L, terra::ncell(fixture)))
  path <- file.path(dir, "asset.tif")
  terra::writeRaster(fixture, path, datatype = "INT2U")
  href <- function(key) {
    paste0("https://fixture.invalid/", key, ".tif?sig=SECRET")
  }
  testthat::local_mocked_bindings(rast = function(x, ...) {
    if (missing(x)) return(original_rast(...))
    if (is.character(x) && length(x) == 1L &&
        startsWith(x, "/vsicurl/https://fixture.invalid/")) {
      if (grepl("bad", x, fixed = TRUE)) {
        stop("HTTP response code: 403 for ", x)
      }
      return(original_rast(path))
    }
    original_rast(x, ...)
  }, .package = "terra")
  items <- NULL
  testthat::local_mocked_bindings(
    .fetchStacItems = function(...) list(features = items),
    .package = "ocular")

  fine_item <- function(date, key) {
    list(properties = list(datetime = paste0(date, "T02:00:00Z"),
                           `s2:processing_baseline` = "03.01", platform = "Sentinel-2A",
                           `s2:mgrs_tile` = "55HED", `sat:relative_orbit` = 73L,
                           `s2:datastrip_id` = "S2A_OPER_MSI_L2A_DS_TEST_S20210601T020000_N03.01"),
         assets = list(B08 = list(href = href(key)),
                       B04 = list(href = href(paste0(key, "-red")))))
  }
  items <- list(fine_item("2021-06-02", "a"), fine_item("2021-06-17", "bad"),
                fine_item("2021-07-02", "c"))
  expect_warning(
    scenes <- ocular:::.fetchStac(bbox, "2021-06-01", "2021-07-31", "NDVI",
                                  50, source = "sentinel-2"),
    "^1 of 3 acquisitions dropped after an error. First error: .*bad\\.tif$")
  expect_length(scenes, 2L)

  tile_item <- function(key) {
    list(properties = list(start_datetime = "2021-06-02T00:00:00Z"),
         assets = list(Nadir_Reflectance_Band2 = list(href = href(key)),
                       Nadir_Reflectance_Band1 = list(href = href(key))))
  }
  items <- list(tile_item("h1"), tile_item("bad"))
  expect_warning(
    scenes <- ocular:::.fetchStac(bbox, "2021-06-01", "2021-06-05", "NDVI",
                                  100, source = "mcd43a4"),
    "^MCD43A4 processing errors: 1\\. First error: .*bad\\.tif$")
  expect_length(scenes, 1L)

  testthat::local_mocked_bindings(
    .finiteIndex = function(x) stop("index step failed"),
    .package = "ocular")
  items <- list(tile_item("h1"))
  expect_warning(
    scenes <- ocular:::.fetchStac(bbox, "2021-06-01", "2021-06-05", "NDVI",
                                  100, source = "mcd43a4"),
    "^MCD43A4 processing errors: 1\\. First error: index step failed$")
  expect_null(scenes)
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

test_that("a failing STAC query stops with its error after the last attempt", {
  attempts <- 0L
  failing <- TRUE
  testthat::local_mocked_bindings(
    post_request = function(q, ...){
      attempts <<- attempts + 1L
      if (failing) stop("HTTP status 400: invalid filter")
      list(type = "FeatureCollection", features = list())
    },
    items_fetch = function(items, ...) items,
    items_sign = function(items, ...) items,
    .package = "rstac")
  bbox <- point_to_bbox(117.8, -32, 300, 300)
  bbox_vec <- unname(unclass(bbox)[c("xmin", "ymin", "xmax", "ymax")])
  expect_error(
    ocular:::.fetchStacItems("https://stac.invalid/api", "landsat-c2-l2",
                             bbox_vec, "2021-01-01", "2021-03-01", 50, TRUE,
                             initial_pause = 0, max_pause = 0),
    "^STAC query failed after 10 attempts: HTTP status 400: invalid filter$")
  expect_identical(attempts, 10L)

  ## A query that succeeds with no features is not a failure.
  failing <- FALSE
  attempts <- 0L
  expect_null(ocular:::.fetchStac(bbox, "2021-01-01", "2021-03-01", "NDVI",
                                  50, source = "landsat-8"))
  expect_identical(attempts, 1L)
})

test_that("the STAC query includes every instant of end_date", {
  sent <- NULL
  testthat::local_mocked_bindings(
    cql2_interval = function(start, end){
      sent <<- c(start, end)
      stop("interval recorded")
    },
    .package = "rstac")
  bbox <- point_to_bbox(117.8, -32, 300, 300)
  for (src in c("sentinel-2", "mcd43a4")) {
    sent <- NULL
    expect_error(ocular:::.fetchStac(bbox, "2021-06-01", "2021-06-30", "NDVI",
                                     50, source = src),
                 "interval recorded")
    expect_identical(sent, c("2021-06-01", "2021-06-30T23:59:59.999999Z"))
  }
})

test_that("MCD43A4 scenes take the product date and stay within the period", {
  original_rast <- terra::rast
  dir <- .test_artifact_dir("mcd43a4-dates")
  bbox <- point_to_bbox(117.8, -32, 300, 300)
  e <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(bbox), "EPSG:32750"))
  fixture <- original_rast(xmin = e[["xmin"]] - 100, xmax = e[["xmax"]] + 100,
                           ymin = e[["ymin"]] - 100, ymax = e[["ymax"]] + 100,
                           res = 10, crs = "EPSG:32750")
  fixture <- terra::setValues(fixture, rep(3000L, terra::ncell(fixture)))
  path <- file.path(dir, "asset.tif")
  terra::writeRaster(fixture, path, datatype = "INT2U", overwrite = TRUE)
  testthat::local_mocked_bindings(rast = function(x, ...) {
    if (missing(x)) return(original_rast(...))
    if (is.character(x) && startsWith(x, "/vsicurl/https://fixture.invalid/")) {
      return(original_rast(path))
    }
    original_rast(x, ...)
  }, .package = "terra")
  ## As the service does, return every product whose 16-day window
  ## (start_datetime, 8 days before datetime) overlaps 2021-06-10 to 06-12.
  product_dates <- seq(as.Date("2021-06-03"), as.Date("2021-06-20"), by = "day")
  items <- lapply(product_dates, function(d) {
    href <- paste0("https://fixture.invalid/", format(d), ".tif")
    list(properties = list(datetime = paste0(format(d), "T00:00:00Z"),
                           start_datetime = paste0(format(d - 8), "T00:00:00Z")),
         assets = list(Nadir_Reflectance_Band2 = list(href = href),
                       Nadir_Reflectance_Band1 = list(href = href)))
  })
  testthat::local_mocked_bindings(
    .fetchStacItems = function(...) list(features = items),
    .package = "ocular")
  scenes <- ocular:::.fetchStac(bbox, "2021-06-10", "2021-06-12", "NDVI", 100,
                                source = "mcd43a4")
  expect_identical(vapply(scenes, function(s) format(s$date), character(1L)),
                   c("2021-06-10", "2021-06-11", "2021-06-12"))
})

test_that("the Sentinel-2 offset is removed for each item from baseline 04.00", {
  ## Every item shows one surface: reflectance B04 0.0563, B08 0.3, B05 0.12.
  ## From baseline 04.00 an item stores DN = 10000 x reflectance + 1000.
  original_rast <- terra::rast
  dir <- .test_artifact_dir("s2-offset")
  bbox <- point_to_bbox(117.8, -32, 300, 300)
  e <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(bbox), "EPSG:32750"))
  x0 <- floor(e[["xmin"]] / 10) * 10 - 100
  x1 <- ceiling(e[["xmax"]] / 10) * 10 + 100
  y0 <- floor(e[["ymin"]] / 10) * 10 - 100
  y1 <- ceiling(e[["ymax"]] / 10) * 10 + 100
  xm <- round((x0 + x1) / 20) * 10
  refl <- c(B04 = 563, B08 = 3000, B05 = 1200)
  tiles <- list(west = c(x0, xm), east = c(xm, x1), full = c(x0, x1))
  files <- list()
  for (tn in names(tiles)) {
    r <- original_rast(xmin = tiles[[tn]][1], xmax = tiles[[tn]][2],
                       ymin = y0, ymax = y1, res = 10, crs = "EPSG:32750")
    for (off in c(0, 1000)) {
      for (b in names(refl)) {
        key <- sprintf("%s-%s-%d", tn, b, off)
        files[[key]] <- file.path(dir, paste0(key, ".tif"))
        terra::writeRaster(terra::setValues(r, refl[[b]] + off), files[[key]],
                           datatype = "INT2U", overwrite = TRUE)
      }
    }
    ## SCL class 8 at one cell of each tile; the mask is never shifted.
    scl <- rep(4L, terra::ncell(r))
    centre <- cbind(mean(tiles[[tn]]), (y0 + y1) / 2)
    scl[[terra::cellFromXY(r, centre)]] <- 8L
    files[[paste0(tn, "-SCL")]] <- file.path(dir, paste0(tn, "-SCL.tif"))
    terra::writeRaster(terra::setValues(r, scl), files[[paste0(tn, "-SCL")]],
                       datatype = "INT1U", overwrite = TRUE)
  }
  testthat::local_mocked_bindings(rast = function(x, ...) {
    if (missing(x)) return(original_rast(...))
    if (is.character(x) && startsWith(x, "/vsicurl/https://fixture.invalid/")) {
      return(original_rast(files[[sub("\\.tif$", "", basename(x))]]))
    }
    original_rast(x, ...)
  }, .package = "terra")
  item <- function(date, tile, baseline, off) {
    href <- function(b) {
      key <- if (b == "SCL") paste0(tile, "-SCL") else sprintf("%s-%s-%d", tile, b, off)
      paste0("https://fixture.invalid/", date, "/", key, ".tif")
    }
    props <- list(datetime = paste0(date, "T02:00:00.024Z"), platform = "Sentinel-2A",
                   `s2:mgrs_tile` = tile, `sat:relative_orbit` = 73L,
                   `s2:datastrip_id` = "S2A_OPER_MSI_L2A_DS_TEST_S20220601T020000_N04.00")
    props[["s2:processing_baseline"]] <- baseline
    list(properties = props,
         assets = lapply(c(B04 = "B04", B08 = "B08", B05 = "B05", SCL = "SCL"),
                         function(b) list(href = href(b))))
  }
  items <- list(item("2022-06-03", "west", "04.00", 1000),
                item("2022-06-03", "east", "03.01", 0),
                item("2022-06-05", "full", NULL, 1000),
                item("2022-06-08", "full", "05.10", 1000),
                item("2022-06-10", "full", "03.01", 0))
  testthat::local_mocked_bindings(
    .fetchStacItems = function(...) list(features = items),
    .package = "ocular")
  ## search_index NDRE needs B05, which is read later through .fetchAssets().
  expect_warning(
    rs <- get_rs(117.8, -32, "2022-06-01", "2022-06-11", index_name = "NDVI",
                 source = "sentinel-2", x_metres = 300, scl_classes = 8,
                 search_index = "NDRE"),
    paste("^No offset was removed from 1 Sentinel-2 item\\(s\\), because",
          "the processing baseline could not be read\\.$"))
  vals <- function(r) {
    v <- terra::values(r, mat = FALSE)
    round(unique(v[is.finite(v)]), 6)
  }
  dates <- vapply(rs$scenes, function(sc) format(sc$date), character(1L),
                  USE.NAMES = FALSE)
  expect_identical(dates, c("2022-06-03", "2022-06-05", "2022-06-08", "2022-06-10"))
  b04 <- unname(lapply(rs$scenes, function(sc) vals(sc$bands$B04)))
  expect_equal(b04, list(0.0563, 0.1563, 0.0563, 0.0563))
  b05 <- unname(lapply(rs$internals$detection_scenes,
                      function(sc) vals(sc$bands$B05)))
  expect_equal(b05, list(0.12, 0.22, 0.12, 0.12))

  ## The same-date pair keeps one offset for each href, in href order.
  pair <- rs$scenes[[1L]]
  expect_equal(pair$dn_offset$B05, c(1000, 0))
  expect_match(pair$item_assets$B05[[1L]], "west", fixed = TRUE)
  expect_equal(sum(is.na(terra::values(pair$bands$B04, mat = FALSE))), 2L)
  expect_equal(rs$scenes[[3L]]$dn_offset$B04, 1000)
  ## Nothing new is stored on a scene without offset.
  expect_identical(names(rs$scenes[[4L]]),
                   c("date", "bands", "index", "item_assets", "pixel_mask"))
  expect_null(rs$scenes[[2L]]$dn_offset)
})

test_that("coarse merge fallback redacts URL queries and preserves the first tile", {
  original_rast <- terra::rast
  template <- original_rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2,
                            ymin = 0, ymax = 2, crs = "EPSG:4326")
  bands <- list(nir = terra::setValues(template, c(3000, 4000, 5000, 6000)),
                red = terra::setValues(template, c(1000, 1000, 2000, 2000)))
  item <- list(properties = list(datetime = "2021-06-02T00:00:00Z"),
                assets = list(Nadir_Reflectance_Band2 = list(href = "fixture-nir"),
                              Nadir_Reflectance_Band1 = list(href = "fixture-red")))
  items <- list(item)
  testthat::local_mocked_bindings(rast = function(x, ...){
    if( missing(x) ){
      return(original_rast(...))
    }
    if( is.character(x) && length(x) == 1L && startsWith(x, "/vsicurl/fixture-") ){
      return(bands[[sub("/vsicurl/fixture-", "", x, fixed = TRUE)]])
    }
    return(original_rast(x, ...))
  }, .package = "terra")
  testthat::local_mocked_bindings(
    .fetchStacItems = function(...){ return(list(features = items)) },
    .mergeRasterTiles = function(...){
      stop(paste0("merge read failed for https://fixture.invalid/tile.tif",
                  "?", "sig", "=", "synthetic-credential"))
    }, .package = "ocular")
  run <- function(){
    return(ocular:::.fetchStac(
      sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2, ymax = 2), crs = 4326),
      "2021-06-01", "2021-06-05", "NDVI", 100, source = "mcd43a4"))
  }
  control <- run()
  items <- list(item, item)
  notices <- character(0L)
  fallback <- withCallingHandlers(run(), warning = function(w){
    notices <<- c(notices, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(notices, 2L)
  expect_true(all(grepl("Using the first tile only", notices, fixed = TRUE)))
  expect_false(any(grepl("sig=|synthetic-credential", notices)))
  expect_true(all(grepl("https://fixture.invalid/tile.tif$", notices)))
  expect_length(fallback, 1L)
  expect_identical(fallback[[1L]]$date, control[[1L]]$date)
  expect_identical(fallback[[1L]]$item_assets, control[[1L]]$item_assets)
  expect_identical(names(fallback[[1L]]$bands), names(control[[1L]]$bands))
  for( nm in names(control[[1L]]$bands) ){
    expect_true(terra::compareGeom(fallback[[1L]]$bands[[nm]], control[[1L]]$bands[[nm]]))
    expect_identical(terra::values(fallback[[1L]]$bands[[nm]]),
                     terra::values(control[[1L]]$bands[[nm]]))
  }
  expect_true(terra::compareGeom(fallback[[1L]]$index, control[[1L]]$index))
  expect_identical(terra::values(fallback[[1L]]$index), terra::values(control[[1L]]$index))
})
