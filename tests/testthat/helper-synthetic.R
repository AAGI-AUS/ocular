.test_artifact_dir <- function(label = "case") {
  root <- Sys.getenv("OCULAR_TEST_ARTIFACTS", unset = "")
  if (!nzchar(root)) root <- file.path(getwd(), ".test-artifacts")
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  path <- tempfile(pattern = paste0(label, "-"), tmpdir = root)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

.make_test_rs <- function(n = 21L, field = TRUE, fused = FALSE) {
  longitude <- 117.8
  latitude <- -32
  utm <- "EPSG:32750"
  point <- terra::vect(matrix(c(longitude, latitude), nrow = 1L),
                       crs = "EPSG:4326")
  xy <- terra::geom(terra::project(point, utm))[1L, c("x", "y")]
  half <- n * 10 / 2
  template <- terra::rast(nrows = n, ncols = n,
                          xmin = xy[["x"]] - half,
                          xmax = xy[["x"]] + half,
                          ymin = xy[["y"]] - half,
                          ymax = xy[["y"]] + half,
                          crs = utm)

  values <- matrix(0.1, n, n)
  if (isTRUE(field)) values[5:(n - 4L), 5:(n - 4L)] <- 0.6
  w1 <- terra::setValues(template, as.vector(t(values)))
  w2 <- terra::setValues(template, as.vector(t(values + 0.02)))
  feature_stack <- c(w1, w2)
  names(feature_stack) <- c("w1", "w2")

  poly_wgs <- terra::project(terra::as.polygons(terra::ext(template), crs = utm),
                             "EPSG:4326")
  bbox <- sf::st_bbox(sf::st_as_sf(poly_wgs))
  scene <- function(date, raster) {
    list(date = as.Date(date), bands = list(nir08 = raster, red = raster),
         index = raster, item_assets = list())
  }

  rs <- ocular:::.newOcular()
  rs$spec$point <- list(longitude = longitude, latitude = latitude)
  rs$spec$start_date <- "2021-01-01"
  rs$spec$end_date <- "2021-02-01"
  rs$spec$x_metres <- n * 10
  rs$spec$y_metres <- n * 10
  rs$spec$index_name <- "EVI2"
  rs$spec$max_cloud_cover <- 50
  rs$geom$bbox <- bbox
  rs$geom$source <- "landsat-8"
  rs$geom$pixel_size_m <- 10
  rs$geom$nr <- n
  rs$geom$nc <- n
  rs$geom$centre_rc <- c((n + 1L) %/% 2L, (n + 1L) %/% 2L)
  rs$scenes <- list(scene("2021-01-05", w1), scene("2021-01-20", w2))
  rs$internals$feature_stack <- feature_stack
  rs$state$decided_mat <- matrix(FALSE, n, n)
  rs$params <- rs_params()
  rs$internals$calibration <- ocular:::.calibrateScene(rs)

  if (isTRUE(fused)) {
    rs$scenes[[1L]]$fused <- FALSE
    rs$scenes[[1L]]$anchor_days <- NA_integer_
    rs$scenes[[1L]]$residual <- NULL
    rs$scenes[[1L]]$landsat_t1 <- NULL
    rs$scenes[[2L]]$fused <- TRUE
    rs$scenes[[2L]]$anchor_days <- 2L
    rs$scenes[[2L]]$residual <- abs(w2 - w1)
    rs$scenes[[2L]]$landsat_t1 <- w1
  }
  rs
}
