test_that("FTW calibration branch and diagnose_against_ftw work on a synthetic field", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")

  # Synthetic 30x30 feature grid, 2 windows, metric CRS (UTM 50S).
  r <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 300,
                   ymin = 0, ymax = 300, crs = "EPSG:32750")
  field <- matrix(0.20, 30, 30); field[8:23, 8:23] <- 0.80
  w1 <- terra::setValues(r, as.vector(t(field)))
  w2 <- terra::setValues(r, as.vector(t(field + 0.02)))
  fstack <- c(w1, w2); names(fstack) <- c("w1", "w2")

  poly_sfc <- sf::st_sfc(sf::st_polygon(list(rbind(
    c(80, 80), c(220, 80), c(220, 220), c(80, 220), c(80, 80)))), crs = 32750)
  poly_sf  <- sf::st_sf(geometry = poly_sfc)

  rs <- structure(
    list(geom      = list(nr = 30, nc = 30, pixel_size_m = 10,
                          ftw_prior = list(polygon = poly_sf)),
         internals = list(feature_stack = fstack),
         state     = list(alive_mat = NULL),
         spec      = list(index_name = "NDVI"),
         params    = list()),
    class = "ocular")

  # FTW-derived support measures the in-field signature on the polygon interior.
  ls <- ocular:::.ftwSupport(rs, c_noise = 0.02)
  expect_false(is.null(ls))
  expect_length(ls$sig$mu_w, 2L)
  expect_true(all(ls$sig$mu_w > 0.5))

  # Full calibration dispatches to the FTW branch and tags the source.
  cal <- ocular:::.calibrateScene(rs)
  expect_identical(cal$source, "ftw")

  # A delineation matching the field block agrees strongly with the reference.
  alive <- matrix(FALSE, 30, 30); alive[8:23, 8:23] <- TRUE
  rs$state$alive_mat <- alive
  d <- diagnose_against_ftw(rs)
  expect_true(d$iou > 0.7 && d$f1 > 0.7)

  # A half-covering delineation agrees less.
  alive2 <- matrix(FALSE, 30, 30); alive2[8:23, 8:15] <- TRUE
  rs$state$alive_mat <- alive2
  d2 <- diagnose_against_ftw(rs)
  expect_true(d2$iou < d$iou)

  # Non-empty but disjoint masks have a defined zero F1 score.
  disjoint <- matrix(FALSE, 30, 30); disjoint[1:2, 1:2] <- TRUE
  rs$state$alive_mat <- disjoint
  d0 <- diagnose_against_ftw(rs)
  expect_equal(d0$iou, 0)
  expect_equal(d0$f1, 0)
})

test_that("FTW selection keeps only the field containing the seed", {
  p1 <- sf::st_polygon(list(rbind(c(117.79, -32.01), c(117.80, -32.01),
                                  c(117.80, -31.99), c(117.79, -31.99),
                                  c(117.79, -32.01))))
  p2 <- sf::st_polygon(list(rbind(c(117.80, -32.01), c(117.82, -32.01),
                                  c(117.82, -31.99), c(117.80, -31.99),
                                  c(117.80, -32.01))))
  fields <- sf::st_sf(source_year = c(2019L, 2021L),
                      geometry = sf::st_sfc(p1, p2, crs = 4326))
  rs <- .make_test_rs()
  rs$spec$point <- list(longitude = 117.795, latitude = -32)
  got1 <- ocular:::.selectFtwField(fields, rs)
  expect_equal(got1$source_year, 2019L)

  rs$spec$point <- list(longitude = 117.81, latitude = -32)
  got2 <- ocular:::.selectFtwField(fields, rs)
  expect_equal(got2$source_year, 2021L)

  rs$spec$point <- list(longitude = 117.85, latitude = -32)
  expect_null(ocular:::.selectFtwField(fields, rs))
})

test_that("add_ftw_prior selects the seed field and recalibrates", {
  rs <- .make_test_rs()
  expect_identical(rs$internals$calibration$source, "self")
  bb <- rs$geom$bbox
  midx <- mean(c(bb[["xmin"]], bb[["xmax"]]))
  midy <- mean(c(bb[["ymin"]], bb[["ymax"]]))
  dx <- diff(c(bb[["xmin"]], bb[["xmax"]]))
  dy <- diff(c(bb[["ymin"]], bb[["ymax"]]))
  seed_field <- sf::st_polygon(list(rbind(
    c(midx - dx / 4, midy - dy / 4), c(midx + dx / 4, midy - dy / 4),
    c(midx + dx / 4, midy + dy / 4), c(midx - dx / 4, midy + dy / 4),
    c(midx - dx / 4, midy - dy / 4))))
  neighbour <- sf::st_polygon(list(rbind(
    c(midx + dx / 3, midy - dy / 4), c(midx + dx / 2, midy - dy / 4),
    c(midx + dx / 2, midy + dy / 4), c(midx + dx / 3, midy + dy / 4),
    c(midx + dx / 3, midy - dy / 4))))
  fields <- sf::st_sf(source_year = c(2020L, 2010L),
                      geometry = sf::st_sfc(seed_field, neighbour, crs = 4326))

  testthat::local_mocked_bindings(
    .ftwCacheRead = function(...) NULL,
    .ftwQuery = function(...) fields,
    .ftwCacheWrite = function(...) invisible(NULL),
    .ftwCacheEvictLRU = function(...) invisible(0L),
    .package = "ocular")
  out <- suppressMessages(add_ftw_prior(rs, source = "local.parquet"))
  expect_equal(nrow(out$geom$ftw_prior$polygon), 1L)
  expect_equal(out$geom$ftw_prior$source_year, 2020L)
  expect_identical(out$internals$calibration$source, "ftw")
})

test_that("add_ftw_prior preserves calibration when recalibration throws", {
  rs <- .make_test_rs()
  before <- rs$internals$calibration
  bb <- rs$geom$bbox
  field <- sf::st_as_sf(terra::project(
    terra::as.polygons(terra::ext(rs$internals$feature_stack),
                       crs = terra::crs(rs$internals$feature_stack)),
    "EPSG:4326"))
  field$source_year <- 2020L
  testthat::local_mocked_bindings(
    .ftwCacheRead = function(...) NULL,
    .ftwQuery = function(...) field,
    .ftwCacheWrite = function(...) invisible(NULL),
    .ftwCacheEvictLRU = function(...) invisible(0L),
    .calibrateScene = function(...) stop("synthetic calibration failure"),
    .package = "ocular")
  expect_warning(
    out <- suppressMessages(add_ftw_prior(rs, source = "local.parquet")),
    "preserved")
  expect_identical(out$internals$calibration, before)
  expect_false(is.null(out$geom$ftw_prior$polygon))
})

test_that("an undated FTW boundary calibrates but exerts no soft penalty", {
  rs <- .make_test_rs()
  rs$geom$ftw_prior <- list(
    polygon = sf::st_as_sf(terra::project(
      terra::as.polygons(terra::ext(rs$internals$feature_stack),
                         crs = terra::crs(rs$internals$feature_stack)),
      "EPSG:4326")),
    source_year = NULL)
  expect_equal(ocular:::.ftwPriorStrength(rs), 0)
  expect_identical(ocular:::.calibrateScene(rs)$source, "ftw")
})
