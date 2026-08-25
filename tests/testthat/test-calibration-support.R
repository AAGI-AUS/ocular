.make_support_rs <- function(nr, nc, centre_rc) {
  base <- matrix(seq_len(nr * nc) / (nr * nc + 1), nrow = nr, ncol = nc)
  template <- terra::rast(nrows = nr, ncols = nc, xmin = 0, xmax = nc,
                          ymin = 0, ymax = nr, crs = "EPSG:32750")
  w1 <- terra::setValues(template, as.vector(t(base)))
  w2 <- terra::setValues(template, as.vector(t(base + 0.01)))
  feature_stack <- c(w1, w2)
  names(feature_stack) <- c("w1", "w2")

  rs <- ocular:::.newOcular()
  rs$geom$nr <- as.integer(nr)
  rs$geom$nc <- as.integer(nc)
  rs$geom$centre_rc <- as.integer(centre_rc)
  rs$internals$feature_stack <- feature_stack
  rs$internals$feat_array <- NULL
  rs$spec$index_name <- "EVI2"
  rs$params <- rs_params()
  rs
}

test_that("self-calibration uses the largest complete seed-centred square", {
  centred <- ocular:::.liteSupport(
    .make_support_rs(9L, 15L, c(5L, 8L)), c_noise = 0)

  expect_identical(centred$r, 4L)
  expect_identical(nrow(centred$nb), 81L)
  expect_identical(
    names(centred$support),
    c("kind", "r", "r_max", "centre_row", "centre_col",
      "initial_radius", "width_px", "limit_width_px", "grid_nrow",
      "grid_ncol", "geometric_cell_count", "limit_rule",
      "calibration_schema_version", "stop_reason", "n_w"))
  expect_identical(centred$support$kind, "seed_neighbourhood")
  expect_identical(centred$support$r_max, 4L)
  expect_identical(centred$support$centre_row, 5L)
  expect_identical(centred$support$centre_col, 8L)
  expect_identical(centred$support$initial_radius, 2L)
  expect_identical(centred$support$width_px, 9L)
  expect_identical(centred$support$limit_width_px, 9L)
  expect_identical(centred$support$grid_nrow, 9L)
  expect_identical(centred$support$grid_ncol, 15L)
  expect_equal(centred$support$geometric_cell_count, 81)
  expect_identical(
    centred$support$limit_rule,
    "largest_complete_seed_centred_square")
  expect_identical(
    centred$support$calibration_schema_version,
    ocular:::.CALIBRATION_SCHEMA_VERSION)
  expect_identical(centred$support$stop_reason, "feature_grid_limit")
  expect_identical(centred$support$n_w, c(81L, 81L))

  off_centre <- ocular:::.liteSupport(
    .make_support_rs(9L, 15L, c(3L, 5L)), c_noise = 0)
  expect_identical(off_centre$r, 2L)
  expect_identical(off_centre$support$r_max, 2L)
  expect_identical(off_centre$support$width_px, 5L)
  expect_identical(nrow(off_centre$nb), 25L)
  expect_identical(off_centre$support$stop_reason, "feature_grid_limit")

  precise_rs <- .make_support_rs(9L, 15L, c(5L, 8L))
  precise_stack <- precise_rs$internals$feature_stack
  precise_stack[[1L]] <- terra::setValues(
    precise_stack[[1L]], rep(0.5, terra::ncell(precise_stack)))
  precise_stack[[2L]] <- terra::setValues(
    precise_stack[[2L]], rep(0.6, terra::ncell(precise_stack)))
  precise_rs$internals$feature_stack <- precise_stack
  precise <- ocular:::.liteSupport(precise_rs, c_noise = 0.01)
  expect_identical(precise$r, 2L)
  expect_identical(nrow(precise$nb), 25L)
  expect_true(all(precise$precise_w))
  expect_identical(precise$support$stop_reason, "precision_reached")
})

test_that("self-calibration handles a seed on the feature-grid edge", {
  rs <- .make_support_rs(9L, 15L, c(1L, 5L))
  support <- ocular:::.liteSupport(rs, c_noise = NA_real_)

  expect_identical(support$r, 0L)
  expect_identical(nrow(support$nb), 1L)
  expect_identical(support$support$r_max, 0L)
  expect_identical(support$support$initial_radius, 0L)
  expect_identical(support$support$width_px, 1L)
  expect_equal(support$support$geometric_cell_count, 1)
  expect_identical(
    support$support$stop_reason, "noise_floor_unavailable")

  calibration <- ocular:::.calibrateScene(rs)
  expect_identical(
    calibration$schema_version, ocular:::.CALIBRATION_SCHEMA_VERSION)
  expect_identical(calibration$support$r_max, 0L)
  expect_false("area_sensitivity" %in% names(calibration$baselines))
  expect_true(all(calibration$reason_window == "n_below_min"))
})

test_that("calibration does not use a cached feature subset without a stack", {
  rs <- .make_support_rs(9L, 15L, c(5L, 8L))
  rs$internals$feat_array <- array(0.5, dim = c(9L, 15L, 1L))
  rs$internals$feature_stack <- NULL

  expect_null(ocular:::.seedSupportRadiusLimit(rs))
  expect_null(ocular:::.seedNeighborhood(rs))
})

test_that("calibration support is grid-bounded and resolution-independent", {
  rs <- .make_support_rs(9L, 15L, c(5L, 8L))
  rs$geom$pixel_size_m <- 10
  at_ten_metres <- ocular:::.seedSupportRadiusLimit(rs)
  rs$geom$pixel_size_m <- 30
  at_thirty_metres <- ocular:::.seedSupportRadiusLimit(rs)

  expect_identical(at_ten_metres, 4L)
  expect_identical(at_thirty_metres, at_ten_metres)

  rs$geom$nr <- 10L
  expect_null(ocular:::.seedSupportRadiusLimit(rs))
})

test_that("a stale calibration schema forces recalibration", {
  rs <- .make_test_rs()
  rs$internals$calibration <- list(
    source = "self",
    strictness = rs$params$strictness,
    area_threshold = rs$params$area_threshold,
    index_name = rs$params$search_index,
    baselines = list(area_sensitivity = 0.12))
  replacement <- list(
    schema_version = ocular:::.CALIBRATION_SCHEMA_VERSION,
    source = "self",
    strictness = rs$params$strictness,
    area_threshold = rs$params$area_threshold,
    index_name = rs$params$search_index,
    baselines = list(area_sensitivity = 0.11))
  calls <- 0L
  testthat::local_mocked_bindings(
    .calibrateScene = function(...) {
      calls <<- calls + 1L
      replacement
    },
    .package = "ocular")

  out <- ocular:::.applyScenePriors(rs, rs$params)

  expect_identical(calls, 1L)
  expect_identical(out$internals$calibration, replacement)
})

test_that("explicit area sensitivity retains scalar ownership", {
  rs <- .make_test_rs()
  baseline <- 0.2
  rs$internals$calibration <- list(
    schema_version = ocular:::.CALIBRATION_SCHEMA_VERSION,
    source = "self",
    strictness = rs$params$strictness,
    area_threshold = rs$params$area_threshold,
    index_name = rs$params$search_index,
    baselines = list(area_sensitivity = baseline),
    tolerance_window = c(0.11, baseline),
    valid_window = c(TRUE, TRUE),
    mu_window = c(0.60, 0.62))

  seen <- numeric()
  testthat::local_mocked_bindings(
    .searchAreaWithWidens = function(...) {
      args <- list(...)
      seen <<- c(seen, args$area_sensitivity)
      matrix(FALSE, args$nr_full, args$nc_full)
    },
    .package = "ocular")
  capture_sensitivity <- function(expr) {
    seen <<- numeric()
    suppressWarnings(force(expr))
    seen
  }

  omitted <- capture_sensitivity(segment_area(rs))
  formal <- capture_sensitivity(
    segment_area(rs, area_sensitivity = baseline))
  partial <- capture_sensitivity(
    segment_area(rs, params = list(area_sensitivity = baseline)))
  stored_rs <- rs
  stored_rs$params <- rs_params(area_sensitivity = baseline)
  stored <- capture_sensitivity(segment_area(stored_rs))

  expect_equal(omitted, c(0.11, baseline))
  expect_equal(formal, rep(baseline, 2L))
  expect_equal(partial, rep(baseline, 2L))
  expect_equal(stored, rep(baseline, 2L))
})

test_that("delineation and fill interfaces have the canonical controls", {
  expect_identical(
    names(rs_params()),
    c("area_threshold", "area_sensitivity", "local_alive_density",
      "time_independent_windows", "search_windows", "search_index",
      "min_windows_alive", "search_start_date", "search_end_date",
      "perimeter_margins", "trace_iter", "multiple_areas",
      "merge_separate_areas", "area_separation_strict",
      "interior_sensitivity", "interior_threshold", "split_sensitivity",
      "split_gate", "split_linear", "strictness"))
  expect_identical(
    names(formals(segment_area)),
    c("x", "area_sensitivity", "area_threshold", "strictness",
      "local_alive_density", "time_independent_windows", "search_windows",
      "min_windows_alive", "perimeter_margins", "trace_iter",
      "merge_separate_areas", "params", "..."))
  expect_identical(
    names(formals(trace_perimeter)),
    c("x", "perimeter_margins", "trace_iter", "strictness", "params",
      "..."))
  expect_identical(
    names(formals(segment_interior)),
    c("x", "interior_sensitivity", "interior_threshold", "strictness",
      "params", "..."))
  expect_identical(
    names(formals(split_area)),
    c("x", "split_sensitivity", "split_gate", "strictness",
      "split_linear", "params", "..."))
  expect_identical(
    names(formals(boundary_delineation)),
    c("x", "vi_sensitivity", "vi_threshold", "strictness",
      "cleanup_boundary", "params", "..."))
  expect_identical(
    names(formals(ocular:::.floodFill)),
    c("feat_array", "centre_rc", "mu", "area_sensitivity",
      "area_threshold", "min_windows_alive", "reserve_rows",
      "reserve_cols", "nr", "nc", "local_alive_density", "prior_field",
      "prior_strength", "prior_cross"))
  expect_identical(
    names(formals(ocular:::.searchAreaWithWidens)),
    c("feat_array", "centre_rc", "mu", "area_sensitivity",
      "area_threshold", "min_windows_alive", "local_alive_density",
      "nr_full", "nc_full", "log_msg", "prior_field", "prior_strength",
      "prior_cross"))
})

test_that("flood fill can reach every corner of a matching rectangle", {
  nr <- 7L
  nc <- 11L
  filled <- ocular:::.floodFill(
    feat_array = array(0.5, dim = c(nr, nc, 1L)),
    centre_rc = c(4L, 6L),
    mu = 0.5,
    area_sensitivity = 0.1,
    area_threshold = c(0.2, 0.7),
    min_windows_alive = 1L,
    reserve_rows = c(1L, nr),
    reserve_cols = c(1L, nc),
    nr = nr,
    nc = nc,
    local_alive_density = 0)

  expect_true(all(filled$alive_mat))
  expect_true(all(filled$alive_mat[cbind(c(1L, 1L, nr, nr),
                                         c(1L, nc, 1L, nc))]))
})
