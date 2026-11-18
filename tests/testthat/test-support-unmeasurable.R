.make_sparse_support_rs <- function(usable = 0L, distant = FALSE) {
  nr <- 15L
  template <- terra::rast(nrows = nr, ncols = nr, xmin = 0, xmax = nr,
                          ymin = 0, ymax = nr, crs = "EPSG:32750")
  good <- terra::setValues(template, rep(0.5, nr * nr))
  values <- rep(NA_real_, nr * nr)
  if (usable > 0L) {
    cells <- if (distant) seq_len(usable) else 113L + seq_len(usable) - 1L
    values[cells] <- 0.6
  }
  bad <- terra::setValues(template, values)
  rs <- ocular:::.newOcular()
  rs$geom$nr <- nr
  rs$geom$nc <- nr
  rs$geom$centre_rc <- c(8L, 8L)
  rs$internals$feature_stack <- c(good, bad)
  rs$spec$index_name <- "EVI2"
  rs$params <- rs_params()
  rs
}

test_that("unmeasurable windows do not widen shared seed support", {
  for (usable in c(0L, 1L)) {
    rs <- .make_sparse_support_rs(usable)
    support <- ocular:::.liteSupport(rs, c_noise = 0.01)
    expect_identical(support$r, 2L)
    expect_identical(support$support$stop_reason, "precision_reached")
    expect_identical(support$precise_w, c(TRUE, FALSE))
    expect_identical(support$sig$n_w, c(25L, usable))
    calibration <- ocular:::.calibrateScene(rs)
    expect_identical(calibration$valid_window, c(TRUE, FALSE))
    expect_identical(calibration$reason_window[2L], "n_below_min")
  }
})

test_that("measurable distant windows retain the support growth rule", {
  rs <- .make_sparse_support_rs(2L, distant = TRUE)
  support <- ocular:::.liteSupport(rs, c_noise = 0.01)
  expect_identical(support$r, 7L)
  expect_identical(support$support$stop_reason, "precision_reached")
  expect_identical(support$sig$n_w, c(225L, 2L))
  expect_true(all(support$precise_w))

  unknown_noise <- ocular:::.liteSupport(rs, c_noise = NA_real_)
  expect_identical(unknown_noise$r, 2L)
  expect_identical(unknown_noise$support$stop_reason, "noise_floor_unavailable")
})

test_that("support records when every window is unmeasurable", {
  rs <- .make_sparse_support_rs()
  rs$internals$feature_stack[[1L]] <- terra::setValues(
    rs$internals$feature_stack[[1L]], rep(NA_real_, 225L))
  support <- ocular:::.liteSupport(rs, c_noise = 0.01)
  expect_identical(support$r, 2L)
  expect_identical(support$support$stop_reason, "no_measurable_windows")
  expect_false(any(support$precise_w))
  expect_true(all(ocular:::.calibrateScene(rs)$reason_window == "n_below_min"))
})

test_that("saved calibration from before the support correction is refreshed", {
  rs <- .make_sparse_support_rs()
  expected <- ocular:::.calibrateScene(rs)
  rs$internals$calibration <- expected
  rs$internals$calibration$schema_version <- 2L
  rs$internals$calibration$mu_window <- c(0.1, NA_real_)
  rs$internals$calibration$support$r <- 7L
  out <- ocular:::.applyScenePriors(rs, rs$params)
  expect_identical(out$internals$calibration, expected)
  expect_identical(out$internals$calibration$support$r, 2L)
})
