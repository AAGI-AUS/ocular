.rerun_mask_rs <- function(){

  rs <- .make_test_rs()
  mask <- matrix(FALSE, 21L, 21L)
  mask[7:15, 7:15] <- TRUE
  rs$state$alive_mat <- mask
  rs
}

test_that("area filling uses the current mask with fresh invocation state", {
  rs <- .rerun_mask_rs()
  rs$state$alive_mat[11L, 11L] <- FALSE
  clean <- rs
  rs$state$decided_mat[] <- TRUE
  attr(rs, "pending") <- list(pre_cleanup = list(
    active = TRUE, thresholds = c(0.9, 1), iter = 20L))

  out <- expect_silent(segment_area(rs))
  fresh <- expect_silent(segment_area(clean))
  expect_identical(out$state$alive_mat, clean$state$alive_mat)
  expect_identical(out$state, fresh$state)
  expect_false(any(out$state$decided_mat))
  expect_false(isTRUE(attr(out, "pending")$pre_cleanup$active))
  expect_identical(out$params, rs$params)
  expect_identical(out$geom, rs$geom)
  expect_identical(out$scenes, rs$scenes)
  expect_identical(out$internals$feature_stack, rs$internals$feature_stack)
  expect_true(all(is.finite(out$internals$feat_array)))
  expect_identical(rs$state$alive_mat, clean$state$alive_mat)
  expect_true(all(rs$state$decided_mat))
})

test_that("rerun settings still drive spectral fill and clear empty-result marks", {
  rs <- .rerun_mask_rs()
  rs$state$decided_mat[] <- TRUE
  keep <- segment_area(rs, area_threshold = c(0.5, 0.7))
  expect_identical(keep$state$alive_mat, rs$state$alive_mat)
  expect_warning(
    empty <- segment_area(rs, area_threshold = c(0.8, 1)),
    "no candidate field was identified", fixed = TRUE)
  expect_false(any(empty$state$alive_mat))
  expect_false(any(empty$state$decided_mat))
  expect_identical(empty$params, rs$params)
  expect_warning(again <- segment_area(empty), "no candidate field was identified")
  expect_identical(again$state$alive_mat, empty$state$alive_mat)
})

test_that("mask fill handles separate areas and a missing original seed", {
  rs <- .rerun_mask_rs()
  rs$state$alive_mat[] <- FALSE
  rs$state$alive_mat[6:9, 6:9] <- TRUE
  rs$state$alive_mat[13:16, 13:16] <- TRUE
  for (joint in c(FALSE, TRUE)) {
    out <- segment_area(rs, multiple_areas = TRUE,
                        time_independent_windows = !joint)
    expect_identical(out$state$alive_mat, rs$state$alive_mat)
    one <- segment_area(rs, multiple_areas = FALSE,
                        time_independent_windows = !joint)
    expected <- matrix(FALSE, 21L, 21L)
    expected[6:9, 6:9] <- TRUE
    expect_identical(one$state$alive_mat, expected)
  }
})

test_that("boundary reruns ignore inherited bookkeeping and keep every stage active", {
  rs <- .rerun_mask_rs()
  clean <- rs
  rs$state$decided_mat[] <- TRUE
  attr(rs, "pending") <- list(pre_cleanup = list(
    active = TRUE, thresholds = c(0.9, 1), iter = 20L))
  first <- expect_silent(boundary_delineation(rs))
  fresh <- expect_silent(boundary_delineation(clean))
  expect_identical(first$state, fresh$state)
  second <- expect_silent(boundary_delineation(first))
  third <- expect_silent(boundary_delineation(second))
  expect_identical(second$state$alive_mat, third$state$alive_mat)
  expect_identical(first$params, rs$params)
  expect_identical(first$internals$feature_stack, rs$internals$feature_stack)

  seen <- character()
  wrap <- function(name, fn){
    force(name)
    force(fn)
    function(...) {
      seen <<- c(seen, name)
      fn(...)
    }
  }
  testthat::local_mocked_bindings(
    trace_perimeter = wrap("cleanup", trace_perimeter),
    segment_area = wrap("fill", segment_area),
    split_area = wrap("split", split_area),
    segment_interior = wrap("interior", segment_interior),
    .package = "ocular")
  boundary_delineation(first)
  expect_identical(as.integer(table(factor(seen,
    levels = c("cleanup", "fill", "split", "interior")))), c(8L, 1L, 5L, 1L))
})


test_that("refilling a split mask consumes its cuts without inheriting locks", {
  rs <- .make_test_rs(n = 41L)
  alive <- matrix(FALSE, 41L, 41L)
  alive[13:27, 10:30] <- TRUE
  feature <- matrix(0.1, 41L, 41L)
  feature[13:27, 10:19] <- 0.7
  feature[13:27, 20:30] <- 0.2
  layer <- terra::setValues(rs$internals$feature_stack[[1L]], as.vector(t(feature)))
  rs$internals$feature_stack <- layer
  rs$internals$feat_array <- array(feature, c(41L, 41L, 1L))
  rs$internals$mu <- 0.7
  rs$internals$calibration <- NULL
  rs$geom$centre_rc <- c(20L, 15L)
  rs$state$alive_mat <- alive
  cut <- split_area(rs, split_sensitivity = 0.05, multiple_areas = TRUE)
  removed <- alive & !cut$state$alive_mat
  expect_true(any(removed))
  expect_true(any(cut$state$decided_mat & removed))
  out <- segment_area(cut, area_sensitivity = 0.8, multiple_areas = TRUE)
  expect_false(any(out$state$alive_mat & removed))
  expect_false(any(out$state$decided_mat))
  fresh <- cut
  fresh$state$decided_mat[] <- FALSE
  expect_identical(out$state, segment_area(
    fresh, area_sensitivity = 0.8, multiple_areas = TRUE)$state)
})
