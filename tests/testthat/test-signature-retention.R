.signature_retention_rs <- function(current = TRUE, missing = 2L, values = c(0.4, 0.6, 0.3)){

  n <- 21L
  field <- matrix(FALSE, n, n)
  field[5:17, 5:17] <- TRUE
  template <- terra::rast(nrows = n, ncols = n, xmin = 0, xmax = n,
                          ymin = 0, ymax = n, crs = "EPSG:32750")
  planes <- lapply(values, function(value){

    plane <- matrix(0.1, n, n)
    plane[field] <- value
    return(plane)
  })
  for( window in missing ){
    planes[[window]][11L, 11L] <- NA_real_
  }
  rs <- ocular:::.newOcular()
  rs$spec$index_name <- "EVI2"
  rs$geom$nr <- n
  rs$geom$nc <- n
  rs$geom$centre_rc <- c(11L, 11L)
  rs$params <- ocular::rs_params()
  rs$internals$feature_stack <- do.call(c, lapply(planes, function(plane){

    return(terra::setValues(template, as.vector(t(plane))))
  }))
  rs$internals$calibration <- ocular:::.calibrateScene(rs)
  rs$state$decided_mat <- matrix(FALSE, n, n)
  if( current ){
    rs$state$alive_mat <- field
  }
  return(rs)
}

test_that("eligible calibration retains a signature for a current-mask fill", {
  for( independent in c(TRUE, FALSE) ){
    rs <- .signature_retention_rs()
    values <- terra::values(rs$internals$feature_stack)
    expected <- rs$state$alive_mat
    expected[11L, 11L] <- FALSE
    filled <- segment_area(rs, search_windows = 2L,
                           time_independent_windows = independent)
    expect_identical(filled$state$alive_mat, expected)
    expect_identical(filled$internals$mu, 0.6)
    expect_identical(filled$params, rs$params)
    expect_identical(filled$internals$calibration, rs$internals$calibration)
    expect_identical(terra::values(filled$internals$feature_stack), values)
    expect_false(any(filled$state$decided_mat))
    interior <- ocular::segment_interior(filled)
    expect_identical(interior$state$alive_mat, expected)
    expect_null(attr(interior, ".calibration_context", exact = TRUE))
    expect_identical(segment_area(filled, search_windows = 2L)$state,
                              filled$state)
  }
})

test_that("retained signatures follow selected order and preserve learned finite values", {
  rs <- .signature_retention_rs()
  plane <- terra::as.matrix(rs$internals$feature_stack[[1L]], wide = TRUE)
  plane[10:12, 10:12] <- 0.41
  rs$internals$feature_stack[[1L]] <- terra::setValues(
    rs$internals$feature_stack[[1L]], as.vector(t(plane)))
  rs$internals$calibration <- ocular:::.calibrateScene(rs)
  expect_identical(rs$internals$calibration$mu_window[1L], 0.4)
  selected <- segment_area(rs, search_windows = c(2L, 1L))
  expect_identical(selected$internals$mu, c(0.6, 0.41))
  reversed <- segment_area(rs, search_windows = c(1L, 2L))
  expect_identical(reversed$internals$mu, c(0.41, 0.6))
  expect_identical(reversed$state, selected$state)
})

test_that("explicit scalar sensitivity and invalid windows retain missing signatures", {
  rs <- .signature_retention_rs()
  baseline <- rs$internals$calibration$baselines$area_sensitivity
  for( setting in list(list(area_sensitivity = baseline),
                       list(params = list(area_sensitivity = baseline)),
                       list(params = ocular::rs_params(area_sensitivity = baseline))) ){
    out <- do.call(segment_area, c(list(rs, search_windows = c(2L, 1L)), setting))
    expect_identical(out$internals$mu, c(NA_real_, 0.4))
  }
  rs$params <- ocular::rs_params(area_sensitivity = baseline)
  expect_identical(segment_area(rs, search_windows = c(2L, 1L))$internals$mu,
                            c(NA_real_, 0.4))
  rs <- .signature_retention_rs()
  rs$internals$feature_stack[[2L]] <- terra::setValues(
    rs$internals$feature_stack[[2L]], rep(NA_real_, 441L))
  rs$internals$calibration <- ocular:::.calibrateScene(rs)
  out <- segment_area(rs, search_windows = c(2L, 1L))
  expect_identical(out$internals$mu, c(NA_real_, 0.4))
})

test_that("first fills retain empty masks and use selected signatures in later refinement", {
  rs <- .signature_retention_rs(current = FALSE)
  expect_warning(empty <- segment_area(rs, search_windows = 2L),
                          "no candidate field was identified")
  expect_false(any(empty$state$alive_mat))
  expect_identical(empty$internals$mu, 0.6)
  expected <- matrix(FALSE, 21L, 21L)
  expected[5:17, 5:17] <- TRUE
  filled <- segment_area(rs, search_windows = c(2L, 1L))
  expect_identical(filled$state$alive_mat, expected)
  expect_identical(ocular::segment_interior(
    filled, interior_sensitivity = 0.05)$state$alive_mat, expected)
})

test_that("retained calibration enables the existing merge predicate", {
  rs <- .signature_retention_rs()
  rs$state$alive_mat[, 10L] <- FALSE
  unmerged <- segment_area(rs, search_windows = 2L, multiple_areas = TRUE)
  merged <- segment_area(rs, search_windows = 2L, multiple_areas = TRUE,
                         merge_separate_areas = TRUE)
  bridge <- matrix(FALSE, 21L, 21L)
  bridge[5:17, 10L] <- TRUE
  expect_identical(merged$state$alive_mat,
                            unmerged$state$alive_mat | bridge)
  expect_identical(merged$state$decided_mat, bridge)
  single <- segment_area(rs, search_windows = 2L, multiple_areas = FALSE,
                         merge_separate_areas = TRUE)
  expect_identical(single$state, merged$state)
})

test_that("an unknown validity flag does not enter the replacement indices", {
  rs <- .signature_retention_rs(missing = c(2L, 3L))
  rs$internals$calibration$valid_window[3L] <- NA
  filled <- segment_area(rs, search_windows = c(1L, 2L, 3L),
                         time_independent_windows = FALSE)
  expect_identical(filled$internals$mu, c(0.4, 0.6, NA_real_))
  expect_identical(filled$state$alive_mat, rs$state$alive_mat)
})


test_that("default boundary refinement uses retained selected signatures", {
  for( current in c(TRUE, FALSE) ){
    rs <- .signature_retention_rs(current = current, values = c(0.3, 0.6))
    expected <- matrix(FALSE, 21L, 21L)
    expected[5:17, 5:17] <- TRUE
    if( current ){
      rs$state$alive_mat[] <- FALSE
      rs$state$alive_mat[7:15, 7:15] <- TRUE
      expected <- rs$state$alive_mat
    }
    out <- boundary_delineation(rs)
    expect_identical(out$state$alive_mat, expected)
    expect_identical(out$internals$mu, c(0.3, 0.6))
    expect_null(attr(out, ".calibration_context", exact = TRUE))
  }
})

test_that("signature retention uses refreshed calibration and excludes non-finite values", {
  for( stale in c(TRUE, FALSE) ){
    rs <- .signature_retention_rs()
    expected <- rs$state$alive_mat
    expected[11L, 11L] <- FALSE
    if( stale ){
      rs$internals$calibration$schema_version <- 0L
    }else{
      rs$internals$calibration <- NULL
    }
    out <- segment_area(rs, search_windows = 2L)
    expect_identical(out$internals$mu, 0.6)
    expect_identical(out$state$alive_mat, expected)
  }
  for( missing in c(NA_real_, NaN, Inf, -Inf) ){
    rs <- .signature_retention_rs()
    plane <- terra::as.matrix(rs$internals$feature_stack[[2L]], wide = TRUE)
    plane[11L, 11L] <- missing
    rs$internals$feature_stack[[2L]] <- terra::setValues(
      rs$internals$feature_stack[[2L]], as.vector(t(plane)))
    out <- segment_area(rs, search_windows = c(2L, 1L))
    expect_identical(out$internals$mu, c(0.6, 0.4))
    rs$internals$calibration$mu_window[2L] <- missing
    out <- segment_area(rs, search_windows = c(2L, 1L))
    expect_identical(out$internals$mu, c(NA_real_, 0.4))
  }
})
