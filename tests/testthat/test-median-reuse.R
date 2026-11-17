test_that("splitting computes one stage-local median matrix", {
  original <- ocular:::.windowMedians
  calls <- 0L
  testthat::local_mocked_bindings(
    .windowMedians = function(feat_array, cells_mat) {
      calls <<- calls + 1L
      original(feat_array, cells_mat)
    }, .package = "ocular")
  area <- segment_area(.make_test_rs())
  split_area(area)
  expect_identical(calls, 1L)
  calls <- 0L
  split_area(area, split_sensitivity = NULL)
  expect_identical(calls, 0L)
  area$state$alive_mat[] <- FALSE
  split_area(area)
  expect_identical(calls, 0L)
})

test_that("an explicitly supplied NULL median cache is not recomputed", {
  calls <- 0L
  original <- ocular:::.windowMedians
  testthat::local_mocked_bindings(
    .windowMedians = function(feat_array, cells_mat) {
      calls <<- calls + 1L
      original(feat_array, cells_mat)
    }, .package = "ocular")
  area <- segment_area(.make_test_rs())
  rule <- ocular:::.splitAreaRule(cell_medians = NULL)
  expect_identical(rule$walker_setup(area, area$params), list(cell_medians = NULL))
  expect_identical(calls, 0L)
  direct <- ocular:::.splitAreaRule()
  expect_identical(direct$walker_setup(area, area$params)$cell_medians,
                   original(area$internals$feat_array, area$state$alive_mat))
  expect_identical(calls, 1L)
})
