test_that("the final fill window reaches every grid edge from an off-centre seed", {
  cases <- list(c(34L, 34L, 18L, 18L), c(35L, 35L, 17L, 17L),
                c(20L, 31L, 3L, 25L), c(7L, 9L, 2L, 8L))
  for (case in cases) {
    nr <- case[1L]
    nc <- case[2L]
    out <- ocular:::.searchAreaWithWidens(
      array(0.5, dim = c(nr, nc, 2L)), case[3:4], c(0.5, 0.5),
      area_sensitivity = 0.05, area_threshold = c(0.3, 0.9),
      min_windows_alive = 1L, local_alive_density = 0,
      nr_full = nr, nc_full = nc)
    expect_identical(out, matrix(TRUE, nr, nc), info = paste(case, collapse = ","))
  }
})

test_that("earlier fill windows and their early stopping remain centred", {
  nr <- 34L
  nc <- 34L
  centre <- c(18L, 18L)
  seen <- list()
  testthat::local_mocked_bindings(
    .floodFill = function(feat_array, centre_rc, nr, nc, ...) {
      seen[[length(seen) + 1L]] <<- list(nr = nr, nc = nc, centre = centre_rc)
      list(alive_mat = matrix(TRUE, nr, nc), n_flagged = 0L)
    }, .package = "ocular")
  out <- ocular:::.searchAreaWithWidens(
    array(0.5, dim = c(nr, nc, 2L)), centre, c(0.5, 0.5),
    0.05, c(0.3, 0.9), 1L, 0, nr, nc)
  expect_length(seen, 1L)
  expect_identical(seen[[1L]], list(nr = 11L, nc = 11L, centre = c(6L, 6L)))
  expected <- matrix(FALSE, nr, nc)
  expected[13:23, 13:23] <- TRUE
  expect_identical(out, expected)
})
