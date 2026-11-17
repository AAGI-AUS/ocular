test_that("FTW diagnostics distinguish absent references, missing CRS and no overlap", {
  rs <- .make_test_rs()
  rs$state$alive_mat <- matrix(TRUE, rs$geom$nr, rs$geom$nc)
  template <- rs$internals$feature_stack[[1L]]
  reference <- sf::st_as_sf(terra::as.polygons(terra::ext(template),
                                              crs = terra::crs(template)))
  shifted <- sf::st_set_crs(sf::st_geometry(reference) + c(50000, 50000),
                            sf::st_crs(reference))
  expect_error(diagnose_against_ftw(rs), "no reference polygon")

  no_crs <- sf::st_set_crs(reference, NA)
  vector_no_crs <- terra::vect(reference)
  terra::crs(vector_no_crs) <- ""
  for( polygon in list(no_crs, sf::st_geometry(no_crs), vector_no_crs) ){
    expect_error(diagnose_against_ftw(rs, polygon), "reference polygon has no CRS")
    attached <- rs
    attached$geom$ftw_prior <- list(polygon = polygon)
    expect_error(diagnose_against_ftw(attached), "reference polygon has no CRS")
  }

  for( polygon in list(shifted, sf::st_as_sf(shifted), terra::vect(shifted),
                       sf::st_transform(shifted, 4326), c(shifted, shifted)) ){
    expect_error(diagnose_against_ftw(rs, polygon), "does not overlap the feature grid")
    attached <- rs
    attached$geom$ftw_prior <- list(polygon = polygon)
    expect_error(diagnose_against_ftw(attached), "does not overlap the feature grid")
  }

  expect_error(diagnose_against_ftw(rs, 1), "reference polygon could not be used")
  grid <- as.vector(terra::ext(template))
  tiny <- terra::as.polygons(
    terra::ext(grid[1L] + 0.1, grid[1L] + 0.2,
               grid[3L] + 0.1, grid[3L] + 0.2), crs = terra::crs(template))
  expect_error(diagnose_against_ftw(rs, tiny), "reference polygon could not be used")
})

test_that("FTW diagnostic guards preserve successful agreement calculations", {
  rs <- .make_test_rs()
  rs$state$alive_mat <- matrix(TRUE, rs$geom$nr, rs$geom$nc)
  template <- rs$internals$feature_stack[[1L]]
  reference <- sf::st_as_sf(terra::as.polygons(terra::ext(template),
                                              crs = terra::crs(template)))
  shifted <- sf::st_set_crs(sf::st_geometry(reference) + c(50000, 50000),
                            sf::st_crs(reference))
  expected <- list(iou = 1, f1 = 1, precision = 1, recall = 1,
                   n_pred = 441L, n_ref = 441L, n_intersection = 441L, n_union = 441L)
  state <- rs$state
  for( polygon in list(reference, sf::st_geometry(reference), terra::vect(reference),
                       sf::st_transform(reference, 4326),
                       sf::st_as_sf(c(sf::st_geometry(reference), shifted))) ){
    expect_identical(diagnose_against_ftw(rs, polygon), expected)
    attached <- rs
    attached$geom$ftw_prior <- list(polygon = polygon)
    expect_identical(diagnose_against_ftw(attached), expected)
  }
  expect_identical(rs$state, state)

  no_crs <- rs
  no_crs$internals$feature_stack <- terra::rast(rs$internals$feature_stack)
  terra::crs(no_crs$internals$feature_stack) <- ""
  reference_no_crs <- terra::vect(reference)
  terra::crs(reference_no_crs) <- ""
  expect_identical(diagnose_against_ftw(no_crs, reference_no_crs), expected)

  rs$state$alive_mat[,] <- FALSE
  rs$state$alive_mat[, rs$geom$nc] <- TRUE
  grid <- as.vector(terra::ext(template))
  left <- terra::as.polygons(
    terra::ext(grid[1L], grid[1L] + 100, grid[3L], grid[4L]),
    crs = terra::crs(template))
  expect_identical(diagnose_against_ftw(rs, left),
                   list(iou = 0, f1 = 0, precision = 0, recall = 0,
                        n_pred = 21L, n_ref = 210L, n_intersection = 0L, n_union = 231L))
})

test_that("fusion validation directs unsupported ocular input to supported outputs", {
  rs <- .make_test_rs()
  testthat::local_mocked_bindings(
    .validateBuildRs = function(...) stop("unexpected retrieval"),
    .package = "ocular")
  msg <- "ocular objects are not supported directly.*as_raster.*daily as_time_series"
  expect_error(validate_data_fusion(rs), msg)
  expect_error(validate_data_fusion(
    rs, latitude = -32, start_date = "2021-01-01",
    end_date = "2021-02-01", x_metres = 100), msg)
})
