test_that("public surface is exported and callable", {
  exported <- c(
    "get_rs", "segment_area", "trace_perimeter", "segment_interior",
    "split_area", "boundary_delineation", "add_modis", "as_raster",
    "as_time_series", "validate_data_fusion", "build_pipeline", "rs_params",
    "is_rs", "rs_utm", "metres_to_degrees", "point_to_bbox", "add_ftw_prior",
    "batch_rs", "clear_ftw_cache", "clear_rs_cache",
    "diagnose_against_ftw"
  )
  for (fn in exported) {
    expect_true(is.function(getExportedValue("ocular", fn)),
                info = paste0(fn, " should be an exported function"))
  }
})

test_that("rs_params() returns a non-empty named list of defaults", {
  p <- rs_params()
  expect_type(p, "list")
  expect_gt(length(p), 0L)
  expect_false(is.null(names(p)))
  expect_false(any(names(p) == ""))
  expect_identical(p$search_windows, "all")
})

test_that("the installed namespace and package identity are lowercase ocular", {
  desc <- utils::packageDescription("ocular")
  expect_identical(desc$Package, "ocular")
  expect_identical(desc$Version, "0.1.0")
  expect_identical(desc$URL, "https://github.com/AAGI-AUS/ocular")
  expect_identical(desc$BugReports,
                   "https://github.com/AAGI-AUS/ocular/issues")
  expect_identical(environmentName(environment(get_rs)), "ocular")
  expect_true("ocular" %in% loadedNamespaces())
})
