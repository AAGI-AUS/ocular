.s2_processing_item <- function(id = "high", baseline = "05.10",
                                generation = "2024-07-26T04:59:51Z"){
  list(id = id, properties = list(
    datetime = "2022-09-29T00:11:21.000000Z", platform = "Sentinel-2A",
    `s2:mgrs_tile` = "55HED", `sat:relative_orbit` = 73L,
    `s2:datastrip_id` = paste0("S2A_OPER_MSI_L2A_DS_TEST_20240726T045951_",
                              "S20220929T001116_N", baseline),
    `s2:processing_baseline` = baseline, `s2:generation_time` = generation),
    assets = setNames(lapply(c("B04", "B08", "SCL"), function(b){
      list(href = paste0("fixture-", id, "_", b))
    }), c("B04", "B08", "SCL")))
}

.s2_processing_ids <- function(features, policy = "highest_baseline"){
  got <- suppressWarnings(ocular:::.selectS2Processing(features, policy))
  vapply(got, `[[`, character(1L), "id")
}

test_that("Sentinel-2 baseline selection is validated and retrieval-only", {
  expect_identical(rs_params()$s2_processing, "highest_baseline")
  expect_identical(rs_params(s2_processing = "lowest_baseline")$s2_processing,
                   "lowest_baseline")
  for( bad in list(NULL, NA, "latest", 1, c("highest_baseline", "lowest_baseline")) ){
    expect_error(rs_params(s2_processing = bad), "s2_processing")
  }
  rs <- .make_test_rs()
  expect_error(ocular:::.resolveParams(rs, list(s2_processing = "lowest_baseline"), list()),
               "fetch-time parameter.*s2_processing")
  expect_identical(ocular:::.resolveParams(rs, NULL, list())$s2_processing,
                   "highest_baseline")
})

test_that("highest and lowest baselines are selected independently of provider order", {
  low <- .s2_processing_item("low", "04.00", "2025-01-01T00:00:00Z")
  high <- .s2_processing_item()
  for( items in list(list(low, high), list(high, low)) ){
    expect_identical(.s2_processing_ids(items), "high")
    expect_identical(.s2_processing_ids(items, "lowest_baseline"), "low")
    expect_warning(ocular:::.selectS2Processing(items, "highest_baseline"),
                   'highest available processing baseline.*05.10.*lowest_baseline')
    expect_warning(ocular:::.selectS2Processing(items, "lowest_baseline"),
                   'lowest available processing baseline.*04.00.*highest_baseline')
  }
  decimal <- .s2_processing_item("decimal", "05.9")
  expect_identical(.s2_processing_ids(list(decimal, high)), "high")
})

test_that("generation time and item ID resolve equal-baseline ties", {
  old <- .s2_processing_item("old", "05.10", "2023-01-01T00:00:00Z")
  new <- .s2_processing_item("new", "05.10", "2024-01-01T00:00:00Z")
  expect_identical(.s2_processing_ids(list(old, new)), "new")
  expect_identical(.s2_processing_ids(list(old, new), "lowest_baseline"), "old")
  same <- new
  same$id <- "a"
  for( items in list(list(new, same), list(same, new)) ){
    expect_identical(.s2_processing_ids(items), "a")
    expect_identical(.s2_processing_ids(items, "lowest_baseline"), "a")
  }
  missing <- new
  missing$id <- "missing"
  missing$properties$`s2:generation_time` <- NULL
  expect_identical(.s2_processing_ids(list(missing, old)), "old")
  expect_identical(.s2_processing_ids(list(missing, old), "lowest_baseline"), "old")

  new$properties$`s2:generation_time` <- "2024-01-01T09:30:00+09:30"
  old$properties$`s2:generation_time` <- "2023-01-01T00:00:00+00:00"
  expect_identical(.s2_processing_ids(list(old, new)), "new")
  expect_identical(.s2_processing_ids(list(old, new), "lowest_baseline"), "old")

  old$properties$`s2:generation_time` <- "2024-01-01T00:00:00.000000001Z"
  new$properties$`s2:generation_time` <- "2024-01-01T00:00:00.000000002Z"
  old$id <- "a-old"
  new$id <- "z-new"
  expect_identical(.s2_processing_ids(list(old, new)), "z-new")
  expect_identical(.s2_processing_ids(list(old, new), "lowest_baseline"), "a-old")
})

test_that("same-date tiles, sensing times, platforms and datastrips remain distinct", {
  base <- .s2_processing_item()
  changes <- list(platform = "Sentinel-2B", datetime = "2022-09-29T00:30:00Z",
                  `s2:mgrs_tile` = "55HFD", `sat:relative_orbit` = 74L,
                  `s2:datastrip_id` = "S2A_OPER_MSI_L2A_DS_TEST_S20220929T001200_N05.10")
  for( nm in names(changes) ){
    other <- base
    other$id <- nm
    other$properties[[nm]] <- changes[[nm]]
    expect_identical(ocular:::.selectS2Processing(list(base, other), "highest_baseline"),
                     list(base, other))
  }
})

test_that("incomplete identity is not grouped and unknown baseline does not invent a rank", {
  base <- .s2_processing_item()
  for( nm in c("platform", "datetime", "s2:mgrs_tile", "sat:relative_orbit",
                "s2:datastrip_id") ){
    missing <- base
    missing$properties[[nm]] <- NULL
    expect_warning(got <- ocular:::.selectS2Processing(list(base, missing), "highest_baseline"),
                     "acquisition identity was incomplete")
    expect_identical(got, list(base, missing))
  }
  unknown <- base
  unknown$id <- "unknown"
  unknown$properties$`s2:processing_baseline` <- "unreadable"
  expect_warning(got <- ocular:::.selectS2Processing(list(base, unknown), "highest_baseline"),
                   "retained without processing selection")
  expect_identical(got, list(base, unknown))
  expect_warning(got <- ocular:::.selectS2Processing(list(list()), "highest_baseline"),
                   "acquisition identity was incomplete")
  expect_identical(got, list(list()))
})

test_that("processing selection precedes reads and prevents backfilling a rejected cloud cell", {
  original_rast <- terra::rast
  r <- original_rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2,
                     ymin = 0, ymax = 2, crs = "EPSG:4326")
  rasters <- list(high_B04 = terra::setValues(r, 3000),
                  high_B08 = terra::setValues(r, 7000),
                  high_SCL = terra::setValues(r, c(9, 4, 4, 4)),
                  low_B04 = terra::setValues(r, 2000),
                  low_B08 = terra::setValues(r, 5000),
                  low_SCL = terra::setValues(r, 4))
  items <- list(features = list(.s2_processing_item(),
                                .s2_processing_item("low", "04.00")))
  reads <- character(0L)
  testthat::local_mocked_bindings(rast = function(x, ...){
    if( missing(x) ){
      return(original_rast(...))
    }
    if( is.character(x) && length(x) == 1L && startsWith(x, "/vsicurl/fixture-") ){
      key <- sub("/vsicurl/fixture-", "", x, fixed = TRUE)
      reads <<- c(reads, key)
      return(rasters[[key]])
    }
    return(original_rast(x, ...))
  }, .package = "terra")
  testthat::local_mocked_bindings(.fetchStacItems = function(...){ return(items) },
                                  .package = "ocular")
  run <- function(policy){
    ocular:::.fetchStac(sf::st_bbox(c(xmin = 0, ymin = 0, xmax = 2, ymax = 2), crs = 4326),
                        "2022-09-29", "2022-10-01", "NDVI", 50,
                        scl_classes = 9L, source = "sentinel-2", s2_processing = policy)
  }
  expect_warning(high <- run("highest_baseline"), "1 Sentinel-2 acquisition")
  expect_equal(terra::values(high[[1L]]$index, mat = FALSE), c(NA, 0.5, 0.5, 0.5))
  expect_identical(as.logical(terra::values(high[[1L]]$pixel_mask, mat = FALSE)),
                   c(TRUE, FALSE, FALSE, FALSE))
  expect_identical(reads, c("high_B08", "high_B04", "high_SCL"))
  expect_identical(high[[1L]]$item_assets$B04, "fixture-high_B04")
  expect_identical(high[[1L]]$dn_offset$B04, 1000)
  reads <- character(0L)
  expect_warning(low <- run("lowest_baseline"), "lowest available")
  expect_equal(terra::values(low[[1L]]$index, mat = FALSE), rep(0.6, 4))
  expect_identical(reads, c("low_B08", "low_B04", "low_SCL"))
})

test_that("get_rs forwards processing selection and batch cache keys distinguish it", {
  captured <- NULL
  testthat::local_mocked_bindings(.fetchStac = function(..., s2_processing){
    captured <<- s2_processing
    return(NULL)
  }, .package = "ocular")
  expect_error(get_rs(147.227906, -33.0675167, "2022-08-01", "2022-10-01",
                      source = "sentinel-2", x_metres = 1000,
                      params = rs_params(s2_processing = "lowest_baseline")),
               "No usable scenes")
  expect_identical(captured, "lowest_baseline")
  args <- list(longitude = 147.227906, latitude = -33.0675167,
                start_date = "2022-08-01", end_date = "2022-10-01", index_name = "NDVI",
                source = "sentinel-2", x_metres = 1000, y_metres = 1000, max_cloud_cover = 50,
                retries = 3L, cloud_relax = 70, pipeline_cache_key = "test",
                package_version = "0.1.1.9000", params = rs_params())
  highest <- do.call(ocular:::.rsCacheKey, args)
  args$params <- rs_params(s2_processing = "lowest_baseline")
  expect_false(identical(do.call(ocular:::.rsCacheKey, args), highest))
  args$params <- rs_params()
  args$params$s2_processing <- NULL
  expect_false(identical(do.call(ocular:::.rsCacheKey, args), highest))
})

test_that("acquisition timestamps compare equivalent instants without rounding fractions", {
  high <- .s2_processing_item()
  low <- .s2_processing_item("low", "04.00")
  for( stamp in c("2022-09-29T00:11:21Z", "2022-09-29T00:11:21+00:00",
                  "2022-09-29T09:41:21+09:30", "2022-09-28T18:11:21-06:00") ){
    low$properties$datetime <- stamp
    expect_identical(.s2_processing_ids(list(high, low)), "high")
  }
  high$properties$datetime <- "2022-09-29T00:11:21.123400Z"
  low$properties$datetime <- "2022-09-29T09:41:21.1234+09:30"
  expect_identical(.s2_processing_ids(list(high, low)), "high")
  high$properties$datetime <- "2022-09-29T00:11:21.000000001Z"
  low$properties$datetime <- "2022-09-29T00:11:21.000000002Z"
  expect_identical(ocular:::.selectS2Processing(list(high, low), "highest_baseline"),
                   list(high, low))
})
