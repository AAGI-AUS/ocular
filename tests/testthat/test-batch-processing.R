test_that("batch cache keys are safe and include result-changing inputs", {
  p <- rs_params()
  key <- ocular:::.rsCacheKey(117.8, -32, "2021-01-01", "2021-02-01",
                               "../../NDVI", "auto", 100, 100, 50, p,
                               3L, 80, "default")
  expect_match(key, "^[0-9a-f]{24}\\.rds$")
  expect_false(grepl("[/\\\\]", key))
  expect_false(identical(
    key,
    ocular:::.rsCacheKey(117.8, -32, "2021-01-01", "2021-02-01",
                          "../../NDVI", "auto", 100, 100, 50,
                          rs_params(strictness = "tight"), 3L, 80, "default")))
  expect_false(identical(
    key,
    ocular:::.rsCacheKey(117.8, -32, "2021-01-01", "2021-02-01",
                          "../../NDVI", "auto", 100, 100, 50, p,
                          3L, 90, "default")))
  expect_false(identical(
    key,
    ocular:::.rsCacheKey(117.8, -32, "2021-01-01", "2021-02-01",
                          "../../NDVI", "auto", 100, 100, 50, p,
                          3L, 80, "other-pipeline")))
  expect_false(identical(
    ocular:::.rsCacheKey(117.8, -32, "2021-01-01", "2021-02-01",
                          "NDVI", "auto", 100, 100, 50, p,
                          2L, 80, "default"),
    ocular:::.rsCacheKey(117.8, -32, "2021-01-01", "2021-02-01",
                          "NDVI", "auto", 100, 100, 50, p,
                          3L, 80, "default")))
})

test_that("corrupt or malformed batch cache entries are clean misses", {
  dir <- .test_artifact_dir("batch-cache-read")
  corrupt <- file.path(dir, "corrupt.rds")
  writeLines("not an rds", corrupt)
  expect_null(ocular:::.rsCacheRead(corrupt))
  malformed <- file.path(dir, "malformed.rds")
  saveRDS(list(data = data.frame(x = 1)), malformed)
  expect_null(ocular:::.rsCacheRead(malformed))
})

test_that("cache ownership guards reject broad and unmarked directories", {
  expect_error(clear_rs_cache("/"), "refusing")
  expect_error(clear_rs_cache(getwd()), "refusing")
  unmarked <- .test_artifact_dir("unmarked-cache")
  expect_error(clear_rs_cache(unmarked), "unmarked custom directory")

  marked <- .test_artifact_dir("marked-cache")
  marked <- ocular:::.initialiseCacheDir(
    marked, "rs", ocular:::.RS_CACHE_VERSION, "test")
  saveRDS(list(), file.path(marked, paste0(strrep("a", 24), ".rds")))
  writeLines("unrelated", file.path(marked, "notes.txt"))
  expect_identical(
    basename(ocular:::.cacheEntryFiles(marked, "^[0-9a-f]{24}[.]rds$")),
    paste0(strrep("a", 24), ".rds"))
})

test_that("cache clearing is selective in an owned directory", {
  marked <- .test_artifact_dir("selective-cache")
  marked <- ocular:::.initialiseCacheDir(
    marked, "rs", ocular:::.RS_CACHE_VERSION, "test")
  marker <- ocular:::.cacheMarker(
    marked, "rs", ocular:::.RS_CACHE_VERSION)
  entry <- file.path(marked, paste0(strrep("b", 24), ".rds"))
  unrelated <- file.path(marked, "notes.txt")
  saveRDS(list(), entry)
  writeLines("keep", unrelated)

  expect_true(ocular:::.cacheMarkerValid(
    marked, "rs", ocular:::.RS_CACHE_VERSION))
  expect_true(clear_rs_cache(marked))
  expect_false(file.exists(entry))
  expect_true(file.exists(unrelated))
})

test_that("default batch reduction uses the union for multiple areas", {
  testthat::local_mocked_bindings(
    boundary_delineation = identity,
    as_time_series = function(...)
      data.frame(date = as.Date("2021-01-01"), area_1 = 0.4,
                 area_2 = 0.6, union = 0.5),
    .package = "ocular")
  out <- ocular:::.defaultBatchPipeline(.make_test_rs())
  expect_equal(out$value, out$union)
})

test_that("batch validation does not touch the network", {
  sites <- data.frame(id = "a", longitude = 117.8, latitude = -32,
                      start_date = "2021-01-01", end_date = "2021-02-01")
  expect_error(batch_rs(sites, retries = 1.5), "positive integer")
  expect_error(batch_rs(sites, cache = NA), "cache")
  expect_error(batch_rs(sites, quiet = NA), "quiet")
  expect_error(batch_rs(sites, cloud_relax = Inf), "cloud_relax")
  expect_error(batch_rs(sites, pipeline = 1), "pipeline")

  bad <- sites
  bad$longitude <- Inf
  out <- batch_rs(bad, cache = FALSE, quiet = TRUE)
  expect_equal(nrow(out), 0L)
  expect_equal(nrow(attr(out, "failures")), 1L)
  expect_error(batch_rs(bad, cache = FALSE, quiet = TRUE, on_error = "stop"),
               "longitude")
})

test_that("batch parameter validation preserves explicit ownership", {
  sites <- data.frame(id = "a", longitude = 117.8, latitude = -32,
                      start_date = "2021-01-01", end_date = "2021-02-01")
  seen_user_set <- NULL
  fake_get <- function(..., params) {
    seen_user_set <<- attr(params, "user_set", exact = TRUE)
    rs <- ocular:::.newOcular()
    rs$geom$source <- "landsat-8"
    rs$scenes <- list(list())
    rs
  }
  pipeline <- function(rs)
    data.frame(date = as.Date("2021-01-05"), value = 0.5)
  testthat::local_mocked_bindings(get_rs = fake_get, .package = "ocular")

  batch_rs(sites, params = rs_params(area_sensitivity = 0.2),
           pipeline = pipeline, cache = FALSE, retries = 1L, quiet = TRUE)
  expect_identical(seen_user_set, "area_sensitivity")
})

test_that("batch retries fetching but runs a deterministic pipeline once", {
  fetch_count <- 0L
  pipeline_count <- 0L
  clouds <- numeric()
  fake_get <- function(..., max_cloud_cover) {
    fetch_count <<- fetch_count + 1L
    clouds <<- c(clouds, max_cloud_cover)
    if (fetch_count < 3L) stop("transient")
    rs <- ocular:::.newOcular()
    rs$geom$source <- "landsat-8"
    rs$scenes <- list(list())
    rs
  }
  pipeline <- function(rs) {
    pipeline_count <<- pipeline_count + 1L
    data.frame(date = as.Date("2021-01-05"), value = 0.5)
  }
  testthat::local_mocked_bindings(get_rs = fake_get, .package = "ocular")
  got <- ocular:::.batchFetchOne(
    117.8, -32, "2021-01-01", "2021-02-01", "EVI2", "landsat-8",
    100, 100, 50, rs_params(), pipeline, 3L, 80)
  expect_false(inherits(got, "error"))
  expect_equal(fetch_count, 3L)
  expect_equal(pipeline_count, 1L)
  expect_equal(clouds, c(50, 50, 80))
})

test_that("batch accepts empty results and valid month buckets safely", {
  fake_get <- function(...) {
    rs <- ocular:::.newOcular()
    rs$geom$source <- "landsat-8"
    rs$scenes <- list(list())
    rs
  }
  testthat::local_mocked_bindings(get_rs = fake_get, .package = "ocular")
  empty <- ocular:::.batchFetchOne(
    117.8, -32, "2021-01-01", "2021-02-01", "EVI2", "landsat-8",
    100, 100, 50, rs_params(),
    function(rs) data.frame(date = character(), value = numeric()), 1L, 80)
  expect_equal(nrow(empty$data), 0L)

  monthly <- ocular:::.batchFetchOne(
    117.8, -32, "2021-01-01", "2021-02-01", "EVI2", "landsat-8",
    100, 100, 50, rs_params(),
    function(rs) data.frame(date = "2021-01", value = 0.5), 1L, 80)
  expect_identical(monthly$data$date, "2021-01")

  matrix_value <- ocular:::.batchFetchOne(
    117.8, -32, "2021-01-01", "2021-02-01", "EVI2", "landsat-8",
    100, 100, 50, rs_params(),
    function(rs) data.frame(date = "2021-01-05",
                            value = I(matrix(c(1, 2), nrow = 1L))), 1L, 80)
  expect_s3_class(matrix_value, "error")
})

test_that("one retry never applies cloud relaxation", {
  clouds <- numeric()
  testthat::local_mocked_bindings(
    get_rs = function(..., max_cloud_cover){
      clouds <<- c(clouds, max_cloud_cover)
      stop("transient")
    }, .package = "ocular")
  got <- ocular:::.batchFetchOne(
    117.8, -32, "2021-01-01", "2021-02-01", "EVI2", "landsat-8",
    100, 100, 50, rs_params(), identity, 1L, 80)
  expect_s3_class(got, "error")
  expect_equal(clouds, 50)
})

test_that("batch cache hits avoid a second fetch", {
  cache_dir <- .test_artifact_dir("batch-cache-hit")
  fetch_count <- 0L
  fake_get <- function(...) {
    fetch_count <<- fetch_count + 1L
    rs <- ocular:::.newOcular()
    rs$geom$source <- "landsat-8"
    rs$scenes <- list(list())
    rs
  }
  pipeline <- function(rs)
    data.frame(date = as.Date("2021-01-05"), value = 0.5)
  testthat::local_mocked_bindings(get_rs = fake_get, .package = "ocular")
  sites <- data.frame(id = "a", longitude = 117.8, latitude = -32,
                      start_date = "2021-01-01", end_date = "2021-02-01")
  first <- batch_rs(sites, source = "landsat-8", pipeline = pipeline,
                    pipeline_cache_key = "test-v1", cache_dir = cache_dir,
                    retries = 1L, quiet = TRUE)
  second <- batch_rs(sites, source = "landsat-8", pipeline = pipeline,
                     pipeline_cache_key = "test-v1", cache_dir = cache_dir,
                     retries = 1L, quiet = TRUE)
  expect_equal(fetch_count, 1L)
  expect_equal(first, second)
})
