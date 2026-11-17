test_that(".ftwQuantiseBbox snaps to the 0.01-degree grid", {
  qb <- ocular:::.ftwQuantiseBbox(
    c(xmin = 117.801, ymin = -32.004, xmax = 117.811, ymax = -31.996))
  expect_equal(unname(qb[["xmin"]]), 117.80, tolerance = 1e-9)
  expect_equal(unname(qb[["xmax"]]), 117.82, tolerance = 1e-9)  # ceiling
  expect_equal(unname(qb[["ymin"]]), -32.01, tolerance = 1e-9)  # floor
  expect_equal(unname(qb[["ymax"]]), -31.99, tolerance = 1e-9)
})

test_that(".ftwCacheKey is deterministic, 16 hex, source-sensitive", {
  qb <- c(xmin = 117.80, ymin = -32.01, xmax = 117.82, ymax = -31.99)
  k1 <- ocular:::.ftwCacheKey("dataset-A", qb)
  k2 <- ocular:::.ftwCacheKey("dataset-A", qb)
  expect_identical(k1, k2)
  expect_equal(nchar(k1), 16L)
  expect_false(identical(k1, ocular:::.ftwCacheKey("dataset-B", qb)))
  expect_false(identical(
    ocular:::.ftwCacheKey("dataset-A", qb, list(crs = 4326)),
    ocular:::.ftwCacheKey("dataset-A", qb, list(crs = 2193))))
  expect_false(identical(
    ocular:::.ftwCacheKey("dataset-A", qb, list(geometry_col = "geometry")),
    ocular:::.ftwCacheKey("dataset-A", qb, list(geometry_col = "geom"))))
})

test_that("FTW cache round-trips a payload, misses cleanly, and evicts LRU", {
  dir <- .test_artifact_dir("ftw-cache")

  payload <- sf::st_sf(
    a = 1:3,
    geometry = sf::st_sfc(sf::st_point(c(0, 0)), sf::st_point(c(1, 1)),
                          sf::st_point(c(2, 2)), crs = 4326))
  ocular:::.ftwCacheWrite("k1", payload, cache_dir = dir)
  got <- ocular:::.ftwCacheRead("k1", cache_dir = dir)
  expect_equal(got$polygons, payload)
  expect_null(ocular:::.ftwCacheRead("absent", cache_dir = dir))

  keys <- sprintf("%016x", seq_len(6L))
  for (key in keys) ocular:::.ftwCacheWrite(key, payload, cache_dir = dir)
  unrelated <- file.path(dir, "unrelated.rds")
  saveRDS("not an ocular cache key", unrelated)
  dropped <- ocular:::.ftwCacheEvictLRU(cache_dir = dir, max_entries = 3L)
  expect_gt(dropped, 0L)
  expect_lte(length(ocular:::.ftwCacheList(dir)), 3L)
  expect_true(file.exists(unrelated))

  atomic <- file.path(dir, "atomic.rds")
  saveRDS(1, atomic)
  expect_null(ocular:::.ftwCacheRead("atomic", cache_dir = dir))
  malformed <- file.path(dir, "malformed.rds")
  saveRDS(list(polygons = 1, meta = list()), malformed)
  expect_null(ocular:::.ftwCacheRead("malformed", cache_dir = dir))
  expect_silent(ocular:::.ftwCacheEvictLRU(cache_dir = dir,
                                            max_entries = 1L))
})

test_that("clear_ftw_cache removes only exact entries from its owned cache", {
  dir <- .test_artifact_dir("ftw-clear")
  dir <- ocular:::.initialiseCacheDir(
    dir, "ftw", ocular:::.FTW_CACHE_VERSION, "test")
  entry <- file.path(dir, paste0(strrep("c", 16), ".rds"))
  unrelated <- file.path(dir, "notes.txt")
  saveRDS(list(), entry)
  writeLines("keep", unrelated)

  testthat::local_mocked_bindings(
    .ftwCacheDir = function(create = TRUE) dir,
    .package = "ocular")
  expect_true(clear_ftw_cache())
  expect_false(file.exists(entry))
  expect_true(file.exists(unrelated))
})

test_that("FTW source labels redact credentials and full paths", {
  source <- "https://user:password@example.test/private/fields.parquet?token=secret"
  label <- ocular:::.safeFtwSourceLabel(source)
  expect_identical(label, "https://.../fields.parquet")
  expect_false(grepl("user|password|token|secret", label))
  message <- ocular:::.redactFtwCondition(simpleError(source), source)
  expect_false(grepl("user|password|token|secret", message))
})

test_that("FTW error redaction removes the truncated DuckDB query excerpt", {
  source <- "gs://svc-user:SECRETKEY@my-bucket/private/fields.parquet"
  duckdb_message <- paste0(
    "Missing Extension Error: File ", source,
    " requires the extension httpfs to be loaded\n\n",
    "LINE 2:        FROM read_parquet('gs://svc-user:SECRETKEY@my-bucket/private...\n",
    "                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n",
    "Context: rapi_prepare")
  message <- ocular:::.redactFtwCondition(simpleError(duckdb_message), source)
  expect_false(grepl("svc-user|SECRETKEY|my-bucket|LINE 2", message))
  expect_identical(
    message,
    paste0("Missing Extension Error: File gs://.../fields.parquet requires ",
           "the extension httpfs to be loaded\n\nContext: rapi_prepare"))

  plain <- "IO Error: No files found that match the pattern \"x.parquet\""
  expect_identical(
    ocular:::.redactFtwCondition(simpleError(plain), "x.parquet"), plain)
})

test_that("add_ftw_prior validates its inputs", {
  expect_error(add_ftw_prior(list(), source = "x.parquet"), "ocular object")
  rs <- .make_test_rs()
  expect_error(
    add_ftw_prior(rs, source = ""),
    "single path or URL to a field boundary GeoParquet",
    fixed = TRUE)
  expect_error(add_ftw_prior(rs, source = "x.parquet", refresh = NA),
               "refresh")
  expect_error(add_ftw_prior(rs, source = "x.parquet", crs = "not-a-crs"),
               "crs")
})

test_that("add_ftw_prior uses field boundary terminology in no-match messages", {
  rs <- .make_test_rs()

  testthat::local_mocked_bindings(
    .ftwCacheRead = function(key) NULL,
    .ftwQuery = function(...) NULL,
    .package = "ocular")
  expect_message(
    unchanged <- add_ftw_prior(rs, source = "fields.parquet"),
    "no field boundaries intersect the analysis area. The ocular object is unchanged",
    fixed = TRUE)
  expect_identical(unchanged, rs)

  testthat::local_mocked_bindings(
    .ftwQuery = function(...) data.frame(id = 1L),
    .ftwCacheWrite = function(...) invisible(NULL),
    .ftwCacheEvictLRU = function(...) invisible(NULL),
    .selectFtwField = function(polygons, rs) NULL,
    .package = "ocular")
  expect_message(
    unchanged <- add_ftw_prior(rs, source = "fields.parquet"),
    "no field boundary contains the supplied point. The ocular object is unchanged",
    fixed = TRUE)
  expect_identical(unchanged, rs)
})

test_that("add_ftw_prior keys the cache of a local file on the file itself", {
  rs <- .make_test_rs()
  fs <- rs$internals$feature_stack
  field <- sf::st_as_sf(terra::project(
    terra::as.polygons(terra::ext(fs), crs = terra::crs(fs)), "EPSG:4326"))
  dir <- .test_artifact_dir("ftw-local-key")
  cache <- file.path(dir, "cache")
  for( sub in c("a", "b", "cache") ){
    dir.create(file.path(dir, sub))
  }
  writeLines("2018", file.path(dir, "a", "fields.parquet"))
  writeLines("2023", file.path(dir, "b", "fields.parquet"))

  testthat::local_mocked_bindings(
    .ftwCacheDir = function(create = TRUE) cache,
    .ftwQuery = function(bbox_wgs84, source, ...){
      out <- field
      out$source_year <- as.integer(readLines(source))
      out
    },
    .package = "ocular")
  year_of <- function(source){
    suppressMessages(add_ftw_prior(rs, source = source))$geom$ftw_prior$source_year
  }

  old <- getwd()
  tryCatch({
    setwd(file.path(dir, "a"))
    expect_identical(year_of("fields.parquet"), 2018L)
    setwd(file.path(dir, "b"))
    expect_identical(year_of("fields.parquet"), 2023L)
  }, finally = setwd(old))

  path <- file.path(dir, "a", "fields.parquet")
  expect_identical(year_of(path), 2018L)
  writeLines("2024", path)
  Sys.setFileTime(path, Sys.time() + 60)
  expect_identical(year_of(path), 2024L)
})

test_that("FTW cache key of a URL or a glob is the source string", {
  qb <- c(xmin = 117.80, ymin = -32.01, xmax = 117.82, ymax = -31.99)
  for( source in c("https://example.test/fields.parquet",
                   "s3://bucket/fields/*.parquet",
                   "no-such-dir/*.parquet") ){
    expect_identical(ocular:::.ftwCacheKey(ocular:::.ftwSourceId(source), qb),
                     ocular:::.ftwCacheKey(source, qb))
  }
})

test_that("add_ftw_prior does not cache a result without any source_year", {
  rs <- .make_test_rs()
  fs <- rs$internals$feature_stack
  field <- sf::st_as_sf(terra::project(
    terra::as.polygons(terra::ext(fs), crs = terra::crs(fs)), "EPSG:4326"))
  cache <- .test_artifact_dir("ftw-no-year")
  year <- NA_integer_

  testthat::local_mocked_bindings(
    .ftwCacheDir = function(create = TRUE) cache,
    .ftwQuery = function(...){
      out <- field
      out$source_year <- year
      out
    },
    .package = "ocular")
  source <- "/data/fields.parquet"

  first <- suppressMessages(add_ftw_prior(rs, source = source))
  expect_null(first$geom$ftw_prior$source_year)
  expect_length(ocular:::.ftwCacheList(cache), 0L)
  year <- 2020L
  second <- suppressMessages(add_ftw_prior(rs, source = source))
  expect_identical(second$geom$ftw_prior$source_year, 2020L)
  expect_length(ocular:::.ftwCacheList(cache), 1L)

  year <- NA_integer_
  suppressMessages(add_ftw_prior(rs, source = source, datetime_col = NULL))
  expect_length(ocular:::.ftwCacheList(cache), 2L)
})

test_that("FTW query never installs DuckDB extensions", {
  body_text <- paste(deparse(body(ocular:::.ftwQuery)), collapse = "\n")
  expect_false(grepl('dbExecute\\(con, "INSTALL', body_text, fixed = FALSE))
  expect_match(body_text, "autoinstall_known_extensions", fixed = TRUE)
  expect_match(body_text, "autoload_known_extensions", fixed = TRUE)
})

test_that("FTW query errors direct users to the input schema", {
  body_text <- paste(deparse(body(ocular:::.ftwQuery)), collapse = "\n")
  expect_false(grepl("assumptions to check", body_text, fixed = TRUE))
  expect_match(body_text, "source schema", fixed = TRUE)
})
