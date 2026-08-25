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
    "no field boundaries intersect the analysis area; ocular object unchanged",
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
    "no field boundary contains the supplied point; ocular object unchanged",
    fixed = TRUE)
  expect_identical(unchanged, rs)
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
