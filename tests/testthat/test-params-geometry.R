test_that("rs_params validates and normalizes every scalar control", {
  expect_identical(attr(rs_params(), "user_set", exact = TRUE),
                   character(0L))
  p <- rs_params(search_windows = c(3, 1), min_windows_alive = 2,
                 trace_iter = 4, split_gate = 2,
                 search_start_date = as.Date("2021-01-01"))
  expect_identical(p$search_windows, c(3L, 1L))
  expect_identical(p$min_windows_alive, 2L)
  expect_identical(p$trace_iter, 4L)
  expect_identical(p$split_gate, 2L)
  expect_identical(p$search_start_date, "2021-01-01")

  bad_calls <- list(
    quote(rs_params(area_sensitivity = NA_real_)),
    quote(rs_params(area_sensitivity = Inf)),
    quote(rs_params(local_alive_density = 2)),
    quote(rs_params(time_independent_windows = NA)),
    quote(rs_params(search_windows = c(1, 1))),
    quote(rs_params(search_windows = 1.5)),
    quote(rs_params(min_windows_alive = 1.5)),
    quote(rs_params(trace_iter = .Machine$integer.max + 1)),
    quote(rs_params(split_gate = Inf)),
    quote(rs_params(multiple_areas = NA)),
    quote(rs_params(perimeter_margins = c(0.8, 0.2))),
    quote(rs_params(search_start_date = "not-a-date")),
    quote(rs_params(search_start_date = "2021-02-01",
                    search_end_date = "2021-01-01")))
  for (call in bad_calls) expect_error(eval(call), info = deparse1(call))
})

test_that("dates that parse to a year below 1000 are rejected", {
  expect_error(rs_params(search_start_date = "01/05/2021"),
               "search_start_date must be a date in YYYY-MM-DD format.",
               fixed = TRUE)
  expect_error(rs_params(search_end_date = "31/12/2021"),
               "search_end_date must be a date in YYYY-MM-DD format.",
               fixed = TRUE)
  for (date in list(" 2021-05-01", "2021/05/01", "2021-5-1",
                    as.Date("2021-05-01"),
                    as.POSIXct("2021-05-01 10:00:00", tz = "UTC"))) {
    expect_identical(rs_params(search_start_date = date)$search_start_date,
                     "2021-05-01")
  }
})

test_that("rs_params records explicit defaults and NULL settings", {
  p <- rs_params(area_sensitivity = 0.2, interior_sensitivity = NULL)
  expect_setequal(attr(p, "user_set", exact = TRUE),
                  c("area_sensitivity", "interior_sensitivity"))
  expect_equal(p$area_sensitivity, 0.2)
  expect_null(p$interior_sensitivity)
})

test_that("parameter ownership distinguishes defaults from direct mutation", {
  rs <- .make_test_rs()
  expect_true(ocular:::.paramIsAtDefault(rs, "area_sensitivity"))

  attr(rs$params, "user_set") <- "area_sensitivity"
  expect_false(ocular:::.paramIsAtDefault(rs, "area_sensitivity"))

  attr(rs$params, "user_set") <- character(0L)
  rs$params$area_sensitivity <- 0.17
  expect_false(ocular:::.paramIsAtDefault(rs, "area_sensitivity"))

  attr(rs$params, "user_set") <- NULL
  expect_error(ocular:::.paramIsAtDefault(rs, "area_sensitivity"),
               "missing explicit parameter ownership metadata")

  rs$params <- list(split_gate = 1L)
  expect_error(ocular:::.paramIsAtDefault(rs, "area_sensitivity"),
               "missing explicit parameter ownership metadata")

  rs$params <- list()
  expect_true(ocular:::.paramIsAtDefault(rs, "area_sensitivity"))

  rs <- .make_test_rs()
  attr(rs$params, "user_set") <- c("area_sensitivity", "area_sensitivity")
  expect_error(ocular:::.paramIsAtDefault(rs, "area_sensitivity"),
               "must contain unique")
})

test_that("parameter resolution overlays values without mutating stored params", {
  rs <- .make_test_rs()
  before <- rs$params
  p <- ocular:::.resolveParams(rs, list(strictness = "tight"),
                               list(area_sensitivity = 0.12))
  expect_identical(p$strictness, "tight")
  expect_equal(p$area_sensitivity, 0.12)
  expect_identical(rs$params, before)
  expect_error(
    ocular:::.resolveParams(rs, list(search_index = "NDVI"), list()),
    "cannot be changed after get_rs")
  expect_error(
    ocular:::.resolveParams(rs, list(bogus = 1), list()),
               "unknown params")
})

test_that("a complete params list keeps stored fetch-time settings", {
  rs <- .make_test_rs()
  rs$params <- rs_params(search_index = "NDVI",
                         search_start_date = "2021-01-10")
  p <- ocular:::.resolveParams(rs, rs_params(trace_iter = 5L), list())
  expect_identical(p$search_index, "NDVI")
  expect_identical(p$search_start_date, "2021-01-10")
  expect_identical(p$trace_iter, 5L)
  expect_setequal(attr(p, "user_set", exact = TRUE),
                  c("search_index", "search_start_date", "trace_iter"))

  edited <- rs_params(trace_iter = 5L)
  edited$search_index <- "NBR"
  expect_error(ocular:::.resolveParams(rs, edited, list()),
               "cannot be changed after get_rs")
  expect_error(
    ocular:::.resolveParams(rs, rs_params(search_index = "EVI2"), list()),
    "cannot be changed after get_rs")
})

test_that("calibration baselines respect explicit default ownership", {
  rs <- .make_test_rs()
  rs$internals$calibration <- list(
    schema_version = ocular:::.CALIBRATION_SCHEMA_VERSION,
    baselines = list(area_sensitivity = 0.11))

  unowned <- ocular:::.resolveParams(rs, NULL, list())
  expect_equal(unowned$area_sensitivity, 0.11)

  attr(rs$params, "user_set") <- "area_sensitivity"
  owned <- ocular:::.resolveParams(rs, NULL, list())
  expect_equal(owned$area_sensitivity, 0.2)

  attr(rs$params, "user_set") <- character(0L)
  rs$params$area_sensitivity <- 0.17
  mutated <- ocular:::.resolveParams(rs, NULL, list())
  expect_equal(mutated$area_sensitivity, 0.17)

  attr(rs$params, "user_set") <- NULL
  expect_error(ocular:::.resolveParams(rs, NULL, list()),
               "missing explicit parameter ownership metadata")
})

test_that("resolution provenance contains only intentional ownership", {
  rs <- .make_test_rs()
  rs$internals$calibration <- list(
    schema_version = ocular:::.CALIBRATION_SCHEMA_VERSION,
    baselines = list(area_sensitivity = 0.11))

  unowned <- ocular:::.resolveParams(rs, NULL, list())
  expect_identical(attr(unowned, "user_set", exact = TRUE),
                   character(0L))

  attr(rs$params, "user_set") <- "area_sensitivity"
  supplied <- rs_params(strictness = "tight")
  combined <- ocular:::.resolveParams(
    rs, supplied,
    list(split_gate = 2L, interior_sensitivity = NULL))
  expect_setequal(attr(combined, "user_set", exact = TRUE),
                  c("area_sensitivity", "strictness", "split_gate"))

  missing_ownership <- .make_test_rs()
  attr(missing_ownership$params, "user_set") <- NULL
  expect_error(ocular:::.resolveParams(missing_ownership, NULL, list()),
               "missing explicit parameter ownership metadata")
})

test_that("parameter resolution preserves intentional NULL settings", {
  rs <- .make_test_rs()
  rs$params["interior_sensitivity"] <- list(NULL)

  stored <- ocular:::.resolveParams(rs, NULL, list())
  supplied <- ocular:::.resolveParams(
    rs, list(split_sensitivity = NULL), list())

  expect_null(stored$interior_sensitivity)
  expect_null(supplied$split_sensitivity)
  expect_true("split_sensitivity" %in%
                attr(supplied, "user_set", exact = TRUE))
  expect_error(
    ocular:::.resolveParams(rs, NULL, list(bogus = NULL)),
    "unknown params")
})

test_that("geometry conversion and bounding boxes validate their inputs", {
  d <- metres_to_degrees(-32, distance = 1000, pixel_size = 30)
  expect_true(all(is.finite(d)))
  expect_error(metres_to_degrees(90), "latitude")

  bb <- point_to_bbox(117.8, -32, 1000, 600)
  expect_s3_class(bb, "bbox")
  expect_lt(bb[["xmin"]], 117.8)
  expect_gt(bb[["xmax"]], 117.8)
  expect_error(point_to_bbox(179.999, 0, 1000, 1000),
               paste0("point_to_bbox: the requested area crosses the antimeridian. ",
                      "Split it into two WGS-84 bounding boxes."), fixed = TRUE)
})

test_that("point_to_bbox rejects boxes wider than 180 degrees", {
  half <- metres_to_degrees(-16.8, distance = 2000,
                            pixel_size = 1000)[["x_degrees"]]
  expect_error(point_to_bbox(180 - half * 1.0005, -16.8, 2000, 2000),
               paste0("point_to_bbox: the requested area crosses the antimeridian. ",
                      "Split it into two WGS-84 bounding boxes."), fixed = TRUE)
  expect_error(point_to_bbox(0, 84, 1400000, 1000), "antimeridian")
})

test_that("rs_utm handles ordinary and special UTM zones", {
  make_point <- function(lon, lat) {
    rs <- ocular:::.newOcular()
    rs$spec$point <- list(longitude = lon, latitude = lat)
    rs
  }
  expect_identical(rs_utm(make_point(117.8, -32)), "EPSG:32750")
  expect_identical(rs_utm(make_point(6, 60)), "EPSG:32632")
  expect_identical(rs_utm(make_point(15, 75)), "EPSG:32633")
  expect_identical(rs_utm(make_point(180, 0)), "EPSG:32660")
  expect_error(rs_utm(make_point(0, 85)), "outside")
})

test_that("a complete params list keeps stored and calibrated settings", {
  rs <- .make_test_rs()
  rs$params <- rs_params(strictness = "tight")
  rs$internals$calibration <- list(
    schema_version = ocular:::.CALIBRATION_SCHEMA_VERSION,
    baselines = list(area_sensitivity = 0.11))

  complete <- ocular:::.resolveParams(rs, rs_params(trace_iter = 5L), list())
  partial <- ocular:::.resolveParams(rs, list(trace_iter = 5L), list())
  expect_identical(complete$strictness, "tight")
  expect_equal(complete$area_sensitivity, 0.11)
  expect_identical(complete$trace_iter, 5L)
  expect_identical(complete, partial)

  edited <- rs_params()
  edited$area_sensitivity <- 0.3
  expect_equal(ocular:::.resolveParams(rs, edited, list())$area_sensitivity,
               0.3)
  chosen <- ocular:::.resolveParams(
    rs, rs_params(strictness = "balanced"), list())
  expect_identical(chosen$strictness, "balanced")
})
