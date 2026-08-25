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
               "Unknown params")
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
    "Unknown params")
})

test_that("geometry conversion and bounding boxes validate their inputs", {
  d <- metres_to_degrees(-32, distance = 1000, pixel_size = 30)
  expect_true(all(is.finite(d)))
  expect_error(metres_to_degrees(90), "latitude")

  bb <- point_to_bbox(117.8, -32, 1000, 600)
  expect_s3_class(bb, "bbox")
  expect_lt(bb[["xmin"]], 117.8)
  expect_gt(bb[["xmax"]], 117.8)
  expect_error(point_to_bbox(179.999, 0, 1000, 1000), "antimeridian")
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
