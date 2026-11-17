test_that("empty detection windows identify their index, dates and available windows", {
  rs <- .make_test_rs()
  rs$internals$feature_stack <- NULL
  rs$spec$end_date <- "2021-07-01"
  rs$scenes[[1L]]$date <- as.Date("2021-01-05")
  rs$scenes[[2L]]$date <- as.Date("2021-06-20")
  expected <- paste0(
    "Cannot prepare EVI2 detection features: no scenes in window 2 ",
    "(2021-03-03 to 2021-05-01). 2 of 3 windows have scenes.")
  expect_error(ocular:::.fetchFeatureStack(rs), expected, fixed = TRUE)

  rs$spec$x_metres <- rs$spec$y_metres <- 5000
  expect_error(ocular:::.fetchFeatureStack(rs), expected, fixed = TRUE)
  rs$scenes <- rs$scenes[1L]
  expect_error(ocular:::.fetchFeatureStack(rs), paste0(
    "no scenes in window 2 (2021-03-03 to 2021-05-01), ",
    "window 3 (2021-05-02 to 2021-07-01). 1 of 3 windows have scenes."),
    fixed = TRUE)
  rs$scenes <- list()
  expect_error(ocular:::.fetchFeatureStack(rs), paste0(
    "Cannot prepare EVI2 detection features: no scenes available from ",
    "2021-01-01 to 2021-07-01 (0 of 3 windows have scenes)."), fixed = TRUE)
})

test_that("reported detection dates follow the existing window boundary rule", {
  rs <- .make_test_rs()
  rs$internals$feature_stack <- NULL
  rs$spec$end_date <- "2021-06-30"
  rs$scenes[[1L]]$date <- as.Date("2021-01-05")
  rs$scenes[[2L]]$date <- as.Date("2021-03-02")
  expect_error(ocular:::.fetchFeatureStack(rs), paste0(
    "no scenes in window 3 (2021-05-01 to 2021-06-30). ",
    "2 of 3 windows have scenes."), fixed = TRUE)

  rs$scenes[[3L]] <- rs$scenes[[2L]]
  rs$scenes[[3L]]$date <- as.Date("2021-06-30")
  result <- ocular:::.fetchFeatureStack(rs)
  expect_identical(names(result$internals$feature_stack), c("w1", "w2", "w3"))
  expect_identical(result$params, rs$params)
})

test_that("window selections report the supplied values and available count", {
  rs <- .make_test_rs()
  expect_error(segment_area(rs, search_windows = 3L), paste0(
    "segment_area(): search_windows out of range: 3. ",
    "Available windows: 1 to 2 (2 total)."), fixed = TRUE)
  expect_error(segment_area(rs, search_windows = c(1L, 4L)), paste0(
    "segment_area(): search_windows out of range: 1, 4. ",
    "Available windows: 1 to 2 (2 total)."), fixed = TRUE)
  expect_error(segment_area(rs, search_windows = c(3L, 3L)),
               'search_windows must be "all" or a unique vector of positive integers.',
               fixed = TRUE)
})

test_that("stage parameter errors identify the public caller and supplied setting", {
  rs <- .make_test_rs()
  callers <- c("segment_area", "trace_perimeter", "split_area", "segment_interior",
               "as_raster", "as_time_series")
  for( caller in callers ){
    stage <- get(caller, envir = asNamespace("ocular"))
    expect_error(stage(rs, params = 1), paste0(
      caller, "(): params must be a uniquely named list."),
      fixed = TRUE)
    expect_error(stage(rs, params = list(search_widnows = 1)),
                 paste0(caller, "(): unknown params: search_widnows."), fixed = TRUE)
    expect_error(stage(rs, search_widnows = 1),
                 paste0(caller, "(): unknown params: search_widnows."), fixed = TRUE)
  }
})

test_that("plain parameter lists and genuine ownership metadata are distinguished", {
  typo <- list(search_widnows = 1)
  expect_error(get_rs(NA, params = typo), paste0(
    "get_rs(): params contains unknown parameter(s): search_widnows."), fixed = TRUE)
  sites <- data.frame(id = "test-site", longitude = 117.8, latitude = -32,
                      start_date = "2021-01-01", end_date = "2021-02-01")
  expect_error(batch_rs(sites, params = typo), paste0(
    "batch_rs(): params contains unknown parameter(s): search_widnows."), fixed = TRUE)

  metadata <- list(search_windows = "all")
  attr(metadata, "user_set") <- "search_widnows"
  expect_error(get_rs(NA, params = metadata), paste0(
    'get_rs(): attr(params, "user_set") contains unknown parameter(s): search_widnows.'),
    fixed = TRUE)
  expect_error(segment_area(.make_test_rs(), params = metadata), paste0(
    'segment_area(): attr(params, "user_set") contains unknown parameter(s): ',
    "search_widnows."), fixed = TRUE)
  rs <- .make_test_rs()
  attr(rs$params, "user_set") <- "search_widnows"
  expect_error(segment_area(rs), paste0(
    'segment_area(): attr(rs$params, "user_set") contains unknown parameter(s): ',
    "search_widnows."), fixed = TRUE)
})

test_that("diagnostic context preserves validation order and genuine user errors", {
  rs <- .make_test_rs()
  broken <- rs
  broken$params <- 1
  expect_error(segment_area(broken, params = 1),
               "segment_area(): rs$params must be a uniquely named list.", fixed = TRUE)
  expect_error(segment_area(rs, params = 1, search_windows = 3L),
               "segment_area(): params must be a uniquely named list", fixed = TRUE)
  expect_error(segment_area(rs, params = list(area_sensitivity = -1)),
               "area_sensitivity must", fixed = TRUE)
  genuine <- structure(list(message = "user expression failed", call = NULL, detail = 7),
                        class = c("user_error", "error", "condition"))
  caught <- tryCatch(segment_area(rs, params = stop(genuine)), error = identity)
  expect_identical(caught, genuine)

  original <- rs$params
  direct <- ocular:::.resolveParams(rs, list(trace_iter = 3), list())
  contextual <- ocular:::.resolveParams(rs, list(trace_iter = 3), list(),
                                        caller = "segment_area()")
  expect_identical(contextual, direct)
  expect_identical(rs$params, original)
})
