test_that("a pipeline uses its supplied object instead of a captured head", {
  rs <- .make_test_rs()
  other <- .make_test_rs(n = 25L)
  other$scenes <- lapply(other$scenes, function(scene){
    scene$index <- scene$index + 0.2
    scene
  })
  expected <- as_time_series(segment_area(other))
  pipeline <- build_pipeline(rs |> segment_area() |> as_time_series())
  alias <- rs
  aliased <- build_pipeline(alias |> segment_area() |> as_time_series())

  expect_identical(pipeline(other), expected)
  expect_identical(aliased(other), expected)
  expect_false(identical(expected, as_time_series(segment_area(rs))))
})

test_that("the rs placeholder is not evaluated when a pipeline is built", {
  env <- new.env(parent = environment())
  reads <- 0L
  makeActiveBinding("rs", function(){
    reads <<- reads + 1L
    stop("placeholder was evaluated")
  }, env)
  pipeline <- eval(quote(build_pipeline(rs |> as_time_series())), env)
  input <- .make_test_rs()

  expect_identical(pipeline(input), as_time_series(input))
  expect_identical(reads, 0L)
})

test_that("positional settings are captured after the pipeline input", {
  input <- .make_test_rs()
  bucket <- "monthly"
  named_setting <- build_pipeline(as_time_series(bucket))
  literal_setting <- build_pipeline(as_time_series("monthly"))
  outer_setting <- build_pipeline(segment_area() |> as_time_series("monthly"))
  bucket <- "yearly"

  expect_identical(named_setting(input), as_time_series(input, "monthly"))
  expect_identical(literal_setting(input), as_time_series(input, "monthly"))
  expect_identical(outer_setting(input),
                   as_time_series(segment_area(input), "monthly"))
  expect_identical(build_pipeline(segment_area(0.3))(input),
                   segment_area(input, 0.3))
})

test_that("an unbound first setting or input alias is rejected at build time", {
  env <- new.env(parent = asNamespace("ocular"))
  expect_error(eval(quote(build_pipeline(as_time_series(montly))), env),
               "'montly' has no value", fixed = TRUE)
  expect_error(eval(quote(build_pipeline(field |> segment_area())), env),
               "'field' has no value", fixed = TRUE)
  pipeline <- eval(quote(build_pipeline(rs |> as_time_series())), env)
  input <- .make_test_rs()
  expect_identical(pipeline(input), as_time_series(input))
})

test_that("pipelines reject standalone retrieval and captured rasters", {
  raster <- .make_test_rs()$scenes[[1L]]$index
  expect_error(build_pipeline(as_time_series(
    117.8, latitude = -32, start_date = "2021-01-01", end_date = "2021-02-01")),
    "retrieval arguments are not supported", fixed = TRUE)
  expect_error(build_pipeline(raster |> trace_perimeter()),
               "a raster cannot be the pipeline input", fixed = TRUE)
})
