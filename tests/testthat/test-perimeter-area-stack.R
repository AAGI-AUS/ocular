.d16_area_stack <- function(){

  rs <- .make_test_rs(n = 41L)
  alive <- matrix(FALSE, 41L, 41L)
  alive[10:30, 4:17] <- TRUE
  alive[10:30, 24:38] <- TRUE
  alive[9, 10] <- TRUE
  alive[31, 30] <- TRUE
  rs$state$alive_mat <- alive
  return(as_raster(rs, multiple_areas = TRUE))
}

.d16_expected <- function(x, multiple_areas){

  keep <- rep(FALSE, terra::ncell(x))
  for( name in c("area_1", "area_2") ){
    mask <- x[[name]]
    values <- terra::values(mask, mat = FALSE)
    alive <- !is.na(values) & values != 0
    if( !any(alive) ){
      next
    }
    terra::values(mask) <- as.integer(alive)
    cleaned <- trace_perimeter(mask, perimeter_margins = c(0.5, 0.9),
                               multiple_areas = multiple_areas)
    keep <- keep | !is.na(terra::values(cleaned, mat = FALSE))
  }
  values <- terra::values(x)
  values[!keep, ] <- NA_real_
  return(values)
}

test_that("raster cleanup unions cleaned areas before masking original values", {
  x <- .d16_area_stack()
  original <- terra::values(x)
  expect_equal(colSums(!is.na(original)), c(area_1 = 316, area_2 = 295, union = 611))

  for( multiple in c(FALSE, TRUE) ){
    out <- trace_perimeter(x, perimeter_margins = c(0.5, 0.9),
                           multiple_areas = multiple)
    values <- terra::values(out)
    expect_identical(values, .d16_expected(x, multiple))
    expect_equal(colSums(!is.na(values)), c(area_1 = 315, area_2 = 294, union = 609))
    expect_identical(!is.na(values[, "union"]),
                     !is.na(values[, "area_1"]) | !is.na(values[, "area_2"]))
    expect_identical(values[!is.na(values)], original[!is.na(values)])
    expect_identical(names(out), names(x))
    expect_true(terra::compareGeom(out, x, stopOnError = FALSE))
  }
})

test_that("area cleanup retains layer order, values and per-area settings", {
  x <- .d16_area_stack()
  values <- terra::values(x)
  values[!is.na(values[, 1L]), 1L] <- -seq_len(sum(!is.na(values[, 1L])))
  values[!is.na(values[, 2L]), 2L] <- seq_len(sum(!is.na(values[, 2L]))) / 1000
  values[!is.na(values[, 3L]), 3L] <- seq_len(sum(!is.na(values[, 3L])))
  values[20L * 41L + 30L, 3L] <- 0
  values[2L * 41L + 2L, 1L] <- 7
  values[2L * 41L + 2L, 3L] <- 8
  terra::values(x) <- values

  for( multiple in c(FALSE, TRUE) ){
    expected <- .d16_expected(x, multiple)
    out <- trace_perimeter(x, perimeter_margins = c(0.5, 0.9),
                           multiple_areas = multiple)
    expect_identical(terra::values(out), expected)
    order <- c(3L, 2L, 1L)
    reordered <- trace_perimeter(x[[order]], perimeter_margins = c(0.5, 0.9),
                                 multiple_areas = multiple)
    expect_identical(terra::values(reordered), expected[, order])
    expect_identical(names(reordered), names(x)[order])
  }
  one <- trace_perimeter(x, perimeter_margins = c(0.5, 0.9))
  many <- trace_perimeter(x, perimeter_margins = c(0.5, 0.9), multiple_areas = TRUE)
  expect_true(is.na(terra::values(one)[2L * 41L + 2L, 1L]))
  expect_identical(unname(terra::values(many)[2L * 41L + 2L, 1L]), 7)
})

test_that("the derived union does not supply an independent cleanup mask", {
  x <- .d16_area_stack()
  values <- terra::values(x)
  values[, "area_1"] <- NA_real_
  values[1L, "union"] <- 9
  terra::values(x) <- values
  out <- trace_perimeter(x, perimeter_margins = c(0.5, 0.9))
  expect_identical(terra::values(out), .d16_expected(x, FALSE))
  expect_equal(colSums(!is.na(terra::values(out))), c(area_1 = 0, area_2 = 294, union = 294))

  values[, "area_2"] <- NA_real_
  terra::values(x) <- values
  out <- trace_perimeter(x, perimeter_margins = c(0.5, 0.9))
  expect_true(all(is.na(terra::values(out))))
})

test_that("generic raster stacks retain their first-layer cleanup rule", {
  x <- .d16_area_stack()
  names(x) <- c("date_1", "date_2", "date_3")
  out <- trace_perimeter(x, perimeter_margins = c(0.5, 0.9))
  first <- trace_perimeter(x[[1L]], perimeter_margins = c(0.5, 0.9))
  keep <- !is.na(terra::values(first, mat = FALSE))
  expected <- terra::values(x)
  expected[!keep, ] <- NA_real_
  expect_identical(terra::values(out), expected)

  values <- terra::values(x)
  values[, 1L] <- NA_real_
  terra::values(x) <- values
  expect_identical(terra::values(trace_perimeter(x, perimeter_margins = c(0.5, 0.9))), values)
})
