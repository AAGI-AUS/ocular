.split_test_rectangle <- function(n, rows, columns) {
  out <- matrix(FALSE, n, n)
  out[rows, columns] <- TRUE
  out
}

.split_test_fixture <- function(kind, n = 40L) {
  empty <- matrix(FALSE, n, n)

  if (identical(kind, "balanced")) {
    target <- .split_test_rectangle(n, 13:27, 10:19)
    neighbour <- .split_test_rectangle(n, 13:27, 20:30)
    point <- c(20L, 15L)
  } else if (identical(kind, "minority_3x3")) {
    joined <- .split_test_rectangle(n, 13:27, 10:30)
    neighbour <- .split_test_rectangle(n, 19:21, 19:21)
    target <- joined & !neighbour
    point <- c(20L, 15L)
  } else if (identical(kind, "dominant_4_17")) {
    target <- .split_test_rectangle(n, 13:27, 10:13)
    neighbour <- .split_test_rectangle(n, 13:27, 14:30)
    point <- c(20L, 12L)
  } else if (identical(kind, "diagonal")) {
    target <- .split_test_rectangle(n, 10:15, 10:15)
    neighbour <- .split_test_rectangle(n, 16:21, 16:21)
    point <- c(12L, 12L)
  } else if (identical(kind, "irregular_uniform")) {
    target <- empty
    target[10:22, 10:12] <- TRUE
    target[18:22, 10:24] <- TRUE
    target[14:22, 22:24] <- TRUE
    neighbour <- empty
    point <- c(19L, 11L)
  } else {
    stop("Unknown split test fixture: ", kind, call. = FALSE)
  }

  alive <- target | neighbour
  feature <- matrix(0.1, n, n)
  feature[target] <- 0.7
  feature[neighbour] <- 0.2

  list(
    name = kind,
    n = n,
    point = point,
    target = target,
    neighbour = neighbour,
    alive = alive,
    feature = feature,
    mu = 0.7
  )
}

.split_test_transform <- function(fixture, transpose = FALSE,
                                  reflect = FALSE, label = "identity") {
  transform_matrix <- function(x) {
    if (isTRUE(transpose)) x <- t(x)
    if (isTRUE(reflect)) {
      if (isTRUE(transpose)) {
        x <- x[rev(seq_len(nrow(x))), , drop = FALSE]
      } else {
        x <- x[, rev(seq_len(ncol(x))), drop = FALSE]
      }
    }
    x
  }

  point_marker <- matrix(FALSE, fixture$n, fixture$n)
  point_marker[fixture$point[1L], fixture$point[2L]] <- TRUE
  point_marker <- transform_matrix(point_marker)
  point <- which(point_marker, arr.ind = TRUE)[1L, ]

  fixture$name <- paste(fixture$name, label, sep = "_")
  fixture$point <- as.integer(point)
  fixture$target <- transform_matrix(fixture$target)
  fixture$neighbour <- transform_matrix(fixture$neighbour)
  fixture$alive <- transform_matrix(fixture$alive)
  fixture$feature <- transform_matrix(fixture$feature)
  fixture
}

.split_test_variants <- function(fixture) {
  specs <- list(
    horizontal_point_first = list(transpose = FALSE, reflect = FALSE),
    horizontal_point_last = list(transpose = FALSE, reflect = TRUE),
    vertical_point_first = list(transpose = TRUE, reflect = FALSE),
    vertical_point_last = list(transpose = TRUE, reflect = TRUE)
  )
  out <- vector("list", length(specs))
  names(out) <- names(specs)
  for (i in seq_along(specs)) {
    out[[i]] <- .split_test_transform(
      fixture,
      transpose = specs[[i]]$transpose,
      reflect = specs[[i]]$reflect,
      label = names(specs)[i]
    )
  }
  out
}

.split_test_rs <- function(fixture, multiple = FALSE, linear = FALSE,
                           sensitivity = 0.05, gate = NULL) {
  rs <- ocular:::.newOcular()
  rs$spec$index_name <- "EVI2"
  rs$geom$source <- "synthetic"
  rs$geom$pixel_size_m <- 10
  rs$geom$nr <- fixture$n
  rs$geom$nc <- fixture$n
  rs$geom$centre_rc <- fixture$point
  rs$internals$feat_array <- array(
    fixture$feature, dim = c(fixture$n, fixture$n, 1L)
  )
  rs$internals$mu <- fixture$mu
  rs$state$alive_mat <- fixture$alive
  rs$state$decided_mat <- matrix(FALSE, fixture$n, fixture$n)
  args <- list(
    split_sensitivity = sensitivity,
    split_linear = linear,
    multiple_areas = multiple
  )
  if (!is.null(gate)) args$split_gate <- gate
  rs$params <- do.call(rs_params, args)
  ocular:::.validateOcular(rs)
  rs
}

.expect_safe_split <- function(before, after, fixture, require_cut = FALSE,
                               info = fixture$name) {
  expect_true(is.matrix(after), info = info)
  expect_true(is.logical(after), info = info)
  expect_identical(dim(after), dim(before), info = info)
  expect_false(anyNA(after), info = info)
  expect_false(any(after & !before), info = info)
  expect_true(after[fixture$point[1L], fixture$point[2L]], info = info)
  expect_true(all(after[fixture$target]), info = info)
  expect_true(any(after), info = info)
  if (isTRUE(require_cut))
    expect_true(any(fixture$neighbour & !after), info = info)
  invisible(after)
}

test_that("protected split component is finite, inclusive and 8-connected", {
  alive <- matrix(FALSE, 5L, 5L)
  alive[cbind(c(2L, 3L, 4L, 1L), c(2L, 3L, 4L, 5L))] <- TRUE
  features <- array(0, dim = c(5L, 5L, 2L))
  features[2L, 2L, ] <- c(0.75, 1.25)
  features[3L, 3L, ] <- c(0.5, 1.5)
  features[4L, 4L, ] <- c(0.25, 0.75)
  features[1L, 5L, ] <- c(0.75, 1.25)

  protected <- ocular:::.splitProtectedComponent(
    alive, features, mu_seed = 1, sens = 0.5, centre_rc = c(2L, 2L)
  )
  expected <- matrix(FALSE, 5L, 5L)
  expected[cbind(2:4, 2:4)] <- TRUE
  expect_identical(protected, expected)

  exact <- ocular:::.splitProtectedComponent(
    alive, features, mu_seed = 1, sens = 0, centre_rc = c(2L, 2L)
  )
  expected_zero <- matrix(FALSE, 5L, 5L)
  expected_zero[cbind(2:3, 2:3)] <- TRUE
  expect_identical(exact, expected_zero)

  dead <- alive
  dead[2L, 2L] <- FALSE
  expect_null(ocular:::.splitProtectedComponent(
    dead, features, 1, 0.5, c(2L, 2L)
  ))

  nonfinite <- features
  nonfinite[2L, 2L, ] <- NA_real_
  expect_null(ocular:::.splitProtectedComponent(
    alive, nonfinite, 1, 0.5, c(2L, 2L)
  ))

  nonmatching <- features
  nonmatching[2L, 2L, ] <- 0
  expect_null(ocular:::.splitProtectedComponent(
    alive, nonmatching, 1, 0.5, c(2L, 2L)
  ))
  expect_null(ocular:::.splitProtectedComponent(
    alive, features, Inf, 0.5, c(2L, 2L)
  ))
  expect_null(ocular:::.splitProtectedComponent(
    alive, array(1, dim = c(5L, 5L, 1L, 1L)),
    1, 0.5, c(2L, 2L)
  ))
})

test_that("split gate requires consecutive measurable qualifying profiles", {
  profiles <- list(
    list(median_alive = 0.5, los_rc = NULL, skip = FALSE),
    list(median_alive = 1.0, los_rc = NULL, skip = FALSE),
    list(median_alive = 0.4, los_rc = NULL, skip = FALSE),
    list(median_alive = NA_real_, los_rc = NULL, skip = TRUE),
    list(median_alive = 0.4, los_rc = NULL, skip = FALSE),
    list(median_alive = 0.4, los_rc = NULL, skip = FALSE)
  )
  profile_index <- 0L
  cut <- matrix(c(2L, 2L), nrow = 1L,
                dimnames = list(NULL, c("row", "col")))
  testthat::local_mocked_bindings(
    .scanLOS = function(...) {
      profile_index <<- profile_index + 1L
      profiles[[profile_index]]
    },
    .linearCut = function(...) cut,
    .package = "ocular"
  )

  rule <- ocular:::.splitAreaRule()
  state <- list(
    envelopes = list(y = list(axis = "y"), x = NULL),
    mu_seed = 1,
    sens = 0.5,
    thresh = 2L,
    is_multi = FALSE,
    linear_mode = TRUE,
    nr = 8L,
    nc = 5L,
    feat_array = array(1, dim = c(8L, 5L, 1L)),
    by_axis = list(y = list(consecutive = 0L, flags = list()),
                   x = list(consecutive = 0L, flags = list())),
    log_msg = function(...) invisible()
  )
  walk <- list(trail = cbind(row = seq_len(7L), col = rep(3L, 7L)))
  alive <- matrix(TRUE, 8L, 5L)
  decided <- matrix(FALSE, 8L, 5L)
  before <- alive
  counts <- integer(6L)

  for (i in seq_along(counts)) {
    action <- rule$apply_cell(
      i + 1L, "y", walk, state, alive, decided,
      rs = NULL, params = NULL, walker_state = NULL
    )
    applied <- ocular:::.applyAction(action, alive, decided)
    alive <- applied$alive_mat
    decided <- applied$decided_mat
    if (!is.null(applied$state)) state <- applied$state
    counts[i] <- state$by_axis$y$consecutive
    if (i < length(counts)) expect_identical(alive, before)
  }

  expect_identical(counts, c(1L, 0L, 1L, 0L, 1L, 0L))
  expect_false(alive[2L, 2L])
  expect_true(decided[2L, 2L])
})

test_that("non-finite seed signatures leave the candidate mask unchanged", {
  fixture <- .split_test_fixture("balanced")
  for (mu in list(NA_real_, numeric(0L), c(NA_real_, NA_real_),
                  c(0.7, Inf))) {
    rs <- .split_test_rs(fixture)
    rs$internals$mu <- mu
    rs$state$alive_mat[2:3, 2:3] <- TRUE
    before_alive <- rs$state$alive_mat
    before_decided <- rs$state$decided_mat
    out <- split_area(rs)
    expect_identical(out$state$alive_mat, before_alive,
                     info = paste("mu", paste(mu, collapse = ",")))
    expect_identical(out$state$decided_mat, before_decided,
                     info = paste("mu", paste(mu, collapse = ",")))
  }
})

test_that("unavailable anchor evidence leaves the candidate mask unchanged", {
  fixture <- .split_test_fixture("balanced")

  states <- list()
  dead <- .split_test_rs(fixture)
  dead$state$alive_mat[fixture$point[1L], fixture$point[2L]] <- FALSE
  dead$state$alive_mat[2:3, 2:3] <- TRUE
  states$dead_point <- dead

  nonfinite <- .split_test_rs(fixture)
  nonfinite$internals$feat_array[
    fixture$point[1L], fixture$point[2L],
  ] <- NA_real_
  nonfinite$state$alive_mat[2:3, 2:3] <- TRUE
  states$nonfinite_point <- nonfinite

  nonmatching <- .split_test_rs(fixture)
  nonmatching$internals$feat_array[
    fixture$point[1L], fixture$point[2L],
  ] <- 0.2
  nonmatching$state$alive_mat[2:3, 2:3] <- TRUE
  states$nonmatching_point <- nonmatching

  out_of_bounds <- .split_test_rs(fixture)
  out_of_bounds$geom$centre_rc <- c(0L, fixture$point[2L])
  out_of_bounds$state$alive_mat[2:3, 2:3] <- TRUE
  states$out_of_bounds_point <- out_of_bounds

  malformed <- .split_test_rs(fixture)
  malformed$internals$feat_array <- array(
    0.7, dim = c(fixture$n, fixture$n, 1L, 1L)
  )
  malformed$state$alive_mat[2:3, 2:3] <- TRUE
  states$four_dimensional_feature <- malformed

  mismatched <- .split_test_rs(fixture)
  mismatched$internals$feat_array <- array(
    0.7, dim = c(fixture$n + 1L, fixture$n, 1L)
  )
  mismatched$state$alive_mat[2:3, 2:3] <- TRUE
  states$mismatched_feature <- mismatched

  for (label in names(states)) {
    rs <- states[[label]]
    before_alive <- rs$state$alive_mat
    before_decided <- rs$state$decided_mat
    out <- split_area(rs)
    expect_identical(out$state$alive_mat, before_alive, info = label)
    expect_identical(out$state$decided_mat, before_decided, info = label)
  }
})

test_that("overall LOS scan handles odd, even and sparse finite support", {
  scan <- function(values) {
    n <- length(values)
    envelope <- list(
      axis = "y",
      near_rc = matrix(c(1L, 1L), nrow = 1L,
                       dimnames = list(NULL, c("row", "col"))),
      far_rc = matrix(c(1L, n), nrow = 1L,
                      dimnames = list(NULL, c("row", "col")))
    )
    ocular:::.scanLOS(
      envelope, c(1L, 1L), "S", matrix(TRUE, 1L, n),
      array(values, dim = c(1L, n, 1L)), 1L, n
    )
  }

  odd <- scan(1:5)
  expect_false(odd$skip)
  expect_identical(nrow(odd$los_rc), 5L)
  expect_equal(odd$median_alive, 3)

  even <- scan(1:4)
  expect_false(even$skip)
  expect_identical(nrow(even$los_rc), 4L)
  expect_equal(even$median_alive, 2.5)

  one_finite <- scan(c(NA_real_, NA_real_, 7, Inf, NA_real_))
  expect_false(one_finite$skip)
  expect_equal(one_finite$median_alive, 7)

  none_finite <- scan(c(NA_real_, NaN, Inf, -Inf, NA_real_))
  expect_true(none_finite$skip)
  expect_true(is.na(none_finite$median_alive))
})

test_that("independent-axis rejection restores snapshots and rule state", {
  alive <- matrix(FALSE, 5L, 5L)
  alive[2:4, 2:4] <- TRUE
  decided <- matrix(FALSE, 5L, 5L)
  decided[4L, 4L] <- TRUE
  protected <- matrix(FALSE, 5L, 5L)
  protected[3L, 3L] <- TRUE

  rs <- ocular:::.newOcular()
  rs$geom$nr <- 5L
  rs$geom$nc <- 5L
  rs$geom$centre_rc <- c(3L, 3L)
  rs$state$alive_mat <- alive
  rs$state$decided_mat <- decided
  rs$params <- rs_params()

  final_marker <- NA_integer_
  rule <- list(
    axes = c("y", "x"),
    axis_mode = "independent",
    build_walk = function(...) list(
      trail = matrix(c(1L, 1L), nrow = 1L,
                     dimnames = list(NULL, c("row", "col")))
    ),
    setup = function(...) list(
      skip = FALSE, protected_component = protected, marker = 0L
    ),
    apply_cell = function(si, axis, walk, state, alive_mat, decided_mat,
                          rs, params, walker_state) {
      state$marker <- state$marker + if (identical(axis, "y")) 1L else 10L
      cell <- if (identical(axis, "y")) c(3L, 3L) else c(2L, 2L)
      alive_mat[cell[1L], cell[2L]] <- FALSE
      decided_mat[cell[1L], cell[2L]] <- TRUE
      list(action = "both", alive_mat = alive_mat,
           decided_mat = decided_mat, state = state)
    },
    accept_axis = ocular:::.splitAreaRule()$accept_axis,
    finalize = function(state, ...) {
      final_marker <<- state$marker
      NULL
    }
  )

  out <- ocular:::.runStagedWalk(rs, rule, rs$params)
  expect_true(out$state$alive_mat[3L, 3L])
  expect_false(out$state$decided_mat[3L, 3L])
  expect_false(out$state$alive_mat[2L, 2L])
  expect_true(out$state$decided_mat[2L, 2L])
  expect_true(out$state$decided_mat[4L, 4L])
  expect_identical(final_marker, 10L)

  apply_count <- 0L
  rejected_only <- rule
  rejected_only$axes <- "y"
  rejected_only$max_iter <- function(params) 2L
  rejected_only$apply_cell <- function(si, axis, walk, state, alive_mat,
                                       decided_mat, rs, params, walker_state) {
    apply_count <<- apply_count + 1L
    state$marker <- state$marker + 1L
    alive_mat[3L, 3L] <- FALSE
    decided_mat[3L, 3L] <- TRUE
    list(action = "both", alive_mat = alive_mat,
         decided_mat = decided_mat, state = state)
  }
  rejected_only$finalize <- function(...) NULL
  rejected <- ocular:::.runStagedWalk(rs, rejected_only, rs$params)
  expect_identical(apply_count, 1L)
  expect_identical(rejected$state$alive_mat, alive)
  expect_identical(rejected$state$decided_mat, decided)

  no_protection <- rejected_only
  no_protection$max_iter <- function(params) 1L
  no_protection$setup <- function(...) list(
    skip = FALSE, protected_component = NULL, marker = 0L
  )
  no_protection$apply_cell <- function(si, axis, walk, state, alive_mat,
                                       decided_mat, rs, params, walker_state) {
    alive_mat[2L, 2L] <- FALSE
    decided_mat[2L, 2L] <- TRUE
    list(action = "both", alive_mat = alive_mat,
         decided_mat = decided_mat, state = state)
  }
  unavailable <- ocular:::.runStagedWalk(rs, no_protection, rs$params)
  expect_identical(unavailable$state$alive_mat, alive)
  expect_identical(unavailable$state$decided_mat, decided)
})

test_that("balanced joined lobes retain the prompted target and still split", {
  fixtures <- .split_test_variants(.split_test_fixture("balanced"))
  for (fixture in fixtures) {
    for (multiple in c(FALSE, TRUE)) {
      for (linear in c(FALSE, TRUE)) {
        label <- paste(fixture$name, "multiple", multiple, "linear", linear)
        rs <- .split_test_rs(fixture, multiple = multiple, linear = linear)
        out <- split_area(rs)
        .expect_safe_split(
          rs$state$alive_mat, out$state$alive_mat, fixture,
          require_cut = TRUE, info = label
        )
      }
    }
  }
})

test_that("dominant-neighbour 4:17 cases retain the prompted target", {
  fixtures <- .split_test_variants(.split_test_fixture("dominant_4_17"))
  for (fixture in fixtures) {
    for (multiple in c(FALSE, TRUE)) {
      for (linear in c(FALSE, TRUE)) {
        label <- paste(fixture$name, "multiple", multiple, "linear", linear)
        rs <- .split_test_rs(fixture, multiple = multiple, linear = linear)
        out <- split_area(rs)
        .expect_safe_split(
          rs$state$alive_mat, out$state$alive_mat, fixture,
          info = label
        )
      }
    }
  }
})

test_that("minority and irregular controls remain non-destructive", {
  minority <- .split_test_fixture("minority_3x3")
  uniform <- .split_test_fixture("irregular_uniform")
  diagonal <- .split_test_fixture("diagonal")

  for (multiple in c(FALSE, TRUE)) {
    for (linear in c(FALSE, TRUE)) {
      minority_rs <- .split_test_rs(
        minority, multiple = multiple, linear = linear
      )
      minority_out <- split_area(minority_rs)
      expect_identical(
        minority_out$state$alive_mat, minority_rs$state$alive_mat,
        info = paste("minority multiple", multiple, "linear", linear)
      )

      uniform_rs <- .split_test_rs(
        uniform, multiple = multiple, linear = linear
      )
      uniform_out <- split_area(uniform_rs)
      expect_identical(
        uniform_out$state$alive_mat, uniform_rs$state$alive_mat,
        info = paste("uniform multiple", multiple, "linear", linear)
      )

      diagonal_rs <- .split_test_rs(
        diagonal, multiple = multiple, linear = linear
      )
      diagonal_out <- split_area(diagonal_rs)
      .expect_safe_split(
        diagonal_rs$state$alive_mat, diagonal_out$state$alive_mat,
        diagonal,
        info = paste("diagonal multiple", multiple, "linear", linear)
      )
    }
  }
})

test_that("repeated and default-workflow split calls do not erode the target", {
  fixtures <- .split_test_variants(.split_test_fixture("balanced"))
  for (fixture in fixtures) {
    out <- .split_test_rs(fixture)
    for (i in seq_len(5L)) {
      before <- out$state$alive_mat
      out <- split_area(out)
      .expect_safe_split(
        before, out$state$alive_mat, fixture,
        info = paste(fixture$name, "direct pass", i)
      )
    }
  }

  fixture <- .split_test_fixture("balanced")
  rs <- .split_test_rs(fixture)
  real_split <- split_area
  split_states <- list()
  pass <- function(x, ...) x
  capture_split <- function(x, ...) {
    out <- real_split(x, ...)
    split_states[[length(split_states) + 1L]] <<- out$state$alive_mat
    out
  }
  testthat::local_mocked_bindings(
    trace_perimeter = pass,
    segment_area = pass,
    segment_interior = pass,
    split_area = capture_split,
    .package = "ocular"
  )

  out <- boundary_delineation(rs)
  expect_identical(length(split_states), 5L)
  for (i in seq_along(split_states)) {
    .expect_safe_split(
      rs$state$alive_mat, split_states[[i]], fixture,
      info = paste("default workflow split", i)
    )
  }
  expect_true(out$state$alive_mat[fixture$point[1L], fixture$point[2L]])
  expect_true(all(out$state$alive_mat[fixture$target]))
})

test_that("split gate default is the guarded two-profile setting", {
  expect_identical(rs_params()$split_gate, 2L)
})
