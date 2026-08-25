# =========================================================================
# Core parameters -- full populated tuning list
# =========================================================================

#' Validate boundary-cleanup thresholds (internal)
#' @noRd
.validateBoundaryCleanup <- function(cb){

  if( is.null(cb) ) return(invisible(TRUE))
  if( !is.numeric(cb) || length(cb) < 2L || length(cb) > 3L ||
      anyNA(cb) || any(!is.finite(cb)) )
    stop("boundary cleanup must be a numeric vector of length 2 or 3.", call. = FALSE)
  if( cb[1] < 0 || cb[1] > 1 || cb[2] < 0 || cb[2] > 1 )
    stop("boundary cleanup low and high must be in [0, 1].", call. = FALSE)
  if( cb[1] > cb[2] )
    stop("boundary cleanup low must be <= high.", call. = FALSE)
  if( cb[1] == cb[2] && length(cb) != 3 )
    stop("boundary cleanup requires a 3rd element (bias 0 or 1) when low == high.", call. = FALSE)
  if( length(cb) == 3 && !cb[3] %in% c(0, 1) )
    stop("boundary cleanup bias must be 0 or 1.", call. = FALSE)
  return(invisible(TRUE))
}

#' Test for a finite numeric scalar (internal)
#' @noRd
.isScalarNumber <- function(x){
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

#' Test for an integer-valued numeric scalar (internal)
#' @noRd
.isWholeNumber <- function(x){
  .isScalarNumber(x) && x >= -.Machine$integer.max &&
    x <= .Machine$integer.max && x == trunc(x)
}

#' Validate and normalise an optional date (internal)
#' @noRd
.normaliseOptionalDate <- function(x, name){
  if( is.null(x) ) return(NULL)
  if( length(x) != 1L || is.na(x) )
    stop(name, " must be NULL or a single date.", call. = FALSE)
  value <- tryCatch(suppressWarnings(as.Date(x)),
                    error = function(e) as.Date(NA))
  if( is.na(value) )
    stop(name, " could not be parsed as a date.", call. = FALSE)
  as.character(value)
}

#' Read explicit parameter ownership metadata (internal)
#'
#' New parameter lists carry the explicitly supplied canonical names in a
#' user_set attribute. Plain named lists have no such metadata, so their
#' names are the explicit settings. This fallback keeps partial hand-written
#' lists useful while preserving exact ownership from rs_params().
#'
#' @noRd
.paramUserSet <- function(params, caller = "params"){
  user_set <- attr(params, "user_set", exact = TRUE)
  if( is.null(user_set) ) user_set <- names(params)
  if( is.null(user_set) ) user_set <- character(0L)
  if( !is.character(user_set) || anyNA(user_set) ||
      any(!nzchar(user_set)) || anyDuplicated(user_set) )
    stop(caller, ": attr(params, \"user_set\") must contain unique, ",
         "non-empty parameter names.", call. = FALSE)
  unknown <- setdiff(user_set, names(rs_params()))
  if( length(unknown) )
    stop(caller, ": attr(params, \"user_set\") contains unknown parameter(s): ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  user_set
}

#' Configure an ocular workflow
#'
#' Returns a complete named list of parameters for retrieval and delineation.
#' \code{get_rs(..., params = ...)} stores this list on the \code{ocular}
#' object and downstream stages use it by default. Most applications can begin
#' with the defaults; individual settings can then be changed for retrieval,
#' initial area growth, boundary refinement, or output requirements.
#'
#' A named argument supplied to a stage call (for example,
#' \code{trace_perimeter(rs, multiple_areas = TRUE)}) overrides the stored
#' setting for that call only and does not mutate \code{rs$params}.
#'
#' @section Parameter reference:
#'
#' \describe{
#'   \item{\strong{Calibration}}{Where support is adequate, calibration derives
#'     per-window signatures and matching tolerances for initial area growth.
#'     Without an attached field boundary polygon, calibration uses a complete
#'     seed-centred square and expands it until its precision condition is met
#'     or the available feature grid limits further expansion. Explicit user
#'     values always take precedence over a calibrated baseline.}
#'   \item{\code{area_threshold}}{Lower and upper vegetation index bounds.
#'     A pixel passes the threshold test only if its VI is within
#'     \code{[area_threshold[1], area_threshold[2]]} (inclusive). The bounds
#'     are used by \code{segment_area()} when learning the seed signature and
#'     when checking the completed fill for vegetation. Default
#'     \code{c(0.2, 0.7)}.}
#'   \item{\code{area_sensitivity}}{Maximum distance from the seed
#'     signature for a pixel to be considered self-consistent.
#'     Default \eqn{0.20}.}
#'   \item{\code{local_alive_density}}{Minimum density of predicate-
#'     passing 8-neighbours required to admit a candidate during fill.
#'     Default \eqn{0.20}.}
#'   \item{\code{time_independent_windows}}{When \code{TRUE} (default),
#'     each search window runs its own fill and outputs are unioned.
#'     When \code{FALSE}, a single joint fill spans all windows.}
#'   \item{\code{search_windows}}{\code{"all"} (default) or an integer
#'     vector selecting which composite windows to use for
#'     detection.}
#'   \item{\code{search_index}}{Vegetation index used internally for
#'     detection. Default \code{"EVI2"}. Independent of the
#'     \code{index_name} used for output.}
#'   \item{\code{min_windows_alive}}{In joint mode, minimum windows in
#'     which a pixel must match the signature. Default \eqn{1L}.}
#'   \item{\code{strictness}}{Sample-free ordinal selector for fixed
#'     robustness multiples \code{k} applied to measured per-window
#'     dispersion: \code{"loose"}, \code{"balanced"} (default), or
#'     \code{"tight"}, mapping respectively to 3.0, 2.2, and 1.5.
#'     Spectral-only by design. These engineering defaults require external
#'     scientific validation.}
#'   \item{\code{search_start_date}, \code{search_end_date}}{Optional
#'     restricted date range for detection (output uses the full
#'     range). Must lie within \code{[start_date, end_date]} and span
#'     at least 2 days; \code{get_rs()} errors otherwise. Both default to
#'     \code{NULL}.}
#'   \item{\code{perimeter_margins}}{Default cleanup low/high (and
#'     optional bias) thresholds for boundary cleanup. Used by
#'     \code{trace_perimeter()} when no inline \code{perimeter_margins}
#'     argument is supplied. Default \code{c(0.5, 0.9)}.}
#'   \item{\code{trace_iter}}{Maximum cleanup iterations. Default
#'     \eqn{20L}.}
#'   \item{\code{multiple_areas}}{\code{TRUE} reports multiple
#'     disconnected areas in the output. \code{FALSE} (default) keeps
#'     only the centre-most rank-1 area.}
#'   \item{\code{area_separation_strict}}{Advanced, experimental
#'     connectivity selector. \code{FALSE} (default) uses 8-neighbour
#'     connectivity; \code{TRUE} uses 4-neighbour connectivity;
#'     \code{c(FALSE, n)} re-evaluates the top \code{n} centre-most areas
#'     under 4-neighbour connectivity.}
#'   \item{\code{merge_separate_areas}}{When \code{TRUE}, attempts to
#'     bridge separable areas via dead-pixel revival before single-
#'     area trim. Default \code{FALSE}.}
#'   \item{\code{interior_sensitivity}}{Maximum permitted difference
#'     between a line-of-sight alive-pixel median and the seed signature during
#'     within-field segmentation. Default \eqn{0.10}; \code{NULL} disables
#'     within-field segmentation.}
#'   \item{\code{interior_threshold}}{Lower and upper vegetation index bounds
#'     for the per-cell predicate used by \code{segment_interior()}.
#'     Alive interior cells outside
#'     \code{[interior_threshold[1], interior_threshold[2]]} are removed.
#'     During a divergent line-of-sight pass, in-range cells are retained or
#'     revived only when they also match the seed signature within
#'     \code{interior_sensitivity}. Default \code{c(0.1, 0.9)}.}
#'   \item{\code{split_sensitivity}}{Distance threshold for the
#'     area-split divergence and protection tests. A line-of-sight median
#'     qualifies as divergent when its difference from the seed signature is
#'     greater than or equal to this value. Finite alive cells whose scalar
#'     values differ by less than or equal to this value can form the protected
#'     8-connected component containing the supplied point. Exact equality
#'     therefore satisfies both conservative comparisons. Default \eqn{0.05};
#'     \code{NULL} disables splitting.}
#'   \item{\code{split_gate}}{Number of consecutive measurable divergent
#'     steps required to commit a split. A measurable
#'     non-divergent or unmeasurable profile clears the streak. Default
#'     \eqn{2L}.}
#'   \item{\code{split_linear}}{Use straight-line cuts instead of
#'     curved walks. Default \code{FALSE}.}
#' }
#'
#' @param area_threshold,area_sensitivity,local_alive_density,time_independent_windows See parameter reference.
#' @param search_windows,search_index,min_windows_alive See parameter reference.
#' @param search_start_date,search_end_date See parameter reference.
#' @param perimeter_margins,trace_iter See parameter reference.
#' @param multiple_areas,merge_separate_areas See parameter reference.
#' @param strictness See parameter reference.
#' @param interior_sensitivity,interior_threshold See parameter reference.
#' @param split_sensitivity,split_gate,split_linear See parameter reference.
#' @param area_separation_strict See parameter reference. This is an advanced,
#'   experimental control rather than a primary tuning parameter.
#' @param ... Catches unknown named arguments; errors if any are passed.
#' @returns A complete named list of all workflow parameters.
#' @examples
#' p <- rs_params(strictness = "tight", search_windows = c(1L, 3L))
#' p[c("strictness", "search_windows")]
#' @export
rs_params <- function(area_threshold                  = c(0.2, 0.7),
                      area_sensitivity                = 0.2,
                      local_alive_density           = 0.2,
                      time_independent_windows      = TRUE,
                      search_windows                = "all",
                      search_index                  = "EVI2",
                      min_windows_alive             = 1L,
                      search_start_date             = NULL,
                      search_end_date               = NULL,
                      perimeter_margins             = c(0.5, 0.9),
                      trace_iter                  = 20L,
                      multiple_areas               = FALSE,
                      merge_separate_areas         = FALSE,
                      ## Advanced connectivity rule for area enumeration:
                      ## FALSE = 8-conn; TRUE = 4-conn; c(0, n) = 8-conn
                      ## baseline with the top-n centre-most areas re-enumerated
                      ## under 4-conn.
                      area_separation_strict       = FALSE,
                      interior_sensitivity = 0.1,
                      interior_threshold   = c(0.1, 0.9),
                      split_sensitivity     = 0.05,
                      split_gate            = 2L,
                      split_linear          = FALSE,
                      strictness                      = "balanced",
                      ...){

  supplied <- names(as.list(match.call(expand.dots = FALSE)))[-1L]
  if( is.null(supplied) ) supplied <- character(0L)
  user_set <- setdiff(supplied, "...")

  ## Catch unknown named arguments.
  if( ...length() > 0L ){
    extra <- names(list(...))
    if( is.null(extra) || any(!nzchar(extra)) )
      stop("rs_params(): all arguments must be named.", call. = FALSE)
    stop("Unknown rs_params: ", paste(extra, collapse = ", "), call. = FALSE)
  }

  ## Validate the calibration strictness selector.
  if( !is.character(strictness) || length(strictness) != 1L ||
      !strictness %in% c("loose", "balanced", "tight") )
    stop("strictness must be one of \"loose\", \"balanced\", \"tight\".",
         call. = FALSE)

  ## Validate each area.
  if( !is.numeric(area_threshold) || length(area_threshold) != 2L ||
      anyNA(area_threshold) || any(!is.finite(area_threshold)) ||
      area_threshold[1L] > area_threshold[2L] )
    stop("area_threshold must be a length-2 numeric c(lower, upper) with lower <= upper.",
         call. = FALSE)
  if( !.isScalarNumber(area_sensitivity) || area_sensitivity < 0 )
    stop("area_sensitivity must be a single non-negative numeric.", call. = FALSE)
  if( !.isScalarNumber(local_alive_density) ||
      local_alive_density < 0 || local_alive_density > 1 )
    stop("local_alive_density must be a numeric in [0, 1].", call. = FALSE)
  if( !is.logical(time_independent_windows) ||
      length(time_independent_windows) != 1L ||
      is.na(time_independent_windows) )
    stop("time_independent_windows must be a single logical.", call. = FALSE)
  valid_windows <- identical(search_windows, "all") ||
    (is.numeric(search_windows) && length(search_windows) > 0L &&
       !anyNA(search_windows) && all(is.finite(search_windows)) &&
       all(search_windows >= 1) &&
       all(search_windows <= .Machine$integer.max) &&
       all(search_windows == trunc(search_windows)) &&
       !anyDuplicated(search_windows))
  if( !valid_windows )
    stop("search_windows must be \"all\" or a unique vector of positive integers.",
         call. = FALSE)
  if( !identical(search_windows, "all") )
    search_windows <- as.integer(search_windows)
  if( !is.character(search_index) || length(search_index) != 1L ||
      is.na(search_index) || !nzchar(search_index) )
    stop("search_index must be a single character string.", call. = FALSE)
  if( is.null(s2_index_list[[search_index]]) &&
      is.null(landsat_index_list[[search_index]]) )
    stop("search_index '", search_index, "' is not a known index. Available: ",
         paste(unique(c(names(s2_index_list), names(landsat_index_list))),
               collapse = ", "), call. = FALSE)
  if( !.isWholeNumber(min_windows_alive) || min_windows_alive < 1L )
    stop("min_windows_alive must be a positive integer.", call. = FALSE)
  min_windows_alive <- as.integer(min_windows_alive)
  search_start_date <- .normaliseOptionalDate(search_start_date,
                                               "search_start_date")
  search_end_date <- .normaliseOptionalDate(search_end_date,
                                             "search_end_date")
  if( !is.null(search_start_date) && !is.null(search_end_date) &&
      as.Date(search_start_date) > as.Date(search_end_date) )
    stop("search_start_date must be on or before search_end_date.", call. = FALSE)
  .validateBoundaryCleanup(perimeter_margins)
  if( !.isWholeNumber(trace_iter) || trace_iter < 1L )
    stop("trace_iter must be a positive integer.", call. = FALSE)
  trace_iter <- as.integer(trace_iter)
  if( !is.logical(multiple_areas) ||
      length(multiple_areas) != 1L || is.na(multiple_areas) )
    stop("multiple_areas must be TRUE or FALSE.", call. = FALSE)
  if( !is.logical(merge_separate_areas) ||
      length(merge_separate_areas) != 1L ||
      is.na(merge_separate_areas) )
    stop("merge_separate_areas must be TRUE or FALSE.", call. = FALSE)
  ok_fss <- FALSE
  if( is.logical(area_separation_strict) &&
      length(area_separation_strict) == 1L &&
      !is.na(area_separation_strict) ){
    ok_fss <- TRUE
  }else if( is.numeric(area_separation_strict) &&
            length(area_separation_strict) == 2L &&
            !anyNA(area_separation_strict) &&
            all(is.finite(area_separation_strict)) &&
            area_separation_strict[1] == 0 &&
            area_separation_strict[2] >= 1L &&
            area_separation_strict[2] <= .Machine$integer.max &&
            area_separation_strict[2] == trunc(area_separation_strict[2]) ){
    ok_fss <- TRUE
    area_separation_strict <- c(0, as.integer(area_separation_strict[2]))
  }
  if( !ok_fss )
    stop("area_separation_strict must be TRUE, FALSE, or c(FALSE, n) for positive integer n.", call. = FALSE)
  if( !is.null(interior_sensitivity) ){
    if( !.isScalarNumber(interior_sensitivity) ||
        interior_sensitivity < 0 || interior_sensitivity > 1 )
      stop("interior_sensitivity must be NULL or a single numeric in [0, 1].",
           call. = FALSE)
  }
  if( !is.numeric(interior_threshold) ||
      length(interior_threshold) != 2L ||
      anyNA(interior_threshold) || any(!is.finite(interior_threshold)) ||
      interior_threshold[1L] > interior_threshold[2L] )
    stop("interior_threshold must be a length-2 numeric c(lower, upper) with lower <= upper.",
         call. = FALSE)
  if( !is.null(split_sensitivity) ){
    if( !.isScalarNumber(split_sensitivity) ||
        split_sensitivity < 0 )
      stop("split_sensitivity must be NULL or a single non-negative numeric.",
           call. = FALSE)
  }
  if( !.isWholeNumber(split_gate) || split_gate < 1L )
    stop("split_gate must be a positive integer.", call. = FALSE)
  split_gate <- as.integer(split_gate)
  if( !is.logical(split_linear) || length(split_linear) != 1L ||
      is.na(split_linear) )
    stop("split_linear must be a single logical.", call. = FALSE)

  out <- list(
    area_threshold                  = area_threshold,
    area_sensitivity                = area_sensitivity,
    local_alive_density           = local_alive_density,
    time_independent_windows      = time_independent_windows,
    search_windows                = search_windows,
    search_index                  = search_index,
    min_windows_alive             = min_windows_alive,
    search_start_date             = search_start_date,
    search_end_date               = search_end_date,
    perimeter_margins             = perimeter_margins,
    trace_iter                  = trace_iter,
    multiple_areas               = multiple_areas,
    merge_separate_areas         = merge_separate_areas,
    area_separation_strict       = area_separation_strict,
    interior_sensitivity = interior_sensitivity,
    interior_threshold   = interior_threshold,
    split_sensitivity     = split_sensitivity,
    split_gate               = split_gate,
    split_linear           = split_linear,
    strictness                      = strictness
  )
  attr(out, "user_set") <- unique(user_set)
  return(out)
}
