# =========================================================================
# Segmentation: boundary delineation
# =========================================================================

#' Delineate a field boundary
#'
#' Combines area growth, perimeter cleanup, splitting and interior refinement
#' in a fixed sequence. \code{cleanup_boundary} sets the number of additional
#' split-and-cleanup pairs. The four stages can also be called individually
#' to build a different sequence.
#'
#' If the object already has a mask, the workflow processes that mask again
#' using the existing spectral data and current settings.
#'
#' The two closing cleanups use perimeter margins \code{c(0.3, 0.9)} and
#' \code{c(0.5, 0.5, 1)}. Supply \code{perimeter_margins} in \code{...} to
#' replace both values.
#'
#' @param x An \code{ocular} object in piped form, or longitude in standalone
#'   form.
#' @param vi_sensitivity Sets \code{area_sensitivity},
#'   \code{interior_sensitivity} and \code{split_sensitivity} together, unless
#'   an individual value is supplied in \code{...}. This shorthand belongs to
#'   this function, not \code{rs_params()}.
#' @param vi_threshold Sets \code{area_threshold} and
#'   \code{interior_threshold} together, unless an individual value is supplied
#'   in \code{...}. This shorthand belongs to this function, not \code{rs_params()}.
#' @param strictness Calibration control forwarded to every stage. It selects
#'   the robustness multiple applied to measured per-window dispersion. See
#'   \code{?rs_params}.
#' @param cleanup_boundary Integer in \code{[0, 5]} giving the number of
#'   additional split-and-cleanup pairs in the middle phase. Default \code{2L}.
#' @param params Optional list of settings for this call. See
#'   \code{\link{rs_params}} for the settings and override rules.
#' @param ... Named \code{rs_params()} overrides forwarded to every stage.
#'   A stage-specific value here takes precedence over the generic
#'   \code{vi_sensitivity}/\code{vi_threshold}. Standalone mode also forwards
#'   spatial and temporal retrieval arguments to \code{get_rs()}.
#' @returns An \code{ocular} object with its delineated field mask stored in
#'   \code{state$alive_mat}.
#' @export
boundary_delineation <- function(x,
                                 vi_sensitivity = NULL,
                                 vi_threshold = NULL,
                                 strictness = NULL,
                                 cleanup_boundary = 2L,
                                 params = NULL,
                                 ...){

  if( !.isWholeNumber(cleanup_boundary) ||
      cleanup_boundary < 0L || cleanup_boundary > 5L ){
    stop("boundary_delineation: cleanup_boundary must be a single integer in [0, 5].",
         call. = FALSE)
  }
  cleanup_boundary <- as.integer(cleanup_boundary)

  dots <- list(...)
  if( length(dots) > 0L &&
      (is.null(names(dots)) || anyNA(names(dots)) ||
       any(!nzchar(names(dots))) || anyDuplicated(names(dots))) ){
    stop("boundary_delineation: all arguments in ... must be uniquely named.",
         call. = FALSE)
  }
  param_names <- names(rs_params())

  if( is_rs(x) ){
    rs <- x
    unknown <- setdiff(names(dots), param_names)
    if( length(unknown) ){
      stop("boundary_delineation: unknown parameter(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    inline <- dots
  }else{
    # Pass retrieval arguments to get_rs() and keep stage overrides for processing.
    fetch_names <- setdiff(names(formals(get_rs)), c("longitude", "params"))
    known <- union(fetch_names, param_names)
    unknown <- setdiff(names(dots), known)
    if( length(unknown) ){
      stop("boundary_delineation: unknown argument(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    fetch_args <- dots[names(dots) %in% fetch_names]
    inline <- dots[names(dots) %in% setdiff(param_names, fetch_names)]
    rs <- do.call(get_rs, c(list(longitude = x, params = params), fetch_args))
    params <- NULL
  }

  # Generic thresholds and sensitivities fill in only the stage overrides
  # that the caller has not supplied explicitly.
  if( !is.null(vi_sensitivity) ){
    if( is.null(inline$area_sensitivity) ){
      inline$area_sensitivity <- vi_sensitivity
    }
    if( !"interior_sensitivity" %in% names(inline) ){
      inline$interior_sensitivity <- vi_sensitivity
    }
    if( !"split_sensitivity" %in% names(inline) ){
      inline$split_sensitivity <- vi_sensitivity
    }
  }
  if( !is.null(vi_threshold) ){
    if( is.null(inline$area_threshold) ){
      inline$area_threshold <- vi_threshold
    }
    if( is.null(inline$interior_threshold) ){
      inline$interior_threshold <- vi_threshold
    }
  }
  # Pass the requested strictness to every stage.
  if( !is.null(strictness) ){
    inline$strictness <- strictness
  }

  stage <- function(fn, rs_in, params){

    return(do.call(fn, c(list(rs_in, params = params), inline)))
  }

  invocation <- .beginCalibrationInvocation(rs)
  rs <- invocation$rs
  if( invocation$owner ){
    on.exit({
      invocation$context$active <- FALSE
      invocation$context$failures <- NULL
    }, add = TRUE)
  }

  if( !is.null(rs$state$alive_mat) ){
    rs$state$decided_mat <- matrix(FALSE, rs$geom$nr, rs$geom$nc)
    attr(rs, "pending") <- list(pre_cleanup = NULL)
  }

  # Initial fill and refinement
  rs <- stage(trace_perimeter,  rs, params)
  rs <- stage(segment_area,     rs, params)
  rs <- stage(trace_perimeter,  rs, params)
  rs <- stage(split_area,       rs, params)
  rs <- stage(trace_perimeter,  rs, params)
  rs <- stage(segment_interior, rs, params)
  rs <- stage(trace_perimeter,  rs, params)

  # Additional split and cleanup pairs
  for( i in seq_len(cleanup_boundary) ){
    rs <- stage(split_area,      rs, params)
    rs <- stage(trace_perimeter, rs, params)
  }

  # Final split and cleanup pairs
  rs <- stage(split_area,      rs, params)
  closing_params <- params %||% list()
  closing_params$perimeter_margins <- c(0.3, 0.9)
  rs <- stage(trace_perimeter, rs,
              closing_params)
  rs <- stage(split_area,      rs, params)
  closing_params$perimeter_margins <- c(0.5, 0.5, 1)
  rs <- stage(trace_perimeter, rs,
              closing_params)
  return(.endCalibrationInvocation(rs, invocation))
}
