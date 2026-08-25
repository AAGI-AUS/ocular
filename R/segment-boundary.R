# =========================================================================
# Segmentation: boundary delineation
# =========================================================================

#' Delineate a field boundary
#'
#' Runs ocular's default field boundary delineation workflow. It grows a
#' candidate field around the query point, refines its perimeter and interior,
#' and repeatedly attempts to separate adjoining fields. The stage order is
#' fixed; \code{cleanup_boundary} controls only the number of additional
#' split-and-cleanup pairs in the middle of the workflow.
#'
#' Parameters stored in \code{rs$params} provide the base settings for every
#' stage. Supply
#' \code{params} to overlay settings for the whole workflow, or named
#' arguments in \code{...} to override individual settings for this call only.
#'
#' @param x An \code{ocular} object in piped form, or longitude in standalone
#'   form.
#' @param vi_sensitivity Generic sensitivity override. When supplied, fans
#'   out to \code{area_sensitivity}, \code{interior_sensitivity} and
#'   \code{split_sensitivity} for any not set explicitly via \code{...}.
#'   Not an \code{rs_params()} parameter.
#' @param vi_threshold Generic threshold override. When supplied, fans out
#'   to \code{area_threshold} and \code{interior_threshold} for any not set
#'   explicitly via \code{...}. Not an \code{rs_params()} parameter.
#' @param strictness Calibration control forwarded to every stage. It selects
#'   the robustness multiple applied to measured per-window dispersion. See
#'   \code{?rs_params}.
#' @param cleanup_boundary Integer in \code{[0, 5]} giving the number of
#'   additional split-and-cleanup pairs in the middle phase. Default \code{2L}.
#' @param params Optional, possibly partial \code{rs_params()} list. Overlays
#'   \code{rs$params} for the whole call. \code{NULL} keeps the stored settings.
#' @param ... Named \code{rs_params()} overrides forwarded to every stage.
#'   A stage-specific value here takes precedence over the generic
#'   \code{vi_sensitivity}/\code{vi_threshold}. Standalone mode also forwards
#'   spatial and temporal retrieval arguments to \code{get_rs()}.
#' @returns An \code{ocular} object with its delineated field mask stored in
#'   \code{state$alive_mat}.
#' @export
boundary_delineation <- function(x,
                                 vi_sensitivity   = NULL,
                                 vi_threshold     = NULL,
                                 strictness       = NULL,
                                 cleanup_boundary = 2L,
                                 params           = NULL,
                                 ...){

  if( !.isWholeNumber(cleanup_boundary) ||
      cleanup_boundary < 0L || cleanup_boundary > 5L )
    stop("boundary_delineation: cleanup_boundary must be a single integer in [0, 5].",
         call. = FALSE)
  cleanup_boundary <- as.integer(cleanup_boundary)

  dots <- list(...)
  if( length(dots) > 0L &&
      (is.null(names(dots)) || anyNA(names(dots)) ||
       any(!nzchar(names(dots))) || anyDuplicated(names(dots))) )
    stop("boundary_delineation: all arguments in ... must be uniquely named.",
         call. = FALSE)
  param_names <- names(rs_params())

  if( is_rs(x) ){
    rs <- x
    unknown <- setdiff(names(dots), param_names)
    if( length(unknown) )
      stop("boundary_delineation: unknown parameter(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    inline <- dots
  }else{
    ## Standalone calls separate retrieval arguments from stage parameter
    ## overrides so only retrieval arguments are passed to get_rs().
    fetch_names <- setdiff(names(formals(get_rs)), c("longitude", "params"))
    known <- union(fetch_names, param_names)
    unknown <- setdiff(names(dots), known)
    if( length(unknown) )
      stop("boundary_delineation: unknown argument(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    fetch_args <- dots[names(dots) %in% fetch_names]
    inline <- dots[names(dots) %in% setdiff(param_names, fetch_names)]
    rs <- do.call(get_rs, c(list(longitude = x, params = params), fetch_args))
    params <- NULL
  }

  ## Generic convenience args fan out to the per-stage areas. An explicit
  ## per-stage value passed via ... wins; the generic fills the rest; both
  ## override rs_params()/calibration. vi_sensitivity and vi_threshold are not
  ## rs_params() areas -- they exist only as these boundary_delineation args.
  if( !is.null(vi_sensitivity) ){
    if( is.null(inline$area_sensitivity) )
      inline$area_sensitivity <- vi_sensitivity
    if( !"interior_sensitivity" %in% names(inline) )
      inline$interior_sensitivity <- vi_sensitivity
    if( !"split_sensitivity" %in% names(inline) )
      inline$split_sensitivity <- vi_sensitivity
  }
  if( !is.null(vi_threshold) ){
    if( is.null(inline$area_threshold) )
      inline$area_threshold <- vi_threshold
    if( is.null(inline$interior_threshold) )
      inline$interior_threshold <- vi_threshold
  }
  ## Forward the call-level calibration setting to every stage; each stage
  ## threads it through .applyScenePriors before resolving params.
  if( !is.null(strictness) ) inline$strictness <- strictness

  stage <- function(fn, rs_in, params){
    do.call(fn, c(list(rs_in, params = params), inline))
  }

  ## Initial phase
  rs <- stage(trace_perimeter,  rs, params) # pre-cleanup
  rs <- stage(segment_area,     rs, params) # detect area
  rs <- stage(trace_perimeter,  rs, params) # post-cleanup
  rs <- stage(split_area,       rs, params) # default split
  rs <- stage(trace_perimeter,  rs, params) # default cleanup
  rs <- stage(segment_interior, rs, params) # in-area object detection
  rs <- stage(trace_perimeter,  rs, params) # default cleanup

  ## Middle phase: cleanup_boundary split + cleanup pairs (default 2)
  for( i in seq_len(cleanup_boundary) ){
    rs <- stage(split_area,      rs, params)
    rs <- stage(trace_perimeter, rs, params)
  }

  ## Closing phase (split x2, cleanup x2)
  rs <- stage(split_area,      rs, params)
  closing_params <- params %||% list()
  closing_params$perimeter_margins <- c(0.3, 0.9)
  rs <- stage(trace_perimeter, rs,
              closing_params)
  rs <- stage(split_area,      rs, params)
  closing_params$perimeter_margins <- c(0.5, 0.5, 1)
  rs <- stage(trace_perimeter, rs,
              closing_params)
  return(rs)
}
