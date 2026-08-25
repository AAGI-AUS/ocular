# =========================================================================
# Segmentation: perimeter tracing
# =========================================================================

#' Cleanup rule module (internal)
#'
#' Rule-module form of boundary-cleanup logic. \code{.classifyTrail} processes
#' each complete perimeter, so the rule uses setup and finalize without a
#' per-cell callback. The driver repeats passes up to max_iter and stops at
#' convergence.
#'
#' \code{walker_setup} creates a shared area matrix and threshold predicates
#' that remain available across iterations and candidate areas.
#'
#' \code{setup} resets the previous area's cells in the shared
#' scratch matrix, paints the current area, and runs the perimeter
#' walk. \code{finalize} runs \code{.classifyTrail} and returns the
#' alive_mat / decided_mat deltas.
#'
#' @noRd
.cleanupRule <- function(){
  list(
    name       = "cleanup",
    axes       = NULL,             ## no per-axis loop
    apply_cell = NULL,              ## finalize-only

    max_iter = function(params){
      as.integer(params$trace_iter)
    },

    walker_setup = function(rs, params){
      preds <- .thresholdPredicates(params$perimeter_margins)
      list(
        area_mat   = matrix(FALSE, rs$geom$nr, rs$geom$nc),
        prev_cells = NULL,
        is_below   = preds$is_below,
        is_above   = preds$is_above
      )
    },

    build_walk = function(area, alive_mat, decided_mat, rs, params, walker_state){
      if( area$size < 3L ) return(NULL)
      perim <- .areaPerimeter(area$cells)
      if( nrow(perim) < 3L ) return(NULL)
      list(perim = perim, area_cells = area$cells)
    },

    setup = function(area, walk, alive_mat, decided_mat, rs, params, walker_state, log_msg){
      if( is.null(walk) ) return(list(skip = TRUE))
      ## Reset previous area's cells in the shared scratch, paint current.
      if( !is.null(walker_state$prev_cells) )
        walker_state$area_mat[walker_state$prev_cells] <- FALSE
      walker_state$area_mat[walk$area_cells] <- TRUE
      walker_state$prev_cells <- walk$area_cells
      wk <- .walkPerimeter(walker_state$area_mat)
      if( length(wk$edges) == 0L )
        return(list(skip = TRUE, .walker_state = walker_state))
      list(
        skip          = FALSE,
        walk          = wk,
        perim         = walk$perim,
        area_cells    = walk$area_cells,
        .walker_state = walker_state
      )
    },

    finalize = function(state, alive_mat, decided_mat, area, walk, rs, params, walker_state, log_msg){
      ## state$walk is the .walkPerimeter output for this area; walker_state
      ## carries the threshold predicates (built once in walker_setup).
      cl <- .classifyTrail(state$walk, alive_mat, decided_mat,
                           state$perim, state$area_cells,
                           walker_state$is_below, walker_state$is_above)
      changed <- FALSE
      if( nrow(cl$to_remove) > 0L ){
        alive_mat[cl$to_remove]   <- FALSE
        decided_mat[cl$to_remove] <- TRUE
        changed <- TRUE
      }
      if( nrow(cl$to_revive) > 0L ){
        alive_mat[cl$to_revive]   <- TRUE
        decided_mat[cl$to_revive] <- TRUE
        changed <- TRUE
      }
      if( !changed ) return(NULL)
      list(alive_mat = alive_mat, decided_mat = decided_mat)
    },

    walker_teardown = NULL,

    postprocess = function(rs, params, log_msg){
      ## Output rank-1 trim: cleanup iterations can erode connecting
      ## cells and split a single component into many. Honours
      ## area_separation_strict. Skipped when multiple_areas = TRUE.
      if( !isTRUE(params$multiple_areas) && any(rs$state$alive_mat) ){
        info <- .applyStrictSeparation(rs$state$alive_mat, rs$geom$centre_rc,
                                       params$area_separation_strict)
        if( info$n_areas > 1L ){
          rs$state$alive_mat <- info$label_mat == info$area_order[1L]
        }
      }
      return(rs)
    }
  )
}

#' Refine a candidate field perimeter
#'
#' Walks the perimeter of a binary candidate mask and removes or revives
#' boundary cells according to local perimeter occupancy. In an \code{ocular}
#' pipeline, a call before \code{segment_area()} stages one cleanup for the
#' subsequently filled mask; a call after segmentation refines the current mask
#' immediately. With a \code{SpatRaster}, cleanup is applied directly.
#'
#' In pipeline form, parameters stored in \code{rs$params} provide the base
#' settings. Named arguments supplied here override the corresponding setting
#' for this call only.
#'
#' @param x A \code{SpatRaster}, an \code{ocular} object, or longitude in
#'   standalone retrieval form.
#' @param perimeter_margins Lower and upper perimeter-occupancy thresholds, with
#'   an optional third value controlling which decision includes equality. In
#'   piped form, \code{NULL} uses the effective
#'   \code{params$perimeter_margins}; a \code{SpatRaster} call requires a value.
#' @param trace_iter Maximum cleanup iterations. In piped form, an omitted value
#'   uses \code{params$trace_iter}; a \code{SpatRaster} call defaults to one
#'   iteration.
#' @param strictness Calibration control for \code{ocular}-object and longitude
#'   forms. It can update the calibration record used by subsequent area growth
#'   but does not directly alter perimeter cleanup. See \code{?rs_params}.
#' @param params Optional, possibly partial \code{rs_params()} list for
#'   \code{ocular}-object or longitude forms. In piped form it overlays
#'   \code{rs$params}; in standalone retrieval form it is passed to
#'   \code{get_rs()}. It is not used for \code{SpatRaster} input.
#' @param ... Piped form accepts further named \code{rs_params()} overrides,
#'   including \code{multiple_areas}. \code{SpatRaster} form accepts only
#'   \code{multiple_areas}. Standalone retrieval form forwards remaining
#'   arguments to \code{get_rs()}.
#' @returns A cleaned \code{SpatRaster} for raster input; otherwise an
#'   \code{ocular} object with cleanup staged or applied.
#' @export
trace_perimeter <- function(x, perimeter_margins = NULL, trace_iter = 1L,
                            strictness = NULL,
                            params = NULL,
                            ...){
  ## Standalone SpatRaster form.
  ## Builds a synthetic rs + params, routes through .runStagedWalk so the
  ## SpatRaster path uses the same execution model as the piped path.
  if( inherits(x, "SpatRaster") ){
    if( is.null(perimeter_margins) )
      stop("trace_perimeter standalone form requires perimeter_margins.",
           call. = FALSE)
    .validateBoundaryCleanup(perimeter_margins)
    if( !.isWholeNumber(trace_iter) || trace_iter < 1L )
      stop("trace_iter must be a positive integer.", call. = FALSE)
    dots <- list(...)
    if( length(dots) > 0L &&
        (is.null(names(dots)) || anyNA(names(dots)) ||
         any(!nzchar(names(dots))) || anyDuplicated(names(dots))) )
      stop("trace_perimeter: all arguments in ... must be uniquely named.",
           call. = FALSE)
    unknown <- setdiff(names(dots), "multiple_areas")
    if( length(unknown) )
      stop("trace_perimeter: unused standalone argument(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    mf <- dots$multiple_areas %||% FALSE
    if( !is.logical(mf) || length(mf) != 1L || is.na(mf) )
      stop("trace_perimeter: multiple_areas must be TRUE or FALSE.",
           call. = FALSE)

    r1 <- x[[1L]]
    base_vals <- terra::values(r1, mat = FALSE)
    is_alive_value <- !is.na(base_vals) & base_vals != 0
    if( !any(is_alive_value) ) return(x)

    nr <- terra::nrow(r1); nc <- terra::ncol(r1)
    centre_rc <- c(as.integer(ceiling(nr / 2)), as.integer(ceiling(nc / 2)))
    alive_mat <- matrix(is_alive_value, nr, nc, byrow = TRUE)

    rs_syn <- .newOcular()
    rs_syn$geom$nr        <- nr
    rs_syn$geom$nc        <- nc
    rs_syn$geom$centre_rc <- centre_rc
    rs_syn$state$alive_mat   <- alive_mat
    rs_syn$state$decided_mat <- matrix(FALSE, nr, nc)
    p_syn  <- list(perimeter_margins       = perimeter_margins,
                   trace_iter            = as.integer(trace_iter),
                   multiple_areas         = mf,
                   area_separation_strict = FALSE)

    res_rs <- .runCleanupOnRs(rs_syn, p_syn,
                              log_msg = function(...) invisible())
    alive_mask <- terra::rast(x[[1L]])
    terra::values(alive_mask) <- as.integer(t(res_rs$state$alive_mat))
    alive_mask <- terra::ifel(alive_mask == 1L, 1L, NA)
    return(terra::mask(x, alive_mask))
  }

  if( is_rs(x) ){
    rs <- x
    params <- .resolveParams(rs, params,
                             c(list(strictness = strictness), list(...)))
  }else{
    rs <- get_rs(longitude = x, params = params, ...)
    params <- .resolveParams(rs, NULL,
                             list(strictness = strictness))
  }
  rs <- .applyScenePriors(rs, params)

  ## Resolve inline overrides against params defaults. Inline
  ## perimeter_margins (or perimeter_margins left NULL but
  ## params$perimeter_margins set) provides the operative values for this
  ## call only -- they don't mutate params. Same for trace_iter: an
  ## inline value wins whenever supplied (missing() test -- value equality
  ## against the default 1L cannot distinguish "user typed 1L" from
  ## "not supplied"); otherwise fall back to params$trace_iter.
  effective_thresholds <- if( !is.null(perimeter_margins) ) perimeter_margins
  else params$perimeter_margins
  if( !is.null(effective_thresholds) )
    .validateBoundaryCleanup(effective_thresholds)
  effective_iter <- if( missing(trace_iter) ) params$trace_iter else trace_iter
  if( !.isWholeNumber(effective_iter) || effective_iter < 1L )
    stop("trace_perimeter: trace_iter must be a positive integer.",
         call. = FALSE)
  effective_iter <- as.integer(effective_iter)

  ## Pre-mode: no alive_mat yet. Cache perimeter_margins/iter onto attr(rs,
  ## "pending") so segment_area runs cleanup with the user's intent at THIS
  ## call site (not whatever rs_params defaults happen to be at segment_area).
  if( is.null(rs$state$alive_mat) ){
    if( !is.null(effective_thresholds) ){
      ## Only the latest cleanup staged before segment_area can be applied.
      ## Warn when a second call replaces the first.
      if( isTRUE(attr(rs, "pending")$pre_cleanup$active) )
        warning("trace_perimeter(): staged cleanup settings are being replaced. ",
                "Only the most recent trace_perimeter() call before ",
                "segment_area() is applied. Move one call after segment_area() ",
                "if two cleanup passes are intended.",
                call. = FALSE)
      attr(rs, "pending")$pre_cleanup$active     <- TRUE
      attr(rs, "pending")$pre_cleanup$thresholds <- effective_thresholds
      attr(rs, "pending")$pre_cleanup$iter       <- effective_iter
    }
    return(rs)
  }
  ## Post-mode
  if( is.null(effective_thresholds) ) return(rs)
  if( !any(rs$state$alive_mat) ) return(rs)

  log_msg <- function(...) invisible()

  ## Synthetic params for the rule, with the call-site-resolved values.
  p_call <- params
  p_call$perimeter_margins <- effective_thresholds
  p_call$trace_iter      <- effective_iter

  rs <- .runCleanupOnRs(rs, p_call, log_msg)
  log_msg(sprintf("After post-cleanup: %d alive", sum(rs$state$alive_mat)))
  return(rs)
}

#' Run cleanup rule on an rs (internal)
#'
#' Wraps the input single-area rank-1 trim and the .runStagedWalk
#' invocation. Used by both the piped and SpatRaster trace_perimeter
#' forms (the SpatRaster form passes a synthetic rs).
#'
#' @noRd
.runCleanupOnRs <- function(rs, params, log_msg = function(...) invisible()){
  if( !any(rs$state$alive_mat) ) return(rs)

  ## Single-area input filter (honours area_separation_strict). This
  ## sits outside the rule's postprocess because it operates on the
  ## input alive set, before any cleanup walks.
  if( !isTRUE(params$multiple_areas) ){
    info <- .applyStrictSeparation(rs$state$alive_mat, rs$geom$centre_rc,
                                   params$area_separation_strict)
    if( info$n_areas > 1L ){
      rs$state$alive_mat <- info$label_mat == info$area_order[1L]
    }
  }

  rule <- .cleanupRule()
  rs   <- .runStagedWalk(rs, rule, params, log_msg)
  rs   <- rule$postprocess(rs, params, log_msg)
  return(rs)
}
