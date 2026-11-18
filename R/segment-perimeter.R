# =========================================================================
# Segmentation: perimeter tracing
# =========================================================================

#' Provide the cleanup callbacks for the perimeter walker
#'
#' Classify a complete trail in \code{finalize()}, without a callback for each
#' cell. Reuse an area mask and threshold tests across areas and iterations.
#' Clear the previous area's cells before drawing the next. The driver repeats
#' cleanup until no changes are returned or the iteration limit is reached.
#' @noRd
.cleanupRule <- function(){

  return(list(
    name = "cleanup",
    axes = NULL,
    apply_cell = NULL,

    max_iter = function(params){
      return(as.integer(params$trace_iter))
    },

    walker_setup = function(rs, params){
      preds <- .thresholdPredicates(params$perimeter_margins)
      return(list(
        area_mat = matrix(FALSE, rs$geom$nr, rs$geom$nc),
        prev_cells = NULL,
        is_below = preds$is_below,
        is_above = preds$is_above
      ))
    },

    build_walk = function(area, alive_mat, decided_mat, rs, params, walker_state){
      if( area$size < 3L ){
        return(NULL)
      }
      perim <- .areaPerimeter(area$cells)
      if( nrow(perim) < 3L ){
        return(NULL)
      }
      return(list(perim = perim, area_cells = area$cells))
    },

    setup = function(area, walk, alive_mat, decided_mat, rs, params, walker_state, log_msg){
      if( is.null(walk) ){
        return(list(skip = TRUE))
      }
      if( !is.null(walker_state$prev_cells) ){
        walker_state$area_mat[walker_state$prev_cells] <- FALSE
      }
      walker_state$area_mat[walk$area_cells] <- TRUE
      walker_state$prev_cells <- walk$area_cells
      wk <- .walkPerimeter(walker_state$area_mat)
      if( length(wk$edges) == 0L ){
        return(list(skip = TRUE, .walker_state = walker_state))
      }
      return(list(
        skip = FALSE,
        walk = wk,
        perim = walk$perim,
        area_cells = walk$area_cells,
        .walker_state = walker_state
      ))
    },

    finalize = function(state, alive_mat, decided_mat, area, walk, rs, params, walker_state,
                        log_msg){
      cl <- .classifyTrail(state$walk, alive_mat, decided_mat,
                           state$perim, state$area_cells,
                           walker_state$is_below, walker_state$is_above)
      changed <- FALSE
      if( nrow(cl$to_remove) > 0L ){
        alive_mat[cl$to_remove] <- FALSE
        decided_mat[cl$to_remove] <- TRUE
        changed <- TRUE
      }
      if( nrow(cl$to_revive) > 0L ){
        alive_mat[cl$to_revive] <- TRUE
        decided_mat[cl$to_revive] <- TRUE
        changed <- TRUE
      }
      if( !changed ){
        return(NULL)
      }
      return(list(alive_mat = alive_mat, decided_mat = decided_mat))
    },

    walker_teardown = NULL,

    postprocess = function(rs, params, log_msg){
      # Cleanup can sever connections. When one area is requested, retain the
      # first ranked component under the requested connectivity.
      if( !isTRUE(params$multiple_areas) && any(rs$state$alive_mat) ){
        info <- .applyStrictSeparation(rs$state$alive_mat, rs$geom$centre_rc,
                                       params$area_separation_strict)
        if( info$n_areas > 1L ){
          rs$state$alive_mat <- info$label_mat == info$area_order[1L]
        }
      }
      return(rs)
    }
  ))
}

#' Refine a candidate field perimeter
#'
#' Walks the perimeter of a binary candidate mask and removes or revives
#' boundary cells according to local perimeter occupancy. In an \code{ocular}
#' pipeline, a call before \code{segment_area()} stages one cleanup for the
#' subsequently filled mask. A call after segmentation refines the current mask
#' immediately. With a \code{SpatRaster}, only the removals of the cleanup are
#' applied. An area stack with layers \code{area_1}, \code{area_2}, and so on,
#' plus \code{union}, is cleaned area by area. The union of the cleaned
#' area masks is applied to the original layers, preserving their values. Other
#' raster inputs use the first layer's mask.
#'
#' @param x A \code{SpatRaster}, an \code{ocular} object, or longitude in
#'   standalone retrieval form.
#' @param perimeter_margins Lower and upper perimeter-occupancy thresholds, with
#'   an optional third value controlling which decision includes equality. In
#'   piped form, \code{NULL} uses the effective
#'   \code{params$perimeter_margins}. A \code{SpatRaster} call requires a value.
#' @param trace_iter Maximum cleanup iterations. In piped form, an omitted value
#'   uses \code{params$trace_iter}. A \code{SpatRaster} call defaults to one
#'   iteration.
#' @param strictness Calibration control for \code{ocular}-object and longitude
#'   forms. It can update the calibration record used by subsequent area growth
#'   but does not directly alter perimeter cleanup. See \code{?rs_params}.
#' @param params Optional list of settings for \code{ocular} input or standalone
#'   retrieval. See \code{\link{rs_params}} for the settings and override rules.
#'   Not used for raster input.
#' @param ... Piped form accepts further named \code{rs_params()} overrides,
#'   including \code{multiple_areas}. \code{SpatRaster} form accepts only
#'   \code{multiple_areas}. Standalone retrieval form forwards remaining
#'   arguments to \code{get_rs()}.
#' @returns For raster input, a \code{SpatRaster} with the removals applied.
#'   Otherwise, an \code{ocular} object with cleanup staged or applied.
#' @export
trace_perimeter <- function(x, perimeter_margins = NULL, trace_iter = 1L,
                            strictness = NULL,
                            params = NULL,
                            ...){

  # Use a temporary ocular object so raster masks receive the same cleanup.
  if( inherits(x, "SpatRaster") ){
    if( is.null(perimeter_margins) ){
      stop("trace_perimeter standalone form requires perimeter_margins.",
           call. = FALSE)
    }
    .validateBoundaryCleanup(perimeter_margins)
    if( !.isWholeNumber(trace_iter) || trace_iter < 1L ){
      stop("trace_iter must be a positive integer.", call. = FALSE)
    }
    dots <- list(...)
    if( length(dots) > 0L &&
        (is.null(names(dots)) || anyNA(names(dots)) ||
         any(!nzchar(names(dots))) || anyDuplicated(names(dots))) ){
      stop("trace_perimeter: all arguments in ... must be uniquely named.",
           call. = FALSE)
    }
    unknown <- setdiff(names(dots), "multiple_areas")
    if( length(unknown) ){
      stop("trace_perimeter: unused standalone argument(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    mf <- dots$multiple_areas %||% FALSE
    if( !is.logical(mf) || length(mf) != 1L || is.na(mf) ){
      stop("trace_perimeter: multiple_areas must be TRUE or FALSE.",
           call. = FALSE)
    }

    n_areas <- terra::nlyr(x) - 1L
    area_names <- paste0("area_", seq_len(n_areas))
    area_stack <- n_areas >= 2L &&
      identical(sort(names(x)), sort(c(area_names, "union")))
    mask_layers <- if( area_stack ) { match(area_names, names(x)) } else { 1L }

    nr <- terra::nrow(x)
    nc <- terra::ncol(x)
    centre_rc <- c(as.integer(ceiling(nr / 2)), as.integer(ceiling(nc / 2)))
    cleaned_union <- matrix(FALSE, nr, nc)
    p_syn <- list(perimeter_margins = perimeter_margins,
                  trace_iter = as.integer(trace_iter),
                  multiple_areas = mf,
                  area_separation_strict = FALSE)

    for( layer in mask_layers ){
      base_vals <- terra::values(x[[layer]], mat = FALSE)
      is_alive_value <- !is.na(base_vals) & base_vals != 0
      if( !any(is_alive_value) ){
        if( !area_stack ){
          return(x)
        }
        next
      }

      rs_syn <- .newOcular()
      rs_syn$geom$nr <- nr
      rs_syn$geom$nc <- nc
      rs_syn$geom$centre_rc <- centre_rc
      rs_syn$state$alive_mat <- matrix(is_alive_value, nr, nc, byrow = TRUE)
      rs_syn$state$decided_mat <- matrix(FALSE, nr, nc)
      res_rs <- .runCleanupOnRs(rs_syn, p_syn,
                                log_msg = function(...){
                                  return(invisible())
                                })
      cleaned_union <- cleaned_union | res_rs$state$alive_mat
    }

    alive_mask <- terra::rast(x[[1L]])
    terra::values(alive_mask) <- as.integer(t(cleaned_union))
    alive_mask <- terra::ifel(alive_mask == 1L, 1L, NA)
    return(terra::mask(x, alive_mask))
  }

  if( is_rs(x) ){
    rs <- x
    params <- .resolveParams(rs, params,
                             c(list(strictness = strictness), list(...)),
                             caller = "trace_perimeter()")
  }else{
    rs <- get_rs(longitude = x, params = params, ...)
    params <- .resolveParams(rs, NULL,
                             list(strictness = strictness), caller = "trace_perimeter()")
  }
  rs <- .applyScenePriors(rs, params)

  # Explicit arguments override stored settings for this call. Use missing()
  # for trace_iter so an explicitly supplied 1L still takes precedence.
  effective_thresholds <- if( !is.null(perimeter_margins) ){
    perimeter_margins
  }else{
    params$perimeter_margins
  }
  if( !is.null(effective_thresholds) ){
    .validateBoundaryCleanup(effective_thresholds)
  }
  effective_iter <- if( missing(trace_iter) ) { params$trace_iter } else { trace_iter }
  if( !.isWholeNumber(effective_iter) || effective_iter < 1L ){
    stop("trace_perimeter: trace_iter must be a positive integer.",
         call. = FALSE)
  }
  effective_iter <- as.integer(effective_iter)

  # Without a mask, save this cleanup request for segment_area(), including
  # the settings supplied here.
  if( is.null(rs$state$alive_mat) ){
    if( !is.null(effective_thresholds) ){
      # Only the latest pending cleanup can run, so warn before replacing it.
      if( isTRUE(attr(rs, "pending")$pre_cleanup$active) ){
        warning("trace_perimeter(): staged cleanup settings are being replaced. ",
                "Only the most recent trace_perimeter() call before ",
                "segment_area() is applied. Move one call after segment_area() ",
                "if two cleanup passes are intended.",
                call. = FALSE)
      }
      attr(rs, "pending")$pre_cleanup$active <- TRUE
      attr(rs, "pending")$pre_cleanup$thresholds <- effective_thresholds
      attr(rs, "pending")$pre_cleanup$iter <- effective_iter
    }
    return(rs)
  }
  if( is.null(effective_thresholds) ){
    return(rs)
  }
  if( !any(rs$state$alive_mat) ){
    return(rs)
  }

  log_msg <- function(...){

    return(invisible())
  }

  p_call <- params
  p_call$perimeter_margins <- effective_thresholds
  p_call$trace_iter <- effective_iter

  rs <- .runCleanupOnRs(rs, p_call, log_msg)
  log_msg(sprintf("After post-cleanup: %d alive", sum(rs$state$alive_mat)))
  return(rs)
}

#' Clean an ocular mask
#'
#' Select the first ranked input area when required, then run the cleanup rule.
#' Raster input reaches this helper through a temporary ocular object.
#' @noRd
.runCleanupOnRs <- function(rs, params, log_msg = function(...) invisible()){

  if( !any(rs$state$alive_mat) ){
    return(rs)
  }

  # Select the input area before walking. The later selection handles any
  # new components created by cleanup.
  if( !isTRUE(params$multiple_areas) ){
    info <- .applyStrictSeparation(rs$state$alive_mat, rs$geom$centre_rc,
                                   params$area_separation_strict)
    if( info$n_areas > 1L ){
      rs$state$alive_mat <- info$label_mat == info$area_order[1L]
    }
  }

  rule <- .cleanupRule()
  rs <- .runStagedWalk(rs, rule, params, log_msg)
  rs <- rule$postprocess(rs, params, log_msg)
  return(rs)
}
