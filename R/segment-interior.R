# =========================================================================
# Segmentation: interior
# =========================================================================

#' Measure a row or column across an area
#'
#' Return the line of cells, their alive median and whether measurement failed.
#' @noRd
.scanLOS <- function(envelope, pos_now, dir_now, alive_mat, feat_array, nr, nc,
                     cell_medians = NULL){

  skip_result <- list(median_alive = NA_real_, los_rc = NULL, skip = TRUE)
  axis <- envelope$axis
  valid_dirs <- if( axis == "y" ) { c("N", "S") } else { c("E", "W") }
  if( !(dir_now %in% valid_dirs) ){
    return(skip_result)
  }

  near <- envelope$near_rc
  far <- envelope$far_rc
  if( axis == "y" ){
    idx <- which(near[, "row"] == pos_now[1L])
    if( length(idx) == 0L ){
      return(skip_result)
    }
    idx <- idx[1L]
    lo <- min(near[idx, "col"], far[idx, "col"])
    hi <- max(near[idx, "col"], far[idx, "col"])
    los_rc <- cbind(row = rep(pos_now[1L], hi - lo + 1L), col = lo:hi)
  }else{
    idx <- which(near[, "col"] == pos_now[2L])
    if( length(idx) == 0L ){
      return(skip_result)
    }
    idx <- idx[1L]
    lo <- min(near[idx, "row"], far[idx, "row"])
    hi <- max(near[idx, "row"], far[idx, "row"])
    los_rc <- cbind(row = lo:hi, col = rep(pos_now[2L], hi - lo + 1L))
  }
  alive_along <- alive_mat[los_rc]
  if( !any(alive_along) ){
    return(skip_result)
  }
  alive_idx <- los_rc[alive_along, , drop = FALSE]
  nw <- dim(feat_array)[3L]
  if( is.null(nw) || is.na(nw) ){
    vals <- feat_array[alive_idx]
  }else if( nw == 1L ){
    vals <- feat_array[cbind(alive_idx[, "row"], alive_idx[, "col"], 1L)]
  }else if( !is.null(cell_medians) ){
    vals <- cell_medians[alive_idx]
    missing_idx <- which(is.na(vals))
    vals[missing_idx] <- vapply(missing_idx, function(i){
      strip <- feat_array[alive_idx[i, "row"], alive_idx[i, "col"], ]
      return(median(strip, na.rm = TRUE))
    }, numeric(1))
  }else{
    per_cell <- vapply(seq_len(nrow(alive_idx)), function(i){
      strip <- feat_array[alive_idx[i, "row"], alive_idx[i, "col"], ]
      return(median(strip, na.rm = TRUE))
    }, numeric(1))
    vals <- per_cell
  }
  vals <- vals[is.finite(vals)]
  if( length(vals) == 0L ){
    return(skip_result)
  }
  return(list(median_alive = as.numeric(median(vals)),
              los_rc = los_rc,
              skip = FALSE))
}

#' Refine an area's interior by scanning rows and columns
#'
#' Evaluate both axes from the same mask, then combine their changes. Within
#' each span, remove alive cells outside \code{interior_threshold}. If its
#' median differs sufficiently from the seed, reclassify interior cells against
#' the seed signature. Each scan leaves its two envelope endpoints unchanged.
#'
#' When only one area is requested, retain the first ranked component using
#' the requested connectivity.
#' @noRd
.runDetectObjects <- function(rs, params, log_msg = function(...) invisible()){

  if( !any(rs$state$alive_mat) ){
    return(rs)
  }

  alive_mat <- rs$state$alive_mat
  decided_mat <- if( is.null(rs$state$decided_mat) ){
    matrix(FALSE, rs$geom$nr, rs$geom$nc)
  }else{
    rs$state$decided_mat
  }
  nr <- rs$geom$nr
  nc <- rs$geom$nc

  tol <- params$interior_sensitivity
  thr_lo <- params$interior_threshold[1L]
  thr_hi <- params$interior_threshold[2L]
  feat_arr <- rs$internals$feat_array
  mu_seed <- if( length(rs$internals$mu) > 1L ){
    median(rs$internals$mu, na.rm = TRUE)
  }else{
    as.numeric(rs$internals$mu)
  }

  info <- .enumerateAreas(alive_mat, rs$geom$centre_rc, strict = FALSE)
  if( info$n_areas == 0L ){
    rs$state$alive_mat <- alive_mat
    rs$state$decided_mat <- decided_mat
    return(rs)
  }
  cell_medians <- .windowMedians(feat_arr, alive_mat)

  scan_and_reclassify <- function(env, alive_mat, decided_mat){

    if( is.null(env) ){
      return(list(alive_mat = alive_mat,
                  decided_mat = decided_mat))
    }
    near <- env$near_rc
    far <- env$far_rc
    n_anchors <- nrow(near)
    for( i in seq_len(n_anchors) ){
      if( env$axis == "y" ){
        rv <- near[i, "row"]
        lo <- min(near[i, "col"], far[i, "col"])
        hi <- max(near[i, "col"], far[i, "col"])
        if( hi - lo < 2L ){
          next
        }
        span_rc <- cbind(row = rep(rv, hi - lo + 1L), col = lo:hi)
      }else{
        cv <- near[i, "col"]
        lo <- min(near[i, "row"], far[i, "row"])
        hi <- max(near[i, "row"], far[i, "row"])
        if( hi - lo < 2L ){
          next
        }
        span_rc <- cbind(row = lo:hi, col = rep(cv, hi - lo + 1L))
      }
      alive_along <- alive_mat[span_rc]
      if( !any(alive_along) ){
        next
      }

      # Remove threshold failures before measuring divergence.
      # Leave this span's endpoints unchanged.
      interior_rc <- span_rc[-c(1L, nrow(span_rc)), , drop = FALSE]
      if( nrow(interior_rc) > 0L ){
        for( j in seq_len(nrow(interior_rc)) ){
          r <- interior_rc[j, "row"]
          c <- interior_rc[j, "col"]
          if( !alive_mat[r, c] ){
            next
          }
          if( decided_mat[r, c] ){
            next
          }
          v <- .viAt(feat_arr, r, c, cell_medians)
          if( !is.finite(v) ){
            next
          }
          if( v < thr_lo || v > thr_hi ){
            alive_mat[r, c] <- FALSE
          }
        }
      }

      alive_along <- alive_mat[span_rc]
      if( !any(alive_along) ){
        next
      }

      alive_idx <- span_rc[alive_along, , drop = FALSE]
      nw <- dim(feat_arr)[3L]
      if( is.null(nw) || is.na(nw) ){
        vals <- feat_arr[alive_idx]
      }else if( nw == 1L ){
        vals <- feat_arr[cbind(alive_idx[, "row"], alive_idx[, "col"], 1L)]
      }else if( !is.null(cell_medians) ){
        vals <- cell_medians[alive_idx]
        missing_idx <- which(is.na(vals))
        vals[missing_idx] <- vapply(missing_idx, function(k){
          strip <- feat_arr[alive_idx[k, "row"], alive_idx[k, "col"], ]
          return(median(strip, na.rm = TRUE))
        }, numeric(1))
      }else{
        vals <- vapply(seq_len(nrow(alive_idx)), function(k){
          strip <- feat_arr[alive_idx[k, "row"], alive_idx[k, "col"], ]
          return(median(strip, na.rm = TRUE))
        }, numeric(1))
      }
      vals <- vals[is.finite(vals)]
      if( length(vals) == 0L ){
        next
      }
      median_alive <- median(vals)
      if( abs(median_alive - mu_seed) <= tol ){
        next
      }

      # Check the range again so reclassification cannot revive an excluded value.
      if( nrow(interior_rc) == 0L ){
        next
      }
      for( j in seq_len(nrow(interior_rc)) ){
        r <- interior_rc[j, "row"]
        c <- interior_rc[j, "col"]
        if( decided_mat[r, c] ){
          next
        }
        v <- .viAt(feat_arr, r, c, cell_medians)
        if( !is.finite(v) ){
          next
        }
        in_tol <- abs(v - mu_seed) <= tol &&
          v >= thr_lo && v <= thr_hi
        cell_alive <- alive_mat[r, c]
        if( in_tol && !cell_alive ){
          alive_mat[r, c] <- TRUE
        }else if( !in_tol && cell_alive ){
          alive_mat[r, c] <- FALSE
        }
      }
    }
    return(list(alive_mat = alive_mat, decided_mat = decided_mat))
  }

  for( ai in seq_len(info$n_areas) ){
    area <- info$areas[[ai]]
    if( area$size < 3L ){
      next
    }
    envelopes <- .firstShotEnvelope(area$cells, nr, nc)

    # Scan both axes from the same mask, then combine their changes.
    alive_snap <- alive_mat
    out_y <- scan_and_reclassify(envelopes$y, alive_snap, decided_mat)
    out_x <- scan_and_reclassify(envelopes$x, alive_snap, decided_mat)
    y_chg <- out_y$alive_mat != alive_snap
    x_chg <- out_x$alive_mat != alive_snap
    alive_mat[y_chg] <- out_y$alive_mat[y_chg]
    alive_mat[x_chg] <- out_x$alive_mat[x_chg]
  }

  if( isFALSE(params$multiple_areas) && any(alive_mat) ){
    sep <- .applyStrictSeparation(alive_mat, rs$geom$centre_rc,
                                  params$area_separation_strict)
    if( sep$n_areas > 1L ){
      keep_label <- sep$area_order[1L]
      alive_mat <- sep$label_mat == keep_label
      log_msg(sprintf("Single-area: keeping rank-1 area (%d px)",
                      sum(alive_mat)))
    }
  }
  log_msg(sprintf("segment_interior: %d alive", sum(alive_mat)))

  rs$state$alive_mat <- alive_mat
  rs$state$decided_mat <- decided_mat
  return(rs)
}

#' Convert a row and column step to a compass direction
#' @noRd
.deltaToDir <- function(dr, dc){

  if( dr == -1L && dc == 0L ){
    return("N")
  }
  if( dr == 1L && dc == 0L ){
    return("S")
  }
  if( dr == 0L && dc == 1L ){
    return("E")
  }
  if( dr == 0L && dc == -1L ){
    return("W")
  }
  return(NA_character_)
}

#' Read a cell's index value, taking the median across windows
#' @noRd
.viAt <- function(feat_array, r, c, cell_medians = NULL){

  if( length(dim(feat_array)) == 3L ){
    if( !is.null(cell_medians) ){
      v <- cell_medians[r, c]
      if( !is.na(v) ){
        return(v)
      }
    }
    strip <- feat_array[r, c, ]
    if( length(strip) == 1L ){
      return(strip)
    }
    return(median(strip, na.rm = TRUE))
  }
  return(feat_array[r, c])
}

#' Cache medians across windows
#'
#' Fill selected cells on a matching grid with at least two double-valued
#' windows. For three windows, the middle of three finite values also fills
#' unselected cells unless it is zero. Other cells stay \code{NA}, so callers
#' calculate their median directly. Return \code{NULL} for unsuitable inputs.
#' @noRd
.windowMedians <- function(feat_array, cells_mat){

  feat_dim <- dim(feat_array)
  if( length(feat_dim) != 3L || feat_dim[3L] < 2L || !is.double(feat_array) ||
      !identical(as.integer(feat_dim[1:2]), as.integer(dim(cells_mat))) ){
    return(NULL)
  }
  nr <- feat_dim[1L]
  out <- matrix(NA_real_, nr, feat_dim[2L])
  todo <- cells_mat
  if( feat_dim[3L] == 3L ){
    v1 <- feat_array[, , 1L]
    v2 <- feat_array[, , 2L]
    v3 <- feat_array[, , 3L]
    middle <- pmax(pmin(v1, v2), pmin(pmax(v1, v2), v3))
    # Use median() for zero to preserve its sign.
    selected <- is.finite(v1) & is.finite(v2) & is.finite(v3) & middle != 0
    out[selected] <- middle[selected]
    todo <- todo & !selected
  }
  idx <- which(todo)
  rows <- (idx - 1L) %% nr + 1L
  cols <- (idx - 1L) %/% nr + 1L
  out[idx] <- vapply(seq_along(idx), function(i){
    return(median(feat_array[rows[i], cols[i], ], na.rm = TRUE))
  }, numeric(1))
  return(out)
}

#' Refine a candidate field interior
#'
#' Examines line-of-sight profiles within each candidate field and revises
#' interior pixels whose spectral behaviour is inconsistent with the seed
#' signature. It can remove pixels on the perimeter, and revives pixels only
#' between the outermost alive pixels of each row and column.
#' \code{trace_perimeter()} refines the boundary itself. If no candidate mask
#' exists, \code{segment_area()} runs first.
#'
#' @param x An \code{ocular} object in piped form, or longitude in standalone
#'   form.
#' @param params Optional list of settings for this call. See
#'   \code{\link{rs_params}} for the settings and override rules.
#' @param interior_sensitivity Maximum permitted difference between a
#'   line-of-sight alive-pixel median and the seed signature before interior
#'   cells are reclassified. An explicit \code{NULL} disables interior
#'   refinement, and an omitted argument uses the stored setting.
#' @param interior_threshold Inclusive lower and upper spectral index bounds.
#'   Alive interior cells outside the range are removed. During a divergent
#'   line-of-sight pass, in-range cells are retained or revived only when they
#'   also match the seed signature within \code{interior_sensitivity}.
#' @param strictness Calibration control. It affects initial area growth when
#'   \code{segment_area()} must run and updates the calibration record for
#'   subsequent stages. It does not directly change the interior predicate. See
#'   \code{?rs_params}.
#' @param ... Standalone mode forwards remaining arguments to \code{get_rs()}.
#'   Piped mode accepts further named \code{rs_params()} overrides for this
#'   call. Unknown names are rejected.
#' @returns An \code{ocular} object with its revised candidate field mask stored
#'   in \code{state$alive_mat}.
#' @export
segment_interior <- function(x,
                             interior_sensitivity = NULL,
                             interior_threshold = NULL,
                             strictness = NULL,
                             params = NULL,
                             ...){

  disable_interior <- !missing(interior_sensitivity) &&
    is.null(interior_sensitivity)
  formal_overrides <- list(
    interior_sensitivity = interior_sensitivity,
    interior_threshold = interior_threshold,
    strictness = strictness
  )
  input_params <- params
  if( is_rs(x) ){
    rs <- x
    params <- .resolveParams(rs, params, c(formal_overrides, list(...)),
                             caller = "segment_interior()")
  }else{
    rs <- get_rs(longitude = x, params = params, ...)
    params <- .resolveParams(rs, NULL, formal_overrides, caller = "segment_interior()")
  }
  if( disable_interior ){
    params["interior_sensitivity"] <- list(NULL)
    attr(params, "user_set") <- union(
      attr(params, "user_set", exact = TRUE), "interior_sensitivity")
  }
  invocation <- .beginCalibrationInvocation(rs)
  rs <- invocation$rs
  if( invocation$owner ){
    on.exit({
      invocation$context$active <- FALSE
      invocation$context$failures <- NULL
    }, add = TRUE)
  }

  rs <- .applyScenePriors(rs, params)
  if( is.null(rs$state$alive_mat) ){
    # Pass the original settings to filling so refreshed calibration can apply.
    # Standalone retrieval has already stored these settings.
    if( is_rs(x) ){
      rs <- segment_area(rs,
                         strictness = strictness,
                         interior_sensitivity = interior_sensitivity,
                         interior_threshold = interior_threshold,
                         params = input_params,
                         ...)
    }else{
      rs <- segment_area(rs,
                         strictness = strictness,
                         interior_sensitivity = interior_sensitivity,
                         interior_threshold = interior_threshold)
    }
  }
  if( is.null(params$interior_sensitivity) ){
    return(.endCalibrationInvocation(rs, invocation))
  }
  if( !any(rs$state$alive_mat) ){
    return(.endCalibrationInvocation(rs, invocation))
  }
  if( is.null(rs$internals$feat_array) || is.null(rs$internals$mu) ){
    stop("segment_interior(): detection state is unavailable. Run ",
         "segment_area() first.", call. = FALSE)
  }

  log_msg <- function(...){

    return(invisible())
  }
  mu_seed <- if( length(rs$internals$mu) > 1L ){
    median(rs$internals$mu, na.rm = TRUE)
  }else{
    as.numeric(rs$internals$mu)
  }
  log_msg(sprintf("segment_interior: tol=%.3f mu_seed=%.3f threshold=[%.3f, %.3f]",
                  params$interior_sensitivity, mu_seed,
                  params$interior_threshold[1L],
                  params$interior_threshold[2L]))

  result <- .runDetectObjects(rs, params, log_msg)
  return(.endCalibrationInvocation(result, invocation))
}
