# =========================================================================
# Segmentation: area growth
# =========================================================================

#' Grow an initial candidate field area
#'
#' Grows a candidate field from the supplied point using spectral similarity
#' across a multi-temporal feature stack. By default, it fills each selected
#' composite window separately and combines the masks. A joint fill across
#' windows is also available.
#'
#' If a mask is already present, the search is confined to that mask. Each
#' connected area starts at the eligible cell nearest the query point.
#' Calibration can supply a signature and matching tolerance for each window.
#' Without suitable calibration, the configured \code{area_sensitivity} is used.
#' When filling uses calibration, an eligible calibrated signature replaces a
#' missing learned signature for merging and later stages. Finite learned
#' signatures are retained.
#'
#' @param x An \code{ocular} object in piped form, or longitude in standalone
#'   form.
#' @param params Optional list of settings for this call. See
#'   \code{\link{rs_params}} for the settings and override rules.
#' @param area_sensitivity Call-level override for the scalar spectral matching
#'   tolerance used when calibrated per-window tolerances are not applied.
#'   \code{NULL} uses the effective stored or calibrated setting.
#' @param area_threshold Inclusive lower and upper spectral index bounds for
#'   this call.
#' @param local_alive_density Minimum proportion of a candidate cell's visited
#'   8-neighbours that must pass the spectral tests before the cell is added.
#' @param strictness Controls calibrated spectral matching tolerance. Changing
#'   it recalibrates the existing feature stack. See \code{\link{rs_params}}.
#' @param time_independent_windows Whether to fill each selected composite
#'   window separately and union the masks (\code{TRUE}, the default), or run
#'   one joint fill across the selected windows (\code{FALSE}).
#' @param search_windows Composite feature-stack windows used for field boundary
#'   delineation: \code{"all"} or a unique integer vector.
#' @param min_windows_alive In joint mode, the minimum number of selected valid
#'   windows in which a pixel must match the seed signature. It is not used by
#'   the default independent-window fills.
#' @param perimeter_margins,trace_iter Cleanup thresholds and maximum
#'   iterations. A cleanup staged by an upstream \code{trace_perimeter()} call
#'   uses the settings given to \code{trace_perimeter()}, and these arguments
#'   do not alter it. Cleanup is applied to the filled or unioned mask.
#' @param merge_separate_areas Whether to attempt to bridge disconnected
#'   candidate components before single-area trimming.
#' @param ... In standalone mode, forwards remaining arguments to
#'   \code{get_rs()}. In piped mode, accepts further named \code{rs_params()}
#'   overrides for this call. Unknown names are rejected.
#' @returns An \code{ocular} object with the candidate field mask stored in
#'   \code{state$alive_mat}.
#' @export
segment_area <- function(x,
                         area_sensitivity = NULL,
                         area_threshold = NULL,
                         strictness = NULL,
                         local_alive_density = NULL,
                         time_independent_windows = NULL,
                         search_windows = NULL,
                         min_windows_alive = NULL,
                         perimeter_margins = NULL,
                         trace_iter = NULL,
                         merge_separate_areas = NULL,
                         params = NULL,
                         ...){

  formal_overrides <- list(
    area_sensitivity = area_sensitivity,
    area_threshold = area_threshold,
    local_alive_density = local_alive_density,
    time_independent_windows = time_independent_windows,
    search_windows = search_windows,
    min_windows_alive = min_windows_alive,
    perimeter_margins = perimeter_margins,
    trace_iter = trace_iter,
    merge_separate_areas = merge_separate_areas,
    strictness = strictness
  )
  input_params <- params
  input_overrides <- c(formal_overrides, list(...))
  if( is_rs(x) ){
    rs <- x
    # Resolve user settings before applying any stale calibration.
    intent <- .resolveParams(rs, input_params, input_overrides,
                             use_calibration = FALSE, caller = "segment_area()")
  }else{
    rs <- get_rs(longitude = x, params = params, ...)
    input_params <- NULL
    input_overrides <- formal_overrides
    intent <- .resolveParams(rs, NULL, input_overrides,
                             use_calibration = FALSE, caller = "segment_area()")
  }
  if( is.null(rs$internals$feature_stack) ){
    rs <- .fetchFeatureStack(rs)
  }
  rs <- .applyScenePriors(rs, intent)
  # Use the refreshed calibration now, while retaining explicit user settings.
  params <- .resolveParams(rs, input_params, input_overrides, caller = "segment_area()")
  # Keep stage calls quiet except when no field is found.
  log_msg <- function(...){

    return(invisible())
  }

  fstack <- rs$internals$feature_stack
  nr <- rs$geom$nr
  nc <- rs$geom$nc

  nw_total <- terra::nlyr(fstack)
  window_idx <- if( identical(params$search_windows, "all") ){
    seq_len(nw_total)
  }else if( is.numeric(params$search_windows) ){
    w <- as.integer(params$search_windows)
    if( any(w < 1L) || any(w > nw_total) || anyDuplicated(w) ){
      stop("segment_area(): search_windows out of range: ", paste(w, collapse = ", "),
           ". Available windows: 1 to ", nw_total, " (", nw_total, " total).",
           call. = FALSE)
    }
    w
  }else{
    stop("search_windows must be \"all\" or integer vector.", call. = FALSE)
  }
  active_rast <- fstack[[window_idx]]
  nw <- length(window_idx)

  feat_array <- array(NA_real_, dim = c(nr, nc, nw))
  for( wi in seq_len(nw) ){
    feat_array[, , wi] <- terra::as.matrix(active_rast[[wi]], wide = TRUE)
  }
  mu <- .learnSignatures(feat_array, rs$geom$centre_rc,
                         params$area_sensitivity, params$area_threshold)

  # Use calibrated signatures and tolerances only for an eligible default
  # sensitivity. Map selected windows to the full calibration and omit invalid
  # ones. Otherwise use the configured scalar sensitivity. Retain features and
  # signatures for later stages.
  cal <- rs$internals$calibration
  use_cal <- FALSE
  if( !is.null(cal) && is.null(cal$error) &&
      identical(cal$schema_version, .CALIBRATION_SCHEMA_VERSION) &&
      !is.null(cal$baselines$area_sensitivity) &&
      is.numeric(cal$tolerance_window) &&
      length(cal$tolerance_window) == nw_total &&
      .paramIsAtDefault(rs, "area_sensitivity") &&
      !"area_sensitivity" %in%
        (attr(params, "user_set", exact = TRUE) %||% character(0L)) &&
      isTRUE(all.equal(unname(params$area_sensitivity),
                       unname(cal$baselines$area_sensitivity))) ){
    sel_valid <- as.logical(cal$valid_window[window_idx])
    if( any(sel_valid) ){
      use_cal <- TRUE
      cal_mu <- cal$mu_window[window_idx]
      cal_tol <- cal$tolerance_window[window_idx]
      missing_mu <- which(!is.finite(mu) & sel_valid & is.finite(cal_mu))
      mu[missing_mu] <- cal_mu[missing_mu]
      log_msg(sprintf("Calibrated fill: %d/%d windows valid. Excluded: %s.",
                      sum(sel_valid), nw,
                      if( all(sel_valid) ){
                        "none"
                      }else{
                        paste(window_idx[!sel_valid], collapse = ",")
                      }))
    }
  }

  if( !is.null(rs$state$alive_mat) ){
    attr(rs, "pending") <- list(pre_cleanup = NULL)
  }

  # A trace_perimeter() call made before filling stores its requested settings
  # for this cleanup.
  pre_cleanup <- isTRUE(attr(rs, "pending")$pre_cleanup$active) &&
    !is.null(attr(rs, "pending")$pre_cleanup$thresholds)
  pre_thresholds <- attr(rs, "pending")$pre_cleanup$thresholds
  pre_max_iter <- attr(rs, "pending")$pre_cleanup$iter
  fill_mode <- if( isTRUE(params$time_independent_windows) ) { "independent" } else { "joint" }
  log_msg(sprintf("Fill mode: %s (W=%d)", fill_mode, nw))

  # The field prior uses its full polygon and an age discount. An unusable
  # polygon supplies no penalty.
  ftw_prior_area <- NULL
  ftw_prior_strength <- 0
  if( !is.null(rs$geom$ftw_prior$polygon) ){
    ftw_prior_area <- .ftwPolygonMask(rs, erode_px = 0L)
    ftw_prior_strength <- .ftwPriorStrength(rs)
    if( is.null(ftw_prior_area) || !isTRUE(ftw_prior_strength > 0) ){
      ftw_prior_area <- NULL
      ftw_prior_strength <- 0
    }else{
      log_msg(sprintf("FTW soft prior: %d in-area px, strength %.3f",
                      sum(ftw_prior_area), ftw_prior_strength))
    }
  }

  input_mask <- rs$state$alive_mat
  rs$state$decided_mat <- matrix(FALSE, nr, nc)
  search_area <- function(...){

    args <- list(...)
    if( is.null(input_mask) ){
      return(do.call(.searchAreaWithWidens, args))
    }
    return(.refillCurrentMask(input_mask, args))
  }
  alive_mat <- matrix(FALSE, nr, nc)

  if( isTRUE(params$time_independent_windows) ){
    for( wi in seq_len(nw) ){
      if( use_cal && !sel_valid[wi] ){
        log_msg(sprintf("--- Window %d (abstained, skipped) ---", window_idx[wi]))
        next
      }
      log_msg(sprintf("--- Window %d ---", window_idx[wi]))
      feat_arr_w <- array(feat_array[, , wi], dim = c(nr, nc, 1L))
      if( use_cal ){
        mu_w <- cal_mu[wi]
        tol_in <- cal_tol[wi]
      }else{
        mu_w <- .learnSignatures(feat_arr_w, rs$geom$centre_rc,
                                 params$area_sensitivity, params$area_threshold)
        tol_in <- params$area_sensitivity
      }
      win_mat <- search_area(feat_array = feat_arr_w,
                             centre_rc = rs$geom$centre_rc,
                             mu = mu_w,
                             area_sensitivity = tol_in,
                             area_threshold = params$area_threshold,
                             min_windows_alive = 1L,
                             local_alive_density = params$local_alive_density,
                             nr_full = nr,
                             nc_full = nc,
                             log_msg = log_msg,
                             prior_field = ftw_prior_area,
                             prior_strength = ftw_prior_strength,
                             prior_cross = .FTW_PRIOR_CROSS)
      alive_mat <- alive_mat | win_mat
      log_msg(sprintf("Window %d: %d alive after merge into union", window_idx[wi], sum(alive_mat)))
    }
  }else{
    if( use_cal ){
      keep <- which(sel_valid)
      feat_jt <- feat_array[, , keep, drop = FALSE]
      mu_jt <- cal_mu[keep]
      tol_jt <- cal_tol[keep]
      if( params$min_windows_alive > length(keep) ){
        stop("min_windows_alive exceeds the number of valid selected windows ",
             "after calibration abstention (", length(keep), ").",
             call. = FALSE)
      }
      min_wa <- as.integer(params$min_windows_alive)
      log_msg(sprintf("Joint calibrated mu (%d window(s)): %s",
                      length(keep), paste(round(mu_jt, 4), collapse = ", ")))
      alive_mat <- search_area(feat_array = feat_jt,
                               centre_rc = rs$geom$centre_rc,
                               mu = mu_jt,
                               area_sensitivity = tol_jt,
                               area_threshold = params$area_threshold,
                               min_windows_alive = min_wa,
                               local_alive_density = params$local_alive_density,
                               nr_full = nr,
                               nc_full = nc,
                               log_msg = log_msg,
                               prior_field = ftw_prior_area,
                               prior_strength = ftw_prior_strength,
                               prior_cross = .FTW_PRIOR_CROSS)
    }else{
      log_msg(sprintf("Joint mu: %s", paste(round(mu, 4), collapse = ", ")))
      if( params$min_windows_alive > nw ){
        stop("min_windows_alive exceeds the number of selected windows (",
             nw, ").", call. = FALSE)
      }
      min_wa <- as.integer(params$min_windows_alive)
      alive_mat <- search_area(feat_array = feat_array,
                               centre_rc = rs$geom$centre_rc,
                               mu = mu,
                               area_sensitivity = params$area_sensitivity,
                               area_threshold = params$area_threshold,
                               min_windows_alive = min_wa,
                               local_alive_density = params$local_alive_density,
                               nr_full = nr,
                               nc_full = nc,
                               log_msg = log_msg,
                               prior_field = ftw_prior_area,
                               prior_strength = ftw_prior_strength,
                               prior_cross = .FTW_PRIOR_CROSS)
    }
  }

  # Clean the combined fill once, using the staged request.
  if( pre_cleanup && any(alive_mat) ){
    rs_pre <- .newOcular()
    rs_pre$geom$nr <- nr
    rs_pre$geom$nc <- nc
    rs_pre$geom$centre_rc <- rs$geom$centre_rc
    rs_pre$state$alive_mat <- alive_mat
    rs_pre$state$decided_mat <- matrix(FALSE, nr, nc)
    p_pre <- list(perimeter_margins = pre_thresholds,
                  trace_iter = pre_max_iter,
                  multiple_areas = params$multiple_areas,
                  area_separation_strict = params$area_separation_strict)
    rs_pre <- .runCleanupOnRs(rs_pre, p_pre, log_msg)
    alive_mat <- rs_pre$state$alive_mat
  }

  if( !any(alive_mat) ){
    rs$state$alive_mat <- alive_mat
    rs$internals$feat_array <- feat_array
    rs$internals$mu <- mu
    attr(rs, "pending")$pre_cleanup$active <- FALSE
    attr(rs, "pending")$pre_cleanup$thresholds <- NULL
    attr(rs, "pending")$pre_cleanup$iter <- 1L
    warning("segment_area(): no candidate field was identified. The returned ",
            "mask is empty.", call. = FALSE)
    return(rs)
  }

  # Bridging is independent of the number of output areas.
  decided_mat <- rs$state$decided_mat
  if( isTRUE(params$merge_separate_areas) ){
    info <- .enumerateAreas(alive_mat, rs$geom$centre_rc, strict = FALSE)
    mr <- .mergeSeparableAreas(alive_mat, decided_mat, info,
                               feat_array, mu,
                               params$area_sensitivity, params$area_threshold)
    alive_mat <- mr$alive_mat
    decided_mat <- mr$decided_mat
    if( mr$revived > 0L ){
      log_msg(sprintf("Merge: %d bridge(s) revived", mr$revived))
    }
  }

  # Select one area after bridging, using the requested connectivity.
  if( isFALSE(params$multiple_areas) ){
    info <- .applyStrictSeparation(alive_mat, rs$geom$centre_rc,
                                   params$area_separation_strict)
    if( info$n_areas > 1L ){
      keep_label <- info$area_order[1L]
      alive_mat <- info$label_mat == keep_label
      log_msg(sprintf("Single-area: keeping rank-1 area (%d px)", sum(alive_mat)))
    }
  }

  rs$state$alive_mat <- alive_mat
  rs$state$decided_mat <- decided_mat
  rs$internals$feat_array <- feat_array
  rs$internals$mu <- mu
  attr(rs, "pending")$pre_cleanup$active <- FALSE
  attr(rs, "pending")$pre_cleanup$thresholds <- NULL
  attr(rs, "pending")$pre_cleanup$iter <- 1L
  log_msg(sprintf("Detected %d alive px", sum(alive_mat)))
  return(rs)
}
