# =========================================================================
# Segmentation: calibration
# =========================================================================
# Measure signatures and tolerances from the seed neighbourhood or a field
# polygon. Flood fill uses the resulting record.

.C_MAD <- 1.166
# Standard error coefficient for the scaled MAD under a Gaussian reference
# (Rousseeuw & Croux 1993): SE(disp) ~= .C_MAD * disp / sqrt(n).

.REPRESENTATIVE_RHO <- c(nir = 0.40, red = 0.05)
# Representative vegetation reflectance for noise propagation: NIR 0.40,
# red 0.05. With sigma_rho = 0.005, NDVI noise is about 0.018–0.020.

.CALIBRATION_SCHEMA_VERSION <- 3L
# Increment the schema when support or its interpretation changes.
# Older records must be recalibrated before use.

#' Propagate reflectance noise through an index formula
#'
#' Use first derivatives at the representative reflectance values.
#' \code{sigma_rho} is a fixed engineering estimate of band noise.
#' Return \code{NA_real_} for unsupported indices, which then have no noise floor.
#' @noRd
.cNoise <- function(index_name, sigma_rho = 0.005,
                    rho = .REPRESENTATIVE_RHO){

  N <- rho[["nir"]]
  R <- rho[["red"]]
  if( identical(index_name, "NDVI") ){
    # Propagate independent NIR and red noise through NDVI.
    return((2 / (N + R)^2) * sqrt(R^2 * sigma_rho^2 + N^2 * sigma_rho^2))
  }
  if( identical(index_name, "EVI2") ){
    # Partial derivatives of EVI2 with respect to NIR and red reflectance.
    D <- N + 2.4 * R + 1
    dN <- 2.5 * (3.4 * R + 1) / D^2
    dR <- -2.5 * (3.4 * N + 1) / D^2
    return(sqrt(dN^2 * sigma_rho^2 + dR^2 * sigma_rho^2))
  }
  return(NA_real_)
}

.MULTIMODAL_SEP_K <- 3.0
# The gap threshold is an engineering default, not fitted or externally validated.

.MULTIMODAL_MIN_N <- 8L
# Smaller samples do not support a decision about separate modes.

#' Detect separated clusters that MAD can miss
#'
#' A dominant cluster can keep MAD small despite a separated minority. Trim
#' five percent from each tail, removing at least one value from each. Scale
#' the largest remaining gap by the larger of dispersion and the noise floor.
#'
#' @returns A nonnegative gap ratio, or \code{NA} when support is too small
#'   or the noise scale is unavailable.
#' @noRd
.modeSeparation <- function(v, disp, c_noise){

  v <- sort(v[is.finite(v)])
  n <- length(v)
  if( n < .MULTIMODAL_MIN_N || !is.finite(c_noise) ){
    return(NA_real_)
  }
  k_trim <- max(1L, as.integer(floor(0.05 * n)))
  central <- v[(k_trim + 1L):(n - k_trim)]
  if( length(central) < 2L ){
    return(0)
  }
  max_gap <- max(diff(central))
  scale <- max(if( is.na(disp) ) { 0 } else { disp }, c_noise)
  return(max_gap / scale)
}

#' Find the largest complete square centred on the seed
#'
#' Return its radius in pixels, or \code{NULL} for missing or inconsistent geometry.
#' @noRd
.seedSupportRadiusLimit <- function(rs){

  centre <- rs$geom$centre_rc
  geom_dims <- c(rs$geom$nr, rs$geom$nc)
  valid_pair <- function(x){

    return(is.numeric(x) && length(x) == 2L && !anyNA(x) &&
             all(is.finite(x)) && all(x == trunc(x)))
  }
  if( !valid_pair(centre) || !valid_pair(geom_dims) ||
      any(geom_dims < 1L) || any(centre < 1L) ||
      any(centre > geom_dims) ){
    return(NULL)
  }

  fstack <- rs$internals$feature_stack
  if( is.null(fstack) ){
    return(NULL)
  }
  grid_dims <- tryCatch({
    if( terra::nlyr(fstack) < 1L ) { NULL } else { c(terra::nrow(fstack), terra::ncol(fstack)) }
  }, error = function(e){
    return(NULL)
  })
  if( !valid_pair(grid_dims) ||
      !identical(as.integer(grid_dims), as.integer(geom_dims)) ){
    return(NULL)
  }

  centre <- as.integer(centre)
  grid_dims <- as.integer(grid_dims)
  return(as.integer(min(centre[1L] - 1L, grid_dims[1L] - centre[1L],
                        centre[2L] - 1L, grid_dims[2L] - centre[2L])))
}

#' Expand seed support until dispersion is measured precisely
#'
#' Start with a complete square of radius two pixels, or the largest that fits.
#' Widen and remeasure until every measurable window's standard error is at most
#' the noise floor, the grid limit is reached, or no cells are added. Zero
#' dispersion passes this test. A missing estimate does not.
#'
#' Exclude windows with fewer than two finite cells in the largest square.
#' Keep the initial support if none can be measured or the noise floor is unknown.
#'
#' @returns A list with \code{sig}, \code{r}, \code{precise_w} and \code{support},
#'   or \code{NULL} when features or the centre are unavailable.
#' @noRd
.liteSupport <- function(rs, c_noise){

  r_max_px <- .seedSupportRadiusLimit(rs)
  if( is.null(r_max_px) ){
    return(NULL)
  }
  initial_r <- min(2L, r_max_px)
  r <- initial_r
  prior_n <- -1L
  stop_reason <- NULL
  cannot_pass <- NULL
  repeat{
    nb <- .seedNeighborhood(rs, radius = r)
    if( is.null(nb) ){
      return(NULL)
    }
    sig <- .signatureFromMatrix(nb)
    se_w <- .C_MAD * sig$disp_w / sqrt(sig$n_w)
    pass <- !is.na(se_w) & se_w <= c_noise
    if( nrow(nb) <= prior_n ){
      stop_reason <- "no_new_cells"
      break
    }
    prior_n <- nrow(nb)
    if( is.na(c_noise) ){
      stop_reason <- "noise_floor_unavailable"
      break
    }
    if( all(pass) ){
      stop_reason <- "precision_reached"
      break
    }
    if( is.null(cannot_pass) ){
      largest_nb <- if( r == r_max_px ) { nb } else { .seedNeighborhood(rs, radius = r_max_px) }
      if( is.null(largest_nb) ){
        return(NULL)
      }
      cannot_pass <- colSums(is.finite(largest_nb)) < 2L
    }
    if( all(pass | cannot_pass) ){
      stop_reason <- if( all(cannot_pass) ) { "no_measurable_windows" } else { "precision_reached" }
      break
    }
    if( r >= r_max_px ){
      stop_reason <- "feature_grid_limit"
      break
    }
    r <- r + 1L
  }
  return(list(sig = sig,
              nb = nb,
              r = r,
              precise_w = if( is.na(c_noise) ) { rep(NA, length(sig$n_w)) } else { pass },
              support = list(kind = "seed_neighbourhood",
                               r = r,
                               r_max = r_max_px,
                               centre_row = as.integer(rs$geom$centre_rc[1L]),
                               centre_col = as.integer(rs$geom$centre_rc[2L]),
                               initial_radius = initial_r,
                               width_px = 2L * r + 1L,
                               limit_width_px = 2L * r_max_px + 1L,
                               grid_nrow = as.integer(rs$geom$nr),
                               grid_ncol = as.integer(rs$geom$nc),
                               geometric_cell_count = (2L * r + 1L)^2L,
                               limit_rule =
                                 "largest_complete_seed_centred_square",
                               calibration_schema_version =
                                 .CALIBRATION_SCHEMA_VERSION,
                               stop_reason = stop_reason,
                               n_w = sig$n_w)))
}

#' Measure signatures and set tolerances for detection windows
#'
#' Use the larger of scaled dispersion and the noise floor as tolerance.
#' Exclude windows with imprecise support, excessive dispersion or separate
#' modes. The largest valid tolerance supplies a scalar baseline. These values
#' apply only when the user has not chosen \code{area_sensitivity}.
#' Perimeter, split and interior controls are unchanged.
#'
#' @returns The calibration record, or a failure record with \code{error}
#'   and \code{warnings}.
#' @noRd
.calibrateScene <- function(rs){

  warnings <- character(0)
  params <- rs$params
  strictness <- params$strictness %||% "balanced"
  k <- .strictnessK(strictness)
  # Detection features can use a different index from the output scenes.
  index_name <- params$search_index %||% rs$spec$index_name
  c_noise <- .cNoise(index_name)
  if( is.na(c_noise) ){
    warnings <- c(warnings, sprintf(
      paste0("Noise propagation is unavailable for index '%s'. Calibration ",
             "did not derive a noise floor or area_sensitivity baseline."),
      index_name))
  }

  if( !is.null(rs$geom$ftw_prior$polygon) ){
    ls <- .ftwSupport(rs, c_noise = c_noise)
    cal_source <- "ftw"
  }else{
    ls <- .liteSupport(rs, c_noise = c_noise)
    cal_source <- "self"
  }
  if( is.null(ls) ){
    return(list(schema_version = .CALIBRATION_SCHEMA_VERSION,
                source = cal_source,
                error = "no support (feature data unavailable)",
                warnings = c(warnings,
                             paste0(cal_source, " support unavailable: ",
                                    "feature data, centre, or polygon missing"))))
  }
  sig <- ls$sig
  tv <- .toleranceFromDispersion(
    disp_w = sig$disp_w,
    n_w = sig$n_w,
    k = k,
    c_noise = if( is.na(c_noise) ) { 0 } else { c_noise },
    n_min = 2,
    area_threshold = params$area_threshold %||% c(0, 1))

  valid_w <- tv$valid_w
  reason_w <- tv$reason_w
  if( !is.na(c_noise) ){
    imprecise <- valid_w & !ls$precise_w
    if( any(imprecise) ){
      valid_w[imprecise] <- FALSE
      reason_w[imprecise] <- "imprecise_support"
      support_detail <- if( identical(cal_source, "self") ){
        sprintf("seed support (r = %d, grid limit %d)",
                ls$r, ls$support$r_max)
      }else{
        "FTW polygon support"
      }
      warnings <- c(warnings, sprintf(
        paste0("window(s) %s: dispersion was not estimated precisely on %s. ",
               "These windows were excluded from calibration."),
        paste(which(imprecise), collapse = ","), support_detail))
    }
  }
  cross <- tv$reason_w == "disp_crossover"
  if( any(cross) ){
    warnings <- c(warnings, sprintf(
      paste0("window(s) %s: dispersion met or exceeded the admissible ",
             "limit (%.4g). These windows were excluded from calibration."),
      paste(which(cross), collapse = ","), tv$disp_star))
  }

  # MAD can miss a separated minority cluster, so check gaps before accepting a window.
  sep_w <- rep(NA_real_, length(sig$mu_w))
  if( !is.na(c_noise) && !is.null(ls$nb) ){
    for( w in seq_along(sep_w) ){
      sep_w[w] <- .modeSeparation(ls$nb[, w], sig$disp_w[w], c_noise)
    }
  }
  multimodal_w <- !is.na(sep_w) & sep_w >= .MULTIMODAL_SEP_K
  if( any(multimodal_w & valid_w) ){
    hit <- which(multimodal_w & valid_w)
    valid_w[hit] <- FALSE
    reason_w[hit] <- "multimodal_support"
    warnings <- c(warnings, sprintf(
      paste0("window(s) %s: mode separation met or exceeded %.2g. ",
             "These windows were excluded from calibration."),
      paste(hit, collapse = ","), .MULTIMODAL_SEP_K))
  }

  baselines <- list()
  if( !is.na(c_noise) ){
    if( any(valid_w) ){
      baselines$area_sensitivity <- max(tv$tolerance_w[valid_w])
    }else{
      warnings <- c(warnings,
                    "no valid detection window: no calibration baseline was derived.")
    }
  }

  return(list(schema_version = .CALIBRATION_SCHEMA_VERSION,
              source = cal_source,
              strictness = strictness,
              area_threshold = params$area_threshold,
              k = k,
              index_name = index_name,
              c_noise = c_noise,
              mu_window = sig$mu_w,
              disp_window = sig$disp_w,
              n_window = sig$n_w,
              tolerance_window = tv$tolerance_w,
              valid_window = valid_w,
              reason_window = reason_w,
              multimodal_window = multimodal_w,
              mode_separation_window = sep_w,
              disp_star = tv$disp_star,
              support = ls$support,
              baselines = baselines,
              warnings = warnings))
}


# =========================================================================
# Shared calibration helpers
# =========================================================================

#' Check whether a parameter is still available for calibration
#'
#' An explicit setting is retained even when it equals the default. A value
#' changed from its default is also retained.
#'
#' @param rs An ocular object.
#' @param area A parameter name accepted by \code{rs_params()}.
#' @noRd
.paramIsAtDefault <- function(rs, area){

  defaults <- rs_params()
  if( !is.character(area) || length(area) != 1L || is.na(area) ||
      !area %in% names(defaults) ){
    stop(".paramIsAtDefault(): unknown parameter.", call. = FALSE)
  }
  stored <- rs$params
  if( is.null(stored) || length(stored) == 0L ){
    return(TRUE)
  }
  user_set <- attr(stored, "user_set", exact = TRUE)
  if( is.null(user_set) ){
    stop(".paramIsAtDefault(): rs$params is missing explicit parameter ",
         "ownership metadata. Recreate the ocular object with get_rs().",
         call. = FALSE)
  }
  user_set <- .paramUserSet(stored, ".paramIsAtDefault()")
  if( is.null(names(stored)) || !area %in% names(stored) ){
    return(TRUE)
  }
  if( area %in% user_set ){
    return(FALSE)
  }
  return(identical(stored[[area]], defaults[[area]]))
}

#' Extract the feature matrix for the seed neighbourhood
#'
#' Returns an \code{(npix, W)} numeric matrix of per-window feature values
#' for the (2*radius + 1) x (2*radius + 1) block centred on the seed pixel
#' (\code{rs$geom$centre_rc}). Reads every layer of the full
#' \code{feature_stack}. \code{NULL} when no centre or feature data is
#' available.
#'
#' @param rs An rs object.
#' @param radius Neighbourhood half-width in pixels (default 1L).
#' @noRd
.seedNeighborhood <- function(rs, radius = 1L){

  centre <- rs$geom$centre_rc
  if( is.null(centre) ){
    return(NULL)
  }
  r0 <- centre[1L]
  c0 <- centre[2L]
  nr <- rs$geom$nr
  nc <- rs$geom$nc
  r_lo <- max(1L, r0 - radius)
  r_hi <- min(nr, r0 + radius)
  c_lo <- max(1L, c0 - radius)
  c_hi <- min(nc, c0 + radius)

  # Use the full detection series. feat_array may hold only windows selected
  # by an earlier segment_area() call.
  fstack <- rs$internals$feature_stack
  if( is.null(fstack) ){
    return(NULL)
  }
  nw <- terra::nlyr(fstack)
  out <- matrix(NA_real_, nrow = (r_hi - r_lo + 1L) * (c_hi - c_lo + 1L),
                ncol = nw)
  for( wi in seq_len(nw) ){
    layer_mat <- terra::as.matrix(fstack[[wi]], wide = TRUE)
    block <- layer_mat[r_lo:r_hi, c_lo:c_hi, drop = FALSE]
    out[, wi] <- as.vector(block)
  }
  return(out)
}
