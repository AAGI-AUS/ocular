# =========================================================================
# Segmentation: calibration
# =========================================================================
## Measures a robust field signature from either the seed neighbourhood or
## an attached FTW polygon. Per-window tolerances, validity checks,
## multimodality abstention, and calibration provenance are computed here;
## the flood-fill implementation consumes the resulting record elsewhere.

.C_MAD <- 1.166
## Asymptotic SE coefficient of the consistency-scaled MAD under the
## Gaussian reference (Rousseeuw & Croux 1993): SE(disp) ~= .C_MAD *
## disp / sqrt(n). This is a literature-derived constant.

.REPRESENTATIVE_RHO <- c(nir = 0.40, red = 0.05)
## Fixed representative vegetated reflectance for c_noise. This
## healthy-vegetation operating point (NIR plateau ~0.40, red absorption
## ~0.05), with sigma_rho = 0.005, gives NDVI c_noise ~0.018-0.020.

.CALIBRATION_SCHEMA_VERSION <- 2L
## Increment when a calibration record's measurement support or interpretation
## changes. Records without the current marker are recomputed before use.

#' Noise-equivalent VI difference for an index
#'
#' First-order propagation of per-band reflectance noise through the
#' index formula at the representative operating point. sigma_rho is a
#' single engineering-order radiometric noise value. Returns \code{NA_real_}
#' for indices without an implemented propagation; callers then omit the
#' noise-derived tolerance floor.
#'
#' @noRd
.cNoise <- function(index_name, sigma_rho = 0.005,
                    rho = .REPRESENTATIVE_RHO){

  N <- rho[["nir"]]; R <- rho[["red"]]
  if( identical(index_name, "NDVI") ){
    ## NDVI = (N - R)/(N + R):
    ## c = (2/(N+R)^2) * sqrt(R^2 sig_N^2 + N^2 sig_R^2)
    return((2 / (N + R)^2) * sqrt(R^2 * sigma_rho^2 + N^2 * sigma_rho^2))
  }
  if( identical(index_name, "EVI2") ){
    ## EVI2 = 2.5 (N - R)/(N + 2.4 R + 1); partials:
    ## d/dN = 2.5 (3.4 R + 1)/D^2 ; d/dR = -2.5 (3.4 N + 1)/D^2
    D  <- N + 2.4 * R + 1
    dN <- 2.5 * (3.4 * R + 1) / D^2
    dR <- -2.5 * (3.4 * N + 1) / D^2
    return(sqrt(dN^2 * sigma_rho^2 + dR^2 * sigma_rho^2))
  }
  return(NA_real_)
}

.MULTIMODAL_SEP_K <- 3.0
## Engineering threshold for the robust max-central-gap statistic returned by
## .modeSeparation. This value is not fitted to observed scenes and requires
## external validation.

.MULTIMODAL_MIN_N <- 8L
## Minimum per-window support for applying the modality check. Smaller samples
## return no modality decision.

#' Robust max-central-gap mode-separation statistic
#'
#' Detects separated clusters that may not increase the median absolute
#' deviation when one cluster contains most observations. The statistic sorts
#' the support, trims 5 percent from each tail, and divides the largest gap
#' between adjacent central values by the larger of the measured dispersion
#' and the noise floor. The comparison is therefore relative to the
#' within-cluster scale.
#'
#' @returns the gap ratio (>= 0), or NA when support is too small
#'   (< .MULTIMODAL_MIN_N) or the noise scale is unavailable.
#' @noRd
.modeSeparation <- function(v, disp, c_noise){
  v <- sort(v[is.finite(v)])
  n <- length(v)
  if( n < .MULTIMODAL_MIN_N || !is.finite(c_noise) )
    return(NA_real_)
  k_trim  <- max(1L, as.integer(floor(0.05 * n)))
  central <- v[(k_trim + 1L):(n - k_trim)]
  if( length(central) < 2L ) return(0)
  max_gap <- max(diff(central))
  scale   <- max(if( is.na(disp) ) 0 else disp, c_noise)
  return(max_gap / scale)
}

#' Largest complete seed-centred square available on the feature grid
#'
#' Returns the Chebyshev radius, in pixels, of the largest complete square
#' centred on \code{rs$geom$centre_rc} that fits inside the feature grid.
#' Invalid, incomplete, or grid-inconsistent geometry returns \code{NULL}.
#'
#' @noRd
.seedSupportRadiusLimit <- function(rs){
  centre <- rs$geom$centre_rc
  geom_dims <- c(rs$geom$nr, rs$geom$nc)
  valid_pair <- function(x){
    is.numeric(x) && length(x) == 2L && !anyNA(x) &&
      all(is.finite(x)) && all(x == trunc(x))
  }
  if( !valid_pair(centre) || !valid_pair(geom_dims) ||
      any(geom_dims < 1L) || any(centre < 1L) ||
      any(centre > geom_dims) )
    return(NULL)

  fstack <- rs$internals$feature_stack
  if( is.null(fstack) ) return(NULL)
  grid_dims <- tryCatch({
    if( terra::nlyr(fstack) < 1L ) NULL
    else c(terra::nrow(fstack), terra::ncol(fstack))
  }, error = function(e) NULL)
  if( !valid_pair(grid_dims) ||
      !identical(as.integer(grid_dims), as.integer(geom_dims)) )
    return(NULL)

  centre <- as.integer(centre)
  grid_dims <- as.integer(grid_dims)
  as.integer(min(centre[1L] - 1L, grid_dims[1L] - centre[1L],
                 centre[2L] - 1L, grid_dims[2L] - centre[2L]))
}

#' Seed-neighbourhood measurement support with a precision stopping rule
#'
#' Uses a complete seed-centred square, beginning at a two-pixel radius when
#' the grid permits and otherwise using the largest available square. The
#' support widens until every window's dispersion estimate is
#' noise-floor-precise
#' (SE(disp_w) = .C_MAD * disp_w / sqrt(n_w) <= c_noise), the feature-grid
#' limit is reached, or expansion adds no cells. Dispersion is re-estimated
#' after each increase in radius. A zero dispersion passes the stopping rule;
#' a missing estimate does not. When c_noise is unavailable, the function
#' returns the initial support and records no precision decision.
#'
#' @returns list(sig, r, precise_w, support) or NULL when no feature
#'   data / centre is available.
#' @noRd
.liteSupport <- function(rs, c_noise){

  r_max_px <- .seedSupportRadiusLimit(rs)
  if( is.null(r_max_px) ) return(NULL)
  initial_r <- min(2L, r_max_px)
  r <- initial_r
  prior_n <- -1L
  stop_reason <- NULL
  repeat{
    nb <- .seedNeighborhood(rs, radius = r)
    if( is.null(nb) ) return(NULL)
    sig  <- .signatureFromMatrix(nb)
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
    if( r >= r_max_px ){
      stop_reason <- "feature_grid_limit"
      break
    }
    r <- r + 1L
  }
  return(list(sig       = sig,
              nb        = nb,
              r         = r,
              precise_w = if( is.na(c_noise) ) rep(NA, length(sig$n_w))
              else pass,
              support   = list(kind = "seed_neighbourhood",
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

#' Calibrate detection-window signatures and tolerances
#'
#' Measures a threshold-free signature on the selected support and derives
#' per-window tolerances as \eqn{max(k * dispersion, c_noise)}. A window is
#' omitted when support precision is insufficient, dispersion crosses the
#' admissible range, or the modality check detects separated clusters. The
#' scalar fallback is the largest tolerance among valid windows. The record
#' includes the per-window values used by \code{segment_area()} when the user
#' has not set \code{area_sensitivity}. Perimeter, split, and interior controls
#' are not derived here.
#'
#' @returns The calibration record; failures return a reduced record containing
#'   \code{error} and \code{warnings}.
#' @noRd
.calibrateScene <- function(rs){

  warnings <- character(0)
  params     <- rs$params
  strictness <- params$strictness %||% "balanced"
  k          <- .strictnessK(strictness)
  ## Calibration measures the detection feature stack, which can deliberately
  ## use a different index from the user-facing output scenes.
  index_name <- params$search_index %||% rs$spec$index_name
  c_noise    <- .cNoise(index_name)
  if( is.na(c_noise) )
    warnings <- c(warnings, sprintf(
      paste0("Noise propagation is unavailable for index '%s'; calibration ",
             "did not derive a noise floor or area_sensitivity baseline."),
      index_name))

  ## A supplied FTW polygon measures the signature on its interior; otherwise
  ## calibration uses the seed neighbourhood. Both paths return the same
  ## measurement structure.
  if( !is.null(rs$geom$ftw_prior$polygon) ){
    ls         <- .ftwSupport(rs, c_noise = c_noise)
    cal_source <- "ftw"
  }else{
    ls         <- .liteSupport(rs, c_noise = c_noise)
    cal_source <- "self"
  }
  if( is.null(ls) )
    return(list(schema_version = .CALIBRATION_SCHEMA_VERSION,
                source   = cal_source,
                error    = "no support (feature data unavailable)",
                warnings = c(warnings,
                             paste0(cal_source, " support unavailable: ",
                                    "feature data, centre, or polygon missing"))))
  sig <- ls$sig
  tv  <- .toleranceFromDispersion(
    disp_w       = sig$disp_w,
    n_w          = sig$n_w,
    k            = k,
    c_noise      = if( is.na(c_noise) ) 0 else c_noise,
    n_min        = 2,
    area_threshold = params$area_threshold %||% c(0, 1))

  valid_w <- tv$valid_w
  reason_w <- tv$reason_w
  if( !is.na(c_noise) ){
    imprecise <- valid_w & !ls$precise_w
    if( any(imprecise) ){
      valid_w[imprecise] <- FALSE
      reason_w[imprecise] <- "imprecise_support"
      support_detail <- if (identical(cal_source, "self"))
        sprintf("seed support (r = %d, grid limit %d)",
                ls$r, ls$support$r_max)
      else "FTW polygon support"
      warnings <- c(warnings, sprintf(
        paste0("window(s) %s: dispersion was not estimated precisely on %s; ",
               "excluded from calibration."),
        paste(which(imprecise), collapse = ","), support_detail))
    }
  }
  cross <- tv$reason_w == "disp_crossover"
  if( any(cross) )
    warnings <- c(warnings, sprintf(
      paste0("window(s) %s: dispersion met or exceeded the admissible ",
             "limit (%.4g); excluded from calibration."),
      paste(which(cross), collapse = ","), tv$disp_star))

  ## The dispersion rule may not detect a smaller cluster separated from a
  ## dominant cluster. Apply the mode-separation statistic before accepting
  ## the window, while keeping the signature estimator limited to median and
  ## median absolute deviation.
  sep_w <- rep(NA_real_, length(sig$mu_w))
  if( !is.na(c_noise) && !is.null(ls$nb) ){
    for( w in seq_along(sep_w) )
      sep_w[w] <- .modeSeparation(ls$nb[, w], sig$disp_w[w], c_noise)
  }
  multimodal_w <- !is.na(sep_w) & sep_w >= .MULTIMODAL_SEP_K
  if( any(multimodal_w & valid_w) ){
    hit <- which(multimodal_w & valid_w)
    valid_w[hit] <- FALSE
    reason_w[hit] <- "multimodal_support"
    warnings <- c(warnings, sprintf(
      paste0("window(s) %s: mode separation met or exceeded %.2g; ",
             "excluded from calibration."),
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

  return(list(schema_version   = .CALIBRATION_SCHEMA_VERSION,
              source           = cal_source,
              strictness       = strictness,
              area_threshold  = params$area_threshold,
              k                = k,
              index_name       = index_name,
              c_noise          = c_noise,
              mu_window        = sig$mu_w,
              disp_window      = sig$disp_w,
              n_window         = sig$n_w,
              tolerance_window = tv$tolerance_w,
              valid_window     = valid_w,
              reason_window    = reason_w,
              multimodal_window     = multimodal_w,
              mode_separation_window = sep_w,
              disp_star        = tv$disp_star,
              support          = ls$support,
              baselines        = baselines,
              warnings         = warnings))
}


# =========================================================================
# Shared calibration helpers
# =========================================================================
# These small helpers are shared with the segmentation layer:
#   .paramIsAtDefault  - used by segment_area() (segmentation)
#   .seedNeighborhood  - used by .liteSupport()  (calibration)

#' Is an rs parameter still at its rs_params() default? (internal)
#'
#' Explicit ownership metadata takes precedence, including when a user chose
#' the canonical default value. A direct non-default mutation is also
#' protected.
#'
#' @param rs An rs object.
#' @param area Character scalar naming an \code{rs_params()} area.
#' @noRd
.paramIsAtDefault <- function(rs, area){
  defaults <- rs_params()
  if( !is.character(area) || length(area) != 1L || is.na(area) ||
      !area %in% names(defaults) )
    stop(".paramIsAtDefault(): unknown parameter.", call. = FALSE)
  stored <- rs$params
  if( is.null(stored) || length(stored) == 0L ) return(TRUE)
  user_set <- attr(stored, "user_set", exact = TRUE)
  if( is.null(user_set) )
    stop(".paramIsAtDefault(): rs$params is missing explicit parameter ",
         "ownership metadata. Recreate the ocular object with get_rs().",
         call. = FALSE)
  user_set <- .paramUserSet(stored, ".paramIsAtDefault()")
  if( is.null(names(stored)) || !area %in% names(stored) ) return(TRUE)
  if( area %in% user_set ) return(FALSE)
  identical(stored[[area]], defaults[[area]])
}

#' Extract the feature matrix for the seed neighborhood (internal)
#'
#' Returns an \code{(npix, W)} numeric matrix of per-window feature values
#' for the (2*radius + 1) x (2*radius + 1) block centred on the seed pixel
#' (\code{rs$geom$centre_rc}). Reads every layer of the full
#' \code{feature_stack}. \code{NULL} when no centre or feature data is
#' available.
#'
#' @param rs An rs object.
#' @param radius Neighborhood half-width in pixels (default 1L).
#' @noRd
.seedNeighborhood <- function(rs, radius = 1L){
  centre <- rs$geom$centre_rc
  if( is.null(centre) ) return(NULL)
  r0 <- centre[1L]; c0 <- centre[2L]
  nr <- rs$geom$nr; nc <- rs$geom$nc
  r_lo <- max(1L, r0 - radius); r_hi <- min(nr, r0 + radius)
  c_lo <- max(1L, c0 - radius); c_hi <- min(nc, c0 + radius)

  ## Calibration must always see the complete detection series. A
  ## segment_area(search_windows = ...) call caches only the selected layers in
  ## feat_array, so preferring that cache would make later recalibration stale.
  fstack <- rs$internals$feature_stack
  if( is.null(fstack) ) return(NULL)
  nw <- terra::nlyr(fstack)
  out <- matrix(NA_real_, nrow = (r_hi - r_lo + 1L) * (c_hi - c_lo + 1L),
                ncol = nw)
  for( wi in seq_len(nw) ){
    layer_mat <- terra::as.matrix(fstack[[wi]], wide = TRUE)
    block <- layer_mat[r_lo:r_hi, c_lo:c_hi, drop = FALSE]
    out[, wi] <- as.vector(block)
  }
  out
}
