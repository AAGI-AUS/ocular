# =========================================================================
# Segmentation: threshold-free calibration signatures
# =========================================================================
## Shared by the seed-neighbourhood and FTW calibration paths. Constants
## (c_noise, n_min, k) are supplied as arguments or package constants and are
## not estimated from validation labels.

#' Tolerance and validity from per-window dispersion
#'
#' The tolerance is \eqn{max(k * disp_w, c_noise)}, where c_noise is the
#' noise-equivalent VI difference. A window is valid when its support contains
#' at least n_min observations and its dispersion is below
#' \eqn{(hi - lo)/(2k)}. Reasons report the first failed condition in the order
#' n_below_min, disp_na, and disp_crossover.
#'
#' @param disp_w,n_w Per-window dispersion estimates and finite support counts.
#' @param k Strictness robustness multiple (dimensionless).
#' @param c_noise Noise-equivalent VI difference for the index/sensor
#'   (for example, the NDVI propagation
#'   sigma_NDVI^2 = 4/(N+R)^4 * (R^2 sigma_N^2 + N^2 sigma_R^2)).
#' @param n_min Minimum support size for a valid window.
#' @param area_threshold c(lo, hi) coarse plausibility band.
#' @returns list(tolerance_w, valid_w, reason_w, disp_star).
#' @noRd
.toleranceFromDispersion <- function(disp_w, n_w, k, c_noise, n_min,
                                     area_threshold){

  if( !.isScalarNumber(k) || k <= 0 )
    stop(".toleranceFromDispersion(): k must be a single positive numeric.",
         call. = FALSE)
  if( !.isScalarNumber(c_noise) || c_noise < 0 )
    stop(".toleranceFromDispersion(): c_noise must be a single non-negative ",
         "numeric.", call. = FALSE)
  if( !.isWholeNumber(n_min) || n_min < 1 )
    stop(".toleranceFromDispersion(): n_min must be a positive integer.",
         call. = FALSE)
  if( !is.numeric(area_threshold) || length(area_threshold) != 2L ||
      anyNA(area_threshold) || any(!is.finite(area_threshold)) ||
      area_threshold[1L] > area_threshold[2L] )
    stop(".toleranceFromDispersion(): area_threshold must be c(lo, hi) with ",
         "lo <= hi.", call. = FALSE)

  disp_star   <- (area_threshold[2L] - area_threshold[1L]) / (2 * k)
  tolerance_w <- pmax(k * disp_w, c_noise)

  reason_w <- rep("", length(disp_w))
  reason_w[n_w < n_min] <- "n_below_min"
  reason_w[reason_w == "" & is.na(disp_w)] <- "disp_na"
  cross <- reason_w == "" & !is.na(disp_w) & disp_w >= disp_star
  reason_w[cross] <- "disp_crossover"
  valid_w <- reason_w == ""

  return(list(tolerance_w = tolerance_w,
              valid_w     = valid_w,
              reason_w    = reason_w,
              disp_star   = disp_star))
}


#' Threshold-free per-window signature from a support matrix
#'
#' Estimates the median and scaled median absolute deviation over a support
#' extracted as a (pixels x windows) matrix, for example the output of
#' \code{.seedNeighborhood()}.
#' Windows require at least two valid pixels because one pixel has no measured
#' dispersion.
#'
#' @param vals_mat Numeric matrix (n_pixels, W).
#' @returns list(mu_w, disp_w, n_w), each length W.
#' @noRd
.signatureFromMatrix <- function(vals_mat){

  if( !is.matrix(vals_mat) )
    stop(".signatureFromMatrix(): vals_mat must be a matrix.", call. = FALSE)
  nw     <- ncol(vals_mat)
  mu_w   <- rep(NA_real_, nw)
  disp_w <- rep(NA_real_, nw)
  n_w    <- integer(nw)
  for( wi in seq_len(nw) ){
    v <- vals_mat[, wi]
    v <- v[is.finite(v)]
    n_w[wi] <- length(v)
    if( n_w[wi] >= 2L ){
      mu_w[wi]   <- stats::median(v)
      disp_w[wi] <- 1.4826 * stats::median(abs(v - mu_w[wi]))
    }
  }
  return(list(mu_w = mu_w, disp_w = disp_w, n_w = n_w))
}

## Strictness -> robustness multiple k. These fixed, spectral-only
## engineering defaults preserve the required ordinal ordering, but their
## spacing needs external validation before it is treated as scientifically
## calibrated.
.STRICTNESS_K <- c(loose = 3.0, balanced = 2.2, tight = 1.5)

#' Strictness ordinal to fixed robustness multiple k (spectral-only)
#' @noRd
.strictnessK <- function(strictness){
  if( !is.character(strictness) || length(strictness) != 1L ||
      !strictness %in% names(.STRICTNESS_K) )
    stop(".strictnessK(): strictness must be one of ",
         paste0('"', names(.STRICTNESS_K), '"', collapse = ", "), ".",
         call. = FALSE)
  return(.STRICTNESS_K[[strictness]])
}
