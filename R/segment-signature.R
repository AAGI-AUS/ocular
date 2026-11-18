# =========================================================================
# Segmentation: calibration signatures
# =========================================================================
# Seed and field polygon support use the same measurements. Constants are
# supplied directly, not fitted to validation labels.

#' Derive tolerances and validity from each window's dispersion
#'
#' Tolerance is \eqn{max(k * disp_w, c_noise)}. Valid windows need at least
#' \code{n_min} observations and dispersion below \eqn{(hi - lo)/(2k)}.
#' Report the first failed test: sample size, missing dispersion, then crossover.
#'
#' @param disp_w,n_w Dispersion and finite observation count for each window.
#' @param k Robustness multiple selected by strictness.
#' @param c_noise Index difference expected from reflectance noise.
#' @param n_min Minimum support size.
#' @param area_threshold Lower and upper plausible index values.
#' @returns A list with \code{tolerance_w}, \code{valid_w}, \code{reason_w}
#'   and the dispersion limit \code{disp_star}.
#' @noRd
.toleranceFromDispersion <- function(disp_w, n_w, k, c_noise, n_min,
                                     area_threshold){

  if( !.isScalarNumber(k) || k <= 0 ){
    stop(".toleranceFromDispersion(): k must be a single positive numeric.",
         call. = FALSE)
  }
  if( !.isScalarNumber(c_noise) || c_noise < 0 ){
    stop(".toleranceFromDispersion(): c_noise must be a single non-negative ",
         "numeric.", call. = FALSE)
  }
  if( !.isWholeNumber(n_min) || n_min < 1 ){
    stop(".toleranceFromDispersion(): n_min must be a positive integer.",
         call. = FALSE)
  }
  if( !is.numeric(area_threshold) || length(area_threshold) != 2L ||
      anyNA(area_threshold) || any(!is.finite(area_threshold)) ||
      area_threshold[1L] > area_threshold[2L] ){
    stop(".toleranceFromDispersion(): area_threshold must be c(lo, hi) with ",
         "lo <= hi.", call. = FALSE)
  }

  disp_star <- (area_threshold[2L] - area_threshold[1L]) / (2 * k)
  tolerance_w <- pmax(k * disp_w, c_noise)

  reason_w <- rep("", length(disp_w))
  reason_w[n_w < n_min] <- "n_below_min"
  reason_w[reason_w == "" & is.na(disp_w)] <- "disp_na"
  cross <- reason_w == "" & !is.na(disp_w) & disp_w >= disp_star
  reason_w[cross] <- "disp_crossover"
  valid_w <- reason_w == ""

  return(list(tolerance_w = tolerance_w,
              valid_w = valid_w,
              reason_w = reason_w,
              disp_star = disp_star))
}

#' Measure each window's signature without threshold filtering
#'
#' Calculate the median and scaled MAD across support pixels. Require two
#' valid pixels because one pixel cannot provide a dispersion estimate.
#'
#' @param vals_mat Numeric matrix with pixels in rows and windows in columns.
#' @returns A list with \code{mu_w}, \code{disp_w} and \code{n_w} for each window.
#' @noRd
.signatureFromMatrix <- function(vals_mat){

  if( !is.matrix(vals_mat) ){
    stop(".signatureFromMatrix(): vals_mat must be a matrix.", call. = FALSE)
  }
  nw <- ncol(vals_mat)
  mu_w <- rep(NA_real_, nw)
  disp_w <- rep(NA_real_, nw)
  n_w <- integer(nw)
  for( wi in seq_len(nw) ){
    v <- vals_mat[, wi]
    v <- v[is.finite(v)]
    n_w[wi] <- length(v)
    if( n_w[wi] >= 2L ){
      mu_w[wi] <- stats::median(v)
      disp_w[wi] <- 1.4826 * stats::median(abs(v - mu_w[wi]))
    }
  }
  return(list(mu_w = mu_w, disp_w = disp_w, n_w = n_w))
}

# These engineering defaults order strictness from loose to tight.
# Their spacing has not been externally validated.
.STRICTNESS_K <- c(loose = 3.0, balanced = 2.2, tight = 1.5)

#' Convert strictness to its robustness multiple
#' @noRd
.strictnessK <- function(strictness){

  if( !is.character(strictness) || length(strictness) != 1L ||
      !strictness %in% names(.STRICTNESS_K) ){
    stop(".strictnessK(): strictness must be one of ",
         paste0('"', names(.STRICTNESS_K), '"', collapse = ", "), ".",
         call. = FALSE)
  }
  return(.STRICTNESS_K[[strictness]])
}
