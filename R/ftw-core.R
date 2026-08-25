# =========================================================================
# FTW: guided field boundary delineation
# =========================================================================
# Fields of The World (FTW) integration has two roles:
#
#   1. FTW as calibration support. When the user supplies
#      a reference field polygon in rs$geom$ftw_prior$polygon, the in-field
#      signature is measured on the polygon interior instead of the seed
#      neighbourhood. .ftwSupport() returns the same shape as .liteSupport(),
#      so it feeds the identical tolerance / validity / baseline machinery in
#      .calibrateScene() -- only the measurement support changes.
#
#   2. diagnose_against_ftw(): pixel-wise agreement (IoU / F1 / precision /
#      recall) between the ocular delineation and a reference polygon. This is
#      an agreement diagnostic, not a tuning signal; calibration does not fit
#      constants to these metrics.
#
# The optional soft-prior penalty is wired into .floodFill, and ftw-query.R
# provides an explicit-source DuckDB query/cache path. A 500 m confidence
# product is not yet incorporated; the current penalty is instead discounted
# by the reference boundary's age. Remote queries remain opt-in and are not
# exercised by the offline package test suite.
#
# Expected ftw_prior shape:
#   rs$geom$ftw_prior <- list(polygon = <sf, sfc, or terra SpatVector>,
#                             source_year = <optional integer>)
# The polygon may be in any CRS; it is reprojected onto the feature grid.

#' Select the FTW polygon containing the rs seed point (internal)
#'
#' A bbox query commonly returns neighbouring fields. Delineation is seeded by
#' one point, so attaching/dissolving every intersecting polygon would turn the
#' reference into a multi-field prior. This helper returns only polygons that
#' intersect the seed and, when boundaries overlap, deterministically chooses
#' the smallest candidate (the most specific enclosing field).
#'
#' @noRd
.selectFtwField <- function(polygons, rs){
  if( is.null(polygons) ) return(NULL)
  p <- tryCatch({
    if( inherits(polygons, "sf") ) polygons
    else if( inherits(polygons, "sfc") ) sf::st_sf(geometry = polygons)
    else if( inherits(polygons, "SpatVector") ) sf::st_as_sf(polygons)
    else NULL
  }, error = function(e) NULL)
  if( is.null(p) || nrow(p) == 0L ) return(NULL)

  pt <- rs$spec$point
  if( !is.list(pt) || !.isScalarNumber(pt$longitude) ||
      !.isScalarNumber(pt$latitude) || is.na(sf::st_crs(p)) )
    return(NULL)
  seed <- sf::st_sfc(sf::st_point(c(pt$longitude, pt$latitude)), crs = 4326)
  seed <- tryCatch(sf::st_transform(seed, sf::st_crs(p)),
                   error = function(e) NULL)
  if( is.null(seed) ) return(NULL)
  hits <- tryCatch(sf::st_intersects(seed, p, sparse = TRUE)[[1L]],
                   error = function(e) integer(0L))
  if( length(hits) == 0L ) return(NULL)
  if( length(hits) == 1L ) return(p[hits, , drop = FALSE])

  candidates <- p[hits, , drop = FALSE]
  areas <- tryCatch(as.numeric(sf::st_area(candidates)),
                    error = function(e) rep(Inf, nrow(candidates)))
  areas[!is.finite(areas)] <- Inf
  candidates[which.min(areas), , drop = FALSE]
}

#' Resolve and validate the FTW reference polygon (internal)
#'
#' Returns a single-geometry terra SpatVector in the feature-stack CRS, or
#' \code{NULL} when no usable polygon is present. Accepts sf, sfc, or
#' SpatVector input.
#'
#' @noRd
.ftwResolvePolygon <- function(rs, polygon = NULL){
  geom_template <- rs$internals$feature_stack
  if( is.null(geom_template) ) return(NULL)
  poly <- polygon %||% rs$geom$ftw_prior$polygon
  if( is.null(poly) ) return(NULL)

  v <- tryCatch({
    if( inherits(poly, "SpatVector") )            poly
    else if( inherits(poly, c("sf", "sfc")) )     terra::vect(poly)
    else                                          NULL
  }, error = function(e) NULL)
  if( is.null(v) || terra::geomtype(v) != "polygons" || nrow(v) == 0L )
    return(NULL)

  ## A query may return neighbours. Keep only the field containing the seed;
  ## never dissolve unrelated fields into one reference unit.
  if( nrow(v) > 1L ){
    selected <- .selectFtwField(poly, rs)
    if( is.null(selected) ) return(NULL)
    v <- tryCatch(terra::vect(selected), error = function(e) NULL)
    if( is.null(v) ) return(NULL)
  }
  ## A selected field can itself be multipart; dissolve that one unit.
  v <- terra::aggregate(v)
  if( !identical(terra::crs(v), terra::crs(geom_template)) )
    v <- terra::project(v, terra::crs(geom_template))
  v
}

#' Rasterise an FTW polygon onto the feature grid (internal)
#'
#' Reprojects, optionally erodes inward by \code{erode_px} pixels (to drop
#' spectrally mixed boundary pixels), and rasterises onto the feature-stack
#' grid. Returns a logical \code{nr x nc} matrix (in-polygon = TRUE) aligned
#' to \code{rs$geom}, or \code{NULL} on failure. Erosion falls back to the
#' un-eroded polygon when it would empty the geometry.
#'
#' @noRd
.ftwPolygonMask <- function(rs, polygon = NULL, erode_px = 1L){
  v        <- .ftwResolvePolygon(rs, polygon)
  if( is.null(v) ) return(NULL)
  template <- rs$internals$feature_stack[[1L]]
  pixel_m  <- as.numeric(terra::res(template)[1L])

  if( erode_px > 0L && is.finite(pixel_m) && pixel_m > 0 ){
    ve <- tryCatch(terra::buffer(v, width = -erode_px * pixel_m),
                   error = function(e) v)
    if( !is.null(ve) && nrow(ve) > 0L && terra::expanse(ve) > 0 ) v <- ve
  }

  rmask <- tryCatch(
    terra::rasterize(v, template, field = 1, background = NA),
    error = function(e) NULL)
  if( is.null(rmask) ) return(NULL)

  m <- matrix(!is.na(terra::values(rmask)[, 1L]),
              nrow = rs$geom$nr, ncol = rs$geom$nc, byrow = TRUE)
  if( !any(m) ) return(NULL)
  m
}

# -------------------------------------------------------------------------
# FTW soft-prior penalty: constants + effective strength
# -------------------------------------------------------------------------
# Fixed a-priori engineering defaults (never IoU-fitted; external scientific
# validation is still required):
#   .FTW_PRIOR_STRENGTH    -- how much an outside-field prior raises the
#                             local-density admission bar at full trust.
#   .FTW_TRUST_HALFLIFE_YR -- exponential half-life (years) over which the
#                             reference boundary's age discounts that strength.
#   .FTW_PRIOR_CROSS       -- per-window evidence fraction (n_match / W) at or
#                             above which strong reflectance evidence crosses
#                             the prior (matches the .floodFill default).
.FTW_PRIOR_STRENGTH    <- 0.35
.FTW_TRUST_HALFLIFE_YR <- 5
.FTW_PRIOR_CROSS       <- 0.75

#' Effective FTW soft-prior strength for a detection period (internal)
#'
#' Returns the age-decayed penalty strength: a conservative base penalty
#' discounted by the age of the reference boundary relative to the midpoint of
#' the detection period via an exponential half-life,
#' trust = 0.5^(|detection_year - source_year| / halflife).
#' Returns 0 (a flood-fill no-op) when no usable FTW polygon or source year is
#' present. An undated boundary remains usable as calibration support, but it
#' is not safe to impose an age-weighted prior when its age is unknown. The
#' detection year is taken from the midpoint of the configured date window
#' (\code{search_*_date} overrides, otherwise \code{rs$spec$*_date}).
#'
#' @noRd
.ftwPriorStrength <- function(rs){
  if( is.null(rs$geom$ftw_prior$polygon) ) return(0)
  src_yr <- suppressWarnings(as.integer(rs$geom$ftw_prior$source_year))
  if( length(src_yr) != 1L || is.na(src_yr) ) return(0)
  sd <- tryCatch(suppressWarnings(as.Date(rs$params$search_start_date %||%
                                          rs$spec$start_date)),
                 error = function(e) as.Date(NA))
  ed <- tryCatch(suppressWarnings(as.Date(rs$params$search_end_date %||%
                                          rs$spec$end_date)),
                 error = function(e) as.Date(NA))
  if( length(sd) != 1L || length(ed) != 1L ) return(0)
  scene_yr <- NA_integer_
  if( !is.na(sd) && !is.na(ed) ){
    mid      <- as.Date((as.numeric(sd) + as.numeric(ed)) / 2,
                        origin = "1970-01-01")
    scene_yr <- as.integer(format(mid, "%Y"))
  }else if( !is.na(sd) ){
    scene_yr <- as.integer(format(sd, "%Y"))
  }
  if( is.na(scene_yr) ) return(0)
  dyr   <- abs(scene_yr - src_yr)
  trust <- 0.5 ^ (dyr / .FTW_TRUST_HALFLIFE_YR)
  .FTW_PRIOR_STRENGTH * trust
}


#' FTW-derived measurement support (internal)
#'
#' Measures the threshold-free per-window signature on the interior of the
#' FTW reference polygon (eroded to exclude boundary-mixed pixels), returning
#' the same structure as \code{.liteSupport()} so \code{.calibrateScene()}
#' consumes it through the identical downstream path.
#'
#' @returns list(sig, nb, precise_w, support), or
#'   \code{NULL} when no usable polygon / feature data is available.
#' @noRd
.ftwSupport <- function(rs, c_noise, erode_px = 1L){
  fstack <- rs$internals$feature_stack
  if( is.null(fstack) ) return(NULL)

  mask <- .ftwPolygonMask(rs, erode_px = erode_px)
  if( is.null(mask) ) return(NULL)

  ## Feature values at the interior pixels (cell order matches the mask,
  ## since the mask was rasterised onto this same grid).
  feat_vals <- terra::values(fstack)
  poly_idx  <- which(as.vector(t(mask)))   ## t(): matrix is row-major (byrow)
  nb        <- feat_vals[poly_idx, , drop = FALSE]
  if( nrow(nb) < 2L ) return(NULL)

  sig  <- .signatureFromMatrix(nb)
  se_w <- .C_MAD * sig$disp_w / sqrt(sig$n_w)
  precise_w <- if( is.na(c_noise) ) rep(NA, length(sig$n_w))
               else !is.na(se_w) & se_w <= c_noise

  list(sig           = sig,
       nb            = nb,
       precise_w     = precise_w,
       support       = list(kind = "ftw", erode_px = erode_px,
                            n_w = sig$n_w))
}

#' Compare an ocular delineation with a field boundary reference
#'
#' Computes pixel-wise agreement between a delineated field and a field
#' boundary reference on the feature grid. By default, it uses the optional
#' prior attached by \code{add_ftw_prior()}. The function does not alter the
#' \code{ocular} object or tune any parameters. When the attached polygon
#' supplied calibration support or a soft prior, agreement with that polygon is
#' an in-sample diagnostic, not independent accuracy. Supply a separate
#' held-out \code{reference} when making accuracy claims.
#'
#' @param rs A delineated \code{ocular} object, for example after
#'   \code{segment_area()} or \code{boundary_delineation()}.
#' @param reference Optional reference polygon supplied as an \code{sf},
#'   \code{sfc}, or \code{SpatVector} object. When \code{NULL}, the attached
#'   field boundary prior is used.
#' @returns A list with \code{iou}, \code{f1}, \code{precision}, \code{recall}
#'   (the ocular mask is treated as the prediction and the reference as the
#'   comparison mask), and the supporting pixel counts \code{n_pred},
#'   \code{n_ref}, \code{n_intersection}, and \code{n_union}.
#' @export
diagnose_against_ftw <- function(rs, reference = NULL){
  if( !is_rs(rs) )
    stop("diagnose_against_ftw(): rs must be an ocular object.", call. = FALSE)
  A <- rs$state$alive_mat
  if( is.null(A) || !is.matrix(A) )
    stop("diagnose_against_ftw(): rs$state$alive_mat is not populated; run ",
         "segment_area() / boundary_delineation() first.", call. = FALSE)
  if( is.null(rs$internals$feature_stack) )
    stop("diagnose_against_ftw(): rs$internals$feature_stack is required to ",
         "rasterise the reference polygon.", call. = FALSE)

  B <- .ftwPolygonMask(rs, polygon = reference, erode_px = 0L)
  if( is.null(B) )
    stop("diagnose_against_ftw(): no usable reference polygon (supply ",
         "`reference` or set rs$geom$ftw_prior$polygon).", call. = FALSE)
  if( any(dim(A) != dim(B)) )
    stop("diagnose_against_ftw(): mask dimensions disagree with the grid.",
         call. = FALSE)

  A <- A & !is.na(A)
  inter <- sum(A & B)
  n_pred <- sum(A); n_ref <- sum(B); union <- sum(A | B)
  precision <- if( n_pred > 0L ) inter / n_pred else NA_real_
  recall    <- if( n_ref  > 0L ) inter / n_ref  else NA_real_
  iou       <- if( union  > 0L ) inter / union  else NA_real_
  f1        <- if( !is.na(precision) && !is.na(recall) ){
    if( precision + recall > 0 )
      2 * precision * recall / (precision + recall) else 0
  }else NA_real_

  list(iou            = iou,
       f1             = f1,
       precision      = precision,
       recall         = recall,
       n_pred         = n_pred,
       n_ref          = n_ref,
       n_intersection = inter,
       n_union        = union)
}
