# =========================================================================
# FTW: field priors and agreement diagnostics
# =========================================================================
# A supplied polygon replaces the seed neighbourhood as calibration support.
# Its age also controls the optional penalty for growing outside the field.
# Agreement metrics compare masks without fitting calibration constants.
#
# ftw_prior holds polygon (sf, sfc or SpatVector) and optional source_year.
# Polygons are projected to the feature grid as needed.

#' Select the field containing the seed point
#'
#' A query can return neighbouring or overlapping fields. Keep the smallest
#' field intersecting the seed instead of combining unrelated fields.
#' @noRd
.selectFtwField <- function(polygons, rs){

  if( is.null(polygons) ){
    return(NULL)
  }
  p <- tryCatch({
    if( inherits(polygons, "sf") ){
      polygons
    }else if( inherits(polygons, "sfc") ){
      sf::st_sf(geometry = polygons)
    }else if( inherits(polygons, "SpatVector") ) { sf::st_as_sf(polygons) } else { NULL }
  }, error = function(e){
    return(NULL)
  })
  if( is.null(p) || nrow(p) == 0L ){
    return(NULL)
  }

  pt <- rs$spec$point
  if( !is.list(pt) || !.isScalarNumber(pt$longitude) ||
      !.isScalarNumber(pt$latitude) || is.na(sf::st_crs(p)) ){
    return(NULL)
  }
  seed <- sf::st_sfc(sf::st_point(c(pt$longitude, pt$latitude)), crs = 4326)
  seed <- tryCatch(sf::st_transform(seed, sf::st_crs(p)),
                   error = function(e){
                     return(NULL)
                   })
  if( is.null(seed) ){
    return(NULL)
  }
  hits <- tryCatch(sf::st_intersects(seed, p, sparse = TRUE)[[1L]],
                   error = function(e){
                     return(.ftwHitsByFeature(seed, p))
                   })
  if( length(hits) == 0L ){
    return(NULL)
  }
  if( length(hits) == 1L ){
    return(p[hits, , drop = FALSE])
  }

  candidates <- p[hits, , drop = FALSE]
  areas <- tryCatch(as.numeric(sf::st_area(candidates)),
                    error = function(e){
                      return(rep(Inf, nrow(candidates)))
                    })
  areas[!is.finite(areas)] <- Inf
  return(candidates[which.min(areas), , drop = FALSE])
}

#' Test polygons separately when the combined spatial test fails
#'
#' An invalid polygon can prevent the combined test. Ignore it only when its
#' bounding box excludes the seed. Otherwise the containing field is uncertain,
#' so return no match.
#' @noRd
.ftwHitsByFeature <- function(seed, p){

  geometry <- sf::st_geometry(p)
  seed_xy <- sf::st_coordinates(seed)[1L, ]
  hits <- integer(0L)
  for( i in seq_along(geometry) ){
    hit <- tryCatch(lengths(sf::st_intersects(seed, geometry[i])) > 0L,
                    error = function(e){
                      return(NA)
                    })
    if( is.na(hit) ){
      # Coordinates still provide a bounding box when the polygon test fails.
      box <- sf::st_bbox(geometry[i])
      in_box <- seed_xy[["X"]] >= box[["xmin"]] && seed_xy[["X"]] <= box[["xmax"]] &&
        seed_xy[["Y"]] >= box[["ymin"]] && seed_xy[["Y"]] <= box[["ymax"]]
      if( isTRUE(in_box) ){
        return(integer(0L))
      }
    }else if( hit ){
      hits <- c(hits, i)
    }
  }
  return(hits)
}

#' Project a usable field polygon to the feature grid's CRS
#'
#' Accept \code{sf}, \code{sfc} or \code{SpatVector} input. Return one
#' \code{SpatVector} geometry, or \code{NULL} if no polygon can be used.
#' @noRd
.ftwResolvePolygon <- function(rs, polygon = NULL){

  geom_template <- rs$internals$feature_stack
  if( is.null(geom_template) ){
    return(NULL)
  }
  poly <- polygon %||% rs$geom$ftw_prior$polygon
  if( is.null(poly) ){
    return(NULL)
  }

  v <- tryCatch({
    if( inherits(poly, "SpatVector") ){
      poly
    }else if( inherits(poly, c("sf", "sfc")) ) { terra::vect(poly) } else { NULL }
  }, error = function(e){
    return(NULL)
  })
  if( is.null(v) || terra::geomtype(v) != "polygons" || nrow(v) == 0L ){
    return(NULL)
  }

  # Select the seed's field before dissolving its parts.
  if( nrow(v) > 1L ){
    selected <- .selectFtwField(poly, rs)
    if( is.null(selected) ){
      return(NULL)
    }
    v <- tryCatch(terra::vect(selected), error = function(e){
      return(NULL)
    })
    if( is.null(v) ){
      return(NULL)
    }
  }
  v <- terra::aggregate(v)
  if( !identical(terra::crs(v), terra::crs(geom_template)) ){
    v <- terra::project(v, terra::crs(geom_template))
  }
  return(v)
}

#' Convert a field polygon to a mask on the feature grid
#'
#' Erode by \code{erode_px} pixels to exclude mixed boundary pixels, unless
#' that would empty the polygon. Return a logical matrix matching the grid,
#' or \code{NULL} if the polygon cannot be rasterised.
#' @noRd
.ftwPolygonMask <- function(rs, polygon = NULL, erode_px = 1L){

  v <- .ftwResolvePolygon(rs, polygon)
  if( is.null(v) ){
    return(NULL)
  }
  template <- rs$internals$feature_stack[[1L]]
  pixel_m <- as.numeric(terra::res(template)[1L])

  if( erode_px > 0L && is.finite(pixel_m) && pixel_m > 0 ){
    ve <- tryCatch(terra::buffer(v, width = -erode_px * pixel_m),
                   error = function(e){
                     return(v)
                   })
    if( !is.null(ve) && nrow(ve) > 0L && terra::expanse(ve) > 0 ){
      v <- ve
    }
  }

  rmask <- tryCatch(
    terra::rasterize(v, template, field = 1, background = NA),
    error = function(e){
      return(NULL)
    })
  if( is.null(rmask) ){
    return(NULL)
  }

  m <- matrix(!is.na(terra::values(rmask)[, 1L]),
              nrow = rs$geom$nr, ncol = rs$geom$nc, byrow = TRUE)
  if( !any(m) ){
    return(NULL)
  }
  return(m)
}

# =========================================================================
# FTW prior strength
# =========================================================================
# These engineering defaults have not been fitted to agreement metrics.
# Strength raises the local density requirement outside the field. Its
# half-life discounts older references. Sufficient spectral agreement can
# still allow growth across the boundary.
.FTW_PRIOR_STRENGTH <- 0.35
.FTW_TRUST_HALFLIFE_YR <- 5
.FTW_PRIOR_CROSS <- 0.75

#' Discount the field prior by its age
#'
#' Use the difference between the reference year and the midpoint of the
#' search period to halve the penalty every \code{.FTW_TRUST_HALFLIFE_YR} years.
#' Without a polygon or source year, return zero. An undated polygon can still
#' supply calibration support.
#' @noRd
.ftwPriorStrength <- function(rs){

  if( is.null(rs$geom$ftw_prior$polygon) ){
    return(0)
  }
  src_yr <- suppressWarnings(as.integer(rs$geom$ftw_prior$source_year))
  if( length(src_yr) != 1L || is.na(src_yr) ){
    return(0)
  }
  sd <- tryCatch(suppressWarnings(as.Date(rs$params$search_start_date %||%
                                          rs$spec$start_date)),
                 error = function(e){
                   return(as.Date(NA))
                 })
  ed <- tryCatch(suppressWarnings(as.Date(rs$params$search_end_date %||%
                                          rs$spec$end_date)),
                 error = function(e){
                   return(as.Date(NA))
                 })
  if( length(sd) != 1L || length(ed) != 1L ){
    return(0)
  }
  scene_yr <- NA_integer_
  if( !is.na(sd) && !is.na(ed) ){
    mid <- as.Date((as.numeric(sd) + as.numeric(ed)) / 2,
                   origin = "1970-01-01")
    scene_yr <- as.integer(format(mid, "%Y"))
  }else if( !is.na(sd) ){
    scene_yr <- as.integer(format(sd, "%Y"))
  }
  if( is.na(scene_yr) ){
    return(0)
  }
  dyr <- abs(scene_yr - src_yr)
  trust <- 0.5 ^ (dyr / .FTW_TRUST_HALFLIFE_YR)
  return(.FTW_PRIOR_STRENGTH * trust)
}

#' Measure calibration support inside a field polygon
#'
#' Exclude boundary pixels by erosion, then measure each window's signature.
#' Return the same structure as \code{.liteSupport()} for shared calibration.
#'
#' @returns A list with \code{sig}, \code{nb}, \code{precise_w} and \code{support},
#'   or \code{NULL} when the polygon or feature data provides no usable support.
#' @noRd
.ftwSupport <- function(rs, c_noise, erode_px = 1L){

  fstack <- rs$internals$feature_stack
  if( is.null(fstack) ){
    return(NULL)
  }

  mask <- .ftwPolygonMask(rs, erode_px = erode_px)
  if( is.null(mask) ){
    return(NULL)
  }

  # The polygon mask and feature values use the same grid and cell order.
  feat_vals <- terra::values(fstack)
  poly_idx <- which(as.vector(t(mask)))   # Match raster cell order.
  nb <- feat_vals[poly_idx, , drop = FALSE]
  if( nrow(nb) < 2L ){
    return(NULL)
  }

  sig <- .signatureFromMatrix(nb)
  se_w <- .C_MAD * sig$disp_w / sqrt(sig$n_w)
  precise_w <- if( is.na(c_noise) ){
    rep(NA, length(sig$n_w))
  }else{
    !is.na(se_w) & se_w <= c_noise
  }

  return(list(sig = sig,
              nb = nb,
              precise_w = precise_w,
              support = list(kind = "ftw", erode_px = erode_px,
                                   n_w = sig$n_w)))
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
#'   field boundary prior is used. A CRS is required for reprojection to the
#'   feature grid.
#' @returns A list with \code{iou}, \code{f1}, \code{precision}, \code{recall}
#'   (the ocular mask is treated as the prediction and the reference as the
#'   comparison mask), and the supporting pixel counts \code{n_pred},
#'   \code{n_ref}, \code{n_intersection}, and \code{n_union}.
#' @export
diagnose_against_ftw <- function(rs, reference = NULL){

  if( !is_rs(rs) ){
    stop("diagnose_against_ftw(): rs must be an ocular object.", call. = FALSE)
  }
  A <- rs$state$alive_mat
  if( is.null(A) || !is.matrix(A) ){
    stop("diagnose_against_ftw(): rs$state$alive_mat is not populated. Run ",
         "segment_area() / boundary_delineation() first.", call. = FALSE)
  }
  if( is.null(rs$internals$feature_stack) ){
    stop("diagnose_against_ftw(): rs$internals$feature_stack is required to ",
         "rasterise the reference polygon.", call. = FALSE)
  }

  reference <- reference %||% rs$geom$ftw_prior$polygon
  if( is.null(reference) ){
    stop("diagnose_against_ftw(): no reference polygon. Supply `reference` ",
         "or attach a field boundary with add_ftw_prior().", call. = FALSE)
  }
  missing_crs <- if( inherits(reference, "SpatVector") ){
    !nzchar(terra::crs(reference))
  }else if( inherits(reference, c("sf", "sfc")) ){
    is.na(sf::st_crs(reference))
  }else{
    FALSE
  }
  mask_error <- NULL
  B <- tryCatch(.ftwPolygonMask(rs, polygon = reference, erode_px = 0L),
                 error = function(e){
                   mask_error <<- e
                   return(NULL)
                 })
  if( is.null(B) ){
    if( missing_crs ){
      stop("diagnose_against_ftw(): the reference polygon has no CRS. ",
           "Set its coordinate reference system before comparing masks.", call. = FALSE)
    }
    if( !is.null(mask_error) ){
      stop(mask_error)
    }
    # Diagnose a disjoint reference after the existing mask rules reject it.
    v <- tryCatch({
      v <- if( inherits(reference, "SpatVector") ){
        reference
      }else if( inherits(reference, c("sf", "sfc")) ) { terra::vect(reference) } else { NULL }
      if( !is.null(v) && terra::geomtype(v) == "polygons" && nrow(v) > 0L ){
        terra::project(v, terra::crs(rs$internals$feature_stack))
      }else{
        NULL
      }
    }, error = function(e){
      return(NULL)
    })
    if( !is.null(v) ){
      bounds <- as.vector(terra::ext(v))
      grid <- as.vector(terra::ext(rs$internals$feature_stack))
      if( bounds[2L] < grid[1L] || bounds[1L] > grid[2L] ||
          bounds[4L] < grid[3L] || bounds[3L] > grid[4L] ){
        stop("diagnose_against_ftw(): the reference polygon does not overlap ",
             "the feature grid.", call. = FALSE)
      }
    }
    stop("diagnose_against_ftw(): the reference polygon could not be used ",
         "on the feature grid. Check its geometry and location.", call. = FALSE)
  }
  if( any(dim(A) != dim(B)) ){
    stop("diagnose_against_ftw(): mask dimensions disagree with the grid.",
         call. = FALSE)
  }

  A <- A & !is.na(A)
  inter <- sum(A & B)
  n_pred <- sum(A)
  n_ref <- sum(B)
  union <- sum(A | B)
  precision <- if( n_pred > 0L ) { inter / n_pred } else { NA_real_ }
  recall <- if( n_ref > 0L ) { inter / n_ref } else { NA_real_ }
  iou <- if( union > 0L ) { inter / union } else { NA_real_ }
  f1 <- if( !is.na(precision) && !is.na(recall) ){
    if( precision + recall > 0 ) { 2 * precision * recall / (precision + recall) } else { 0 }
  }else{
    NA_real_
  }

  return(list(iou = iou,
              f1 = f1,
              precision = precision,
              recall = recall,
              n_pred = n_pred,
              n_ref = n_ref,
              n_intersection = inter,
              n_union = union))
}
