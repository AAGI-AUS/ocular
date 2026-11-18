# =========================================================================
# Landsat-MODIS fusion
# =========================================================================

#' Order bands consistently for fusion
#'
#' Return the requested rasters in \code{assets} order, or \code{NULL} if a band
#' is missing. Fine and coarse bands must represent the same quantities.
#' @noRd
.orderedBands <- function(scene_bands, assets){

  b <- scene_bands[assets]
  if( any(vapply(b, is.null, logical(1L))) ){
    return(NULL)
  }
  names(b) <- assets
  return(b)
}

#' Read cropped remote bands into memory before repeated projection
#'
#' Multiplying by one reads each cropped raster once, avoiding repeated
#' network reads by \code{terra::project()}.
#' @noRd
.bandsToMemory <- function(bands){

  return(lapply(bands, function(r){
    return(r * 1)
  }))
}

#' Fetch coarse scenes near the fine scene dates
#'
#' Validation needs nearby MCD43A4 scenes, rather than the whole daily series.
#' Fetch within \code{pad_days} of each fine date and remove duplicate dates.
#' The default one-day margin tolerates a missing daily scene.
#'
#' @returns A list of coarse scenes, possibly empty.
#' @noRd
.fetchCoarseNearDates <- function(bbox, fine_dates, index_name, max_cloud_cover,
                                  source = "mcd43a4", pad_days = 1L){

  uniq <- sort(unique(as.numeric(fine_dates)))
  out <- list()
  errors <- character(0L)
  for( d in uniq ){
    dd <- as.Date(d, origin = "1970-01-01")
    sc <- tryCatch(
      .fetchStac(bbox = bbox,
                 start_date = as.character(dd - pad_days),
                 end_date = as.character(dd + pad_days),
                 index_name = index_name,
                 max_cloud_cover = max_cloud_cover,
                 source = source),
      error = function(e){
        errors <<- c(errors, paste0(as.character(dd), ": ",
                                    conditionMessage(e)))
        return(NULL)
      })
    if( !is.null(sc) && length(sc) ){
      out <- c(out, sc)
    }
  }
  if( length(errors) > 0L ){
    warning("MCD43A4 retrieval failed for ", length(errors),
            " fine-date window(s). First error: ", errors[[1L]],
            call. = FALSE)
  }
  if( !length(out) ){
    return(list())
  }
  keys <- vapply(out, function(s){
    return(as.character(s$date))
  }, character(1L))
  return(out[!duplicated(keys)])
}

#' Estimate fine reflectance from change in coarse reflectance
#'
#' Add the coarse change from t1 to t2 to the fine bands at t1, weighted by
#' similarity and bounded to physical reflectance. This is a STARFM-style update.
#'
#' @returns The estimated \code{bands} and the projected \code{coarse_t1} bands,
#'   both on the fine grid.
#' @noRd
.fuseDelta <- function(fine_bands, coarse_t1, coarse_t2){

  n <- length(fine_bands)
  out <- vector("list", n)
  c1_out <- vector("list", n)
  for( i in seq_len(n) ){
    f1 <- fine_bands[[i]]
    c1 <- terra::project(coarse_t1[[i]], f1, method = "bilinear")
    c2 <- terra::project(coarse_t2[[i]], f1, method = "bilinear")
    band_resid <- abs(f1 - c1)
    max_resid <- terra::global(band_resid, "max", na.rm = TRUE)[[1L]]
    if( !is.finite(max_resid) || max_resid < 1e-6 ){
      max_resid <- 1
    }
    weight <- terra::ifel((1 - band_resid / max_resid) < 0, 0,
                          1 - band_resid / max_resid)
    pred <- f1 + weight * (c2 - c1)
    pred <- terra::ifel(pred < 0, 0, pred)
    pred <- terra::ifel(pred > 1, 1, pred)
    out[[i]] <- pred
    c1_out[[i]] <- c1
  }
  return(list(bands = out, coarse_t1 = c1_out))
}

#' Calculate an index from fused reflectance bands
#'
#' Fine and coarse bands are paired by position. Apply the fine source's
#' index function to the fused bands, using its band names.
#'
#' @returns The fused index, its anchor index and the absolute fine/coarse
#'   index difference at the anchor. The latter two support fusion diagnostics.
#' @noRd
.fuseLandsatMODIS <- function(fine_bands, coarse_t1, coarse_t2, index_fun){

  fine_names <- names(fine_bands)
  delta <- .fuseDelta(fine_bands, coarse_t1, coarse_t2)
  fused <- delta$bands
  names(fused) <- fine_names

  # Use fine band names so the same index function applies to coarse bands.
  c1_fine <- delta$coarse_t1
  names(c1_fine) <- fine_names

  fine_stk <- do.call(c, unname(fine_bands))
  names(fine_stk) <- fine_names
  c1_stk <- do.call(c, unname(c1_fine))
  names(c1_stk) <- fine_names
  fused_stk <- do.call(c, unname(fused))
  names(fused_stk) <- fine_names

  # Report the anchor and its mismatch in index units.
  vi_fine_t1 <- .finiteIndex(index_fun(fine_stk))
  residual <- .finiteIndex(abs(vi_fine_t1 - index_fun(c1_stk)))
  return(list(index = .finiteIndex(index_fun(fused_stk)),
              bands = fused,
              residual = residual,
              landsat_t1 = vi_fine_t1))
}
