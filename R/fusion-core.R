# =========================================================================
# Landsat-MODIS fusion
# =========================================================================

#' Extract scene bands in a required asset order (internal)
#'
#' Returns a named list of the requested band rasters in \code{assets}
#' order, or NULL if any are missing. This gives fine and coarse reflectance
#' bands the same positional order for fusion.
#' @noRd
.orderedBands <- function(scene_bands, assets){
  b <- scene_bands[assets]
  if( any(vapply(b, is.null, logical(1L))) ) return(NULL)
  names(b) <- assets
  b
}

#' Force AOI-windowed band rasters into memory once (internal)
#'
#' Fetched bands are lazy /vsicurl/ handles cropped to the AOI, so every
#' downstream terra::project re-reads them over the network -- the dominant
#' cost in fusion. Multiplying by 1 reads each cropped raster into memory once;
#' subsequent projections then run locally.
#' @noRd
.bandsToMemory <- function(bands) lapply(bands, function(r) r * 1)

#' Fetch coarse scenes only near a set of fine-scene dates (internal)
#'
#' MCD43A4 is a daily product, so fetching a whole multi-month window retrieves
#' far more coarse scenes than leave-one-out validation consumes, which uses
#' the coarse scene nearest each fine-scene date. This fetches a tight
#' +/- \code{pad_days} window around each unique fine date and dedups by
#' date -- decoupling coarse-fetch cost from window length. Returns a flat
#' list of coarse scenes (possibly empty). \code{pad_days = 1} keeps one
#' coarse scene per fine date while tolerating a single-day product gap.
#' @noRd
.fetchCoarseNearDates <- function(bbox, fine_dates, index_name, max_cloud_cover,
                                  source = "mcd43a4", pad_days = 1L){
  uniq <- sort(unique(as.numeric(fine_dates)))
  out  <- list()
  errors <- character(0L)
  for( d in uniq ){
    dd <- as.Date(d, origin = "1970-01-01")
    sc <- tryCatch(
      .fetchStac(bbox            = bbox,
                 start_date      = as.character(dd - pad_days),
                 end_date        = as.character(dd + pad_days),
                 index_name      = index_name,
                 max_cloud_cover = max_cloud_cover,
                 source          = source),
      error = function(e){
        errors <<- c(errors, paste0(as.character(dd), ": ",
                                    conditionMessage(e)))
        NULL
      })
    if( !is.null(sc) && length(sc) ) out <- c(out, sc)
  }
  if( length(errors) > 0L )
    warning("MCD43A4 retrieval failed for ", length(errors),
            " fine-date window(s). First error: ", errors[[1L]],
            call. = FALSE)
  if( !length(out) ) return(list())
  keys <- vapply(out, function(s) as.character(s$date), character(1L))
  out[!duplicated(keys)]
}

#' Additive-delta band fusion -- STARFM-style fusion core (internal)
#'
#' Per reflectance band: fine at t1 plus a similarity-weighted coarse delta
#' (t1->t2), clamped to physical reflectance bounds (a STARFM-style additive
#' update). Returns a list of fused band rasters on the fine grid.
#' @noRd
.fuseDelta <- function(fine_bands, coarse_t1, coarse_t2){
  n   <- length(fine_bands)
  out <- vector("list", n)
  for( i in seq_len(n) ){
    f1 <- fine_bands[[i]]
    c1 <- terra::project(coarse_t1[[i]], f1, method = "bilinear")
    c2 <- terra::project(coarse_t2[[i]], f1, method = "bilinear")
    band_resid <- abs(f1 - c1)
    max_resid  <- terra::global(band_resid, "max", na.rm = TRUE)[[1L]]
    if( !is.finite(max_resid) || max_resid < 1e-6 ) max_resid <- 1
    weight <- terra::ifel((1 - band_resid / max_resid) < 0, 0,
                          1 - band_resid / max_resid)
    pred <- f1 + weight * (c2 - c1)
    pred <- terra::ifel(pred < 0, 0, pred)
    pred <- terra::ifel(pred > 1, 1, pred)
    out[[i]] <- pred
  }
  out
}

#' Fuse paired reflectance bands, then compute the index (internal)
#'
#' Fuses the fine and coarse reflectance bands with additive-delta
#' (\code{.fuseDelta}), then computes the index from the fused band stack
#' rather than fusing VI values directly. \code{fine_bands} are fine-source reflectance rasters on the
#' reference grid; \code{coarse_t1}/\code{coarse_t2} are coarse-source
#' reflectance rasters paired to \code{fine_bands} by position (callers
#' guarantee matching semantic band order). \code{index_fun} is the fine
#' source's index function, applied to the fused, fine-named stack.
#' \code{residual} and \code{landsat_t1} are returned in index units to
#' preserve the anchor-mismatch diagnostics consumed downstream by
#' \code{as_time_series()}.
#' @noRd
.fuseLandsatMODIS <- function(fine_bands, coarse_t1, coarse_t2, index_fun){

  fine_names <- names(fine_bands)
  fused <- .fuseDelta(fine_bands, coarse_t1, coarse_t2)
  names(fused) <- fine_names

  ## Coarse VI@t1 for residual reporting: coarse bands on the fine grid,
  ## carrying fine band names so the fine index_fun applies.
  c1_fine <- lapply(coarse_t1, function(r)
    terra::project(r, fine_bands[[1L]], method = "bilinear"))
  names(c1_fine) <- fine_names

  fine_stk  <- do.call(c, fine_bands); names(fine_stk)  <- fine_names
  c1_stk    <- do.call(c, c1_fine);    names(c1_stk)    <- fine_names
  fused_stk <- do.call(c, fused);      names(fused_stk) <- fine_names

  ## Reference and residual in index units: the fine VI at t1 and |fine VI at t1 -
  ## coarse VI@t1| (coarse stack carries fine band names so index_fun applies).
  vi_fine_t1 <- .finiteIndex(index_fun(fine_stk))
  residual   <- .finiteIndex(abs(vi_fine_t1 - index_fun(c1_stk)))
  return(list(index      = .finiteIndex(index_fun(fused_stk)),
              bands      = fused,
              residual   = residual,
              landsat_t1 = vi_fine_t1))
}
