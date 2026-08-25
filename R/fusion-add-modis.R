# =========================================================================
# Fusion: add MODIS
# =========================================================================

#' Add MCD43A4-informed estimates for dates not covered by Landsat
#'
#' For a Landsat \code{ocular} object, \code{add_modis()} estimates supported
#' spectral indices for dates represented by daily MCD43A4 NBAR data when no
#' Landsat acquisition falls within three days. Each estimate uses the nearest
#' Landsat acquisition as an anchor, is produced on the Landsat grid, and is
#' appended as an estimated scene. It is not a satellite observation.
#' Sentinel-2 objects are returned unchanged.
#'
#' This pathway is experimental. It uses MCD43A4 NBAR reflectance bands but
#' does not yet apply product-specific MCD43A4 quality layers.
#' \code{validate_data_fusion()} reports leave-one-out agreement with held-out
#' Landsat scenes for the selected data, place, and period. These metrics do not
#' establish general accuracy or replace independent validation. Supported
#' indices are NDVI, EVI, DVI, kNDVI, BNDVI, NIRv, SAVI, MSAVI, SR, EVI2,
#' OSAVI, ARVI, NBR, and BAI.
#'
#' @param x An \code{ocular} object, or a numeric longitude for standalone
#'   retrieval.
#' @param params Optional parameter list created by \code{rs_params()}.
#'   Forwarded to \code{get_rs()} in standalone mode. \code{add_modis()} does
#'   not use these parameters when an \code{ocular} object is supplied.
#' @param ... Additional arguments passed to \code{get_rs()} in standalone
#'   mode. No additional arguments are accepted when an \code{ocular} object
#'   is supplied.
#' @returns An \code{ocular} object, with estimated scenes appended when
#'   fusion is available and applicable.
#' @export
add_modis <- function(x, params = NULL, ...){
  if( is_rs(x) ){
    if( ...length() > 0L )
      stop("add_modis: piped mode does not accept additional arguments: ",
           paste(names(list(...)), collapse = ", "), ".", call. = FALSE)
    rs <- x
  }else{
    rs <- get_rs(longitude = x, params = params, ...)
  }
  if( is.null(rs$scenes) || length(rs$scenes) == 0L )
    stop("add_modis: rs has no scenes.", call. = FALSE)
  if( identical(rs$geom$source, "sentinel-2") ){
    message("add_modis: Sentinel-2 source -- no fusion needed.")
    return(rs)
  }
  is_fused <- vapply(rs$scenes, function(sc) isTRUE(sc$fused), logical(1L))
  fine_scenes <- rs$scenes[!is_fused]
  if( length(fine_scenes) == 0L )
    stop("add_modis: rs has no real fine-resolution anchor scenes.",
         call. = FALSE)
  ## Preserve a fine-only detection list before the idempotence guard returns.
  if( is.null(rs$internals$detection_scenes) )
    rs$internals$detection_scenes <- fine_scenes
  if( any(is_fused) ){
    message("add_modis: fused scenes are already present; rs unchanged.")
    return(rs)
  }
  ## A scenes-only pipe is intentionally lazy. Freeze the real fine-resolution
  ## scenes before appending synthetic fused scenes so subsequent boundary
  ## detection is identical to the eager path.
  index_name <- rs$spec$index_name
  ## Fusion requires an index computable from both the fine (Landsat) and
  ## coarse (MCD43A4) reflectance band sets. The index is validated for
  ## the fine source at get_rs(); here it must also exist in the coarse lookup.
  if( is.null(landsat_index_list[[index_name]]) ||
      is.null(mcd43a4_index_list[[index_name]]) ){
    warning("Coarse-tier fusion not supported for index '", index_name,
            "' (no Landsat/MCD43A4 reflectance mapping).", call. = FALSE)
    return(rs)
  }
  l_entry  <- landsat_index_list[[index_name]]
  m_assets <- mcd43a4_index_list[[index_name]]$assets

  log_msg <- function(...) message("    ", ...)
  ## Coarse tier: MCD43A4 daily 500 m NBAR (BRDF-normalised reflectance).
  coarse_scenes <- tryCatch(
    .fetchStac(bbox            = rs$geom$bbox,
               start_date      = rs$spec$start_date,
               end_date        = rs$spec$end_date,
               index_name      = index_name,
               max_cloud_cover = 100,
               source          = "mcd43a4"),
    error = function(e){
      warning("add_modis: MCD43A4 retrieval failed: ", conditionMessage(e),
              call. = FALSE)
      NULL
    })
  if( is.null(coarse_scenes) || length(coarse_scenes) == 0L ){
    log_msg("add_modis: no MCD43A4 scenes available.")
    return(rs)
  }
  log_msg(sprintf("add_modis: %d MCD43A4 scenes fetched", length(coarse_scenes)))

  ref_r <- fine_scenes[[1L]]$index
  ## Carry coarse reflectance bands in their native CRS for projection during
  ## fusion.
  coarse_b <- lapply(coarse_scenes, function(sc){
    b <- .orderedBands(sc$bands, m_assets)
    if( is.null(b) ) return(NULL)
    list(date = sc$date, bands = .bandsToMemory(b))
  })
  coarse_b <- Filter(Negate(is.null), coarse_b)
  if( length(coarse_b) == 0L ){
    log_msg("add_modis: MCD43A4 scenes lack the required NBAR bands.")
    return(rs)
  }

  ls_dates     <- vapply(fine_scenes, function(sc) as.numeric(sc$date), numeric(1))
  coarse_dates <- vapply(coarse_b, function(sc) as.numeric(sc$date), numeric(1))
  tol_d <- 3L

  fused <- list()
  fusion_errors <- character(0L)
  for( mi in seq_along(coarse_b) ){
    m_date <- coarse_dates[mi]
    if( any(abs(ls_dates - m_date) <= tol_d) ) next
    anchor_idx <- which.min(abs(ls_dates - m_date))
    ## Fine reflectance bands at the anchor, on the reference grid, in the
    ## fine index's band order.
    fine_raw <- .orderedBands(fine_scenes[[anchor_idx]]$bands, l_entry$assets)
    if( is.null(fine_raw) ) next
    fine_bands <- lapply(fine_raw, function(b)
      terra::project(b, ref_r, method = "bilinear"))
    names(fine_bands) <- l_entry$assets
    ci_t1 <- which.min(abs(coarse_dates - ls_dates[anchor_idx]))
    fr <- tryCatch(.fuseLandsatMODIS(fine_bands,
                                     coarse_b[[ci_t1]]$bands,
                                     coarse_b[[mi]]$bands,
                                     l_entry$fun),
                   error = function(e){
                     fusion_errors <<- c(
                       fusion_errors,
                       paste0(as.character(coarse_b[[mi]]$date), ": ",
                              conditionMessage(e)))
                     NULL
                   })
    if( is.null(fr) ) next
    fused[[length(fused) + 1L]] <- list(
      date        = as.Date(m_date, origin = "1970-01-01"),
      index       = fr$index,
      bands       = fr$bands,    ## retained for downstream output and validation
      item_assets = list(),      ## synthetic scene -- no STAC assets to re-fetch
      fused       = TRUE,
      anchor_days = as.integer(abs(m_date - ls_dates[anchor_idx])),
      residual    = fr$residual,
      landsat_t1  = fr$landsat_t1)
  }
  if( length(fusion_errors) > 0L )
    warning("add_modis: fusion failed for ", length(fusion_errors),
            " candidate scene(s). First error: ", fusion_errors[[1L]],
            call. = FALSE)
  log_msg(sprintf("add_modis: %d fused scenes created", length(fused)))

  base_tagged <- lapply(fine_scenes, function(sc){
    sc$fused       <- isTRUE(sc$fused)
    sc$anchor_days <- sc$anchor_days %||% NA_integer_
    sc$residual    <- sc$residual    %||% NULL
    sc$landsat_t1  <- sc$landsat_t1  %||% NULL
    sc
  })
  all_scenes <- c(base_tagged, fused)
  ord <- order(vapply(all_scenes, function(sc) as.numeric(sc$date), numeric(1)))
  rs$scenes <- all_scenes[ord]
  return(rs)
}
