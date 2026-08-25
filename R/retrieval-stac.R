# =========================================================================
# Retrieval: source resolution and STAC fetch
# =========================================================================

#' Source-config lookup (internal)
#' @noRd
.stacCfg <- function(source){

  ## `grid$tiling` and `mask_kind` lift two behaviours out of inline
  ## collection-string checks downstream: tiling ("utm" | "sinusoidal")
  ## drives bbox padding and per-tile mosaicking in .fetchStac; mask_kind
  ## ("scl" | "bit" | NULL) drives cloud-mask interpretation. MCD43A4 and
  ## any future sinusoidal coarse source are represented by configuration
  ## entries rather than source-specific fetch branches.
  return(switch(source,
                "sentinel-2" = list(
                  endpoint = .STAC_ENDPOINT, collection = "sentinel-2-l2a",
                  cloud_key = "eo:cloud_cover", platforms = NULL,
                  scale_factor = 0.0001, scale_offset = 0,
                  valid_range = c(1, 10000), mask_asset = "SCL",
                  grid = list(tiling = "utm"), mask_kind = "scl"),
                "landsat-8" = list(
                  endpoint = .STAC_ENDPOINT, collection = "landsat-c2-l2",
                  cloud_key = "eo:cloud_cover", platforms = c("landsat-8", "landsat-9"),
                  scale_factor = 0.0000275, scale_offset = -0.2,
                  valid_range = c(7273, 43636), mask_asset = "qa_pixel",
                  grid = list(tiling = "utm"), mask_kind = "bit"),
                "landsat-5" = list(
                  endpoint = .STAC_ENDPOINT, collection = "landsat-c2-l2",
                  cloud_key = "eo:cloud_cover", platforms = c("landsat-5", "landsat-7"),
                  scale_factor = 0.0000275, scale_offset = -0.2,
                  valid_range = c(7273, 43636), mask_asset = "qa_pixel",
                  grid = list(tiling = "utm"), mask_kind = "bit"),
                ## MCD43A4 V6.1 daily 500 m NBAR (BRDF-normalised reflectance, bands 1-7;
                ## Band1 red, Band2 NIR, Band3 blue). Preview collection on Planetary
                ## Computer. Fusion combines reflectance bands before computing
                ## an index.
                ## valid_range/fill: confirm against MCD43A4 spec (NBAR fill = 32767).
                "mcd43a4" = list(
                  endpoint = .STAC_ENDPOINT, collection = "modis-43A4-061",
                  cloud_key = NULL, platforms = NULL,
                  scale_factor = 0.0001, scale_offset = 0,
                  valid_range = c(0, 32766), mask_asset = NULL,
                  grid = list(tiling = "sinusoidal"), mask_kind = NULL),
                stop("Unknown source: ", source, call. = FALSE)))
}

#' Resolve source / index / pixel size from start_date and index_name (internal)
#' @noRd
.resolveSource <- function(start_date, index_name, source = "auto"){

  d <- as.Date(start_date)
  if( is.na(d) )
    stop("start_date could not be parsed: '", start_date, "'", call. = FALSE)
  source <- match.arg(source,
                      c("auto", "sentinel-2", "landsat-8", "landsat-5"))
  if( identical(source, "auto") ){
    source <- if( d >= .S2_START ) "sentinel-2"
    else if( d >= .L8_START ) "landsat-8" else "landsat-5"
    if( identical(source, "landsat-8") )
      message("Sentinel-2 not available before ",
              format(.S2_START, "%Y-%m-%d"), "; Landsat 8/9 used.")
    if( identical(source, "landsat-5") )
      message("Landsat 8 unavailable before ",
              format(.L8_START, "%Y-%m-%d"), "; Landsat 5/7 used.")
  }
  if( identical(source, "sentinel-2") ){
    if( d < .S2_START )
      stop("Sentinel-2 is unavailable before ", format(.S2_START), ".",
           call. = FALSE)
    idx <- s2_index_list[[index_name]]
    if( is.null(idx) )
      stop("Unknown index_name '", index_name, "'. Available: ",
           paste(names(s2_index_list), collapse = ", "), call. = FALSE)
    return(list(source = "sentinel-2", cfg = .stacCfg("sentinel-2"), idx = idx))
  }
  if( identical(source, "landsat-8") ){
    if( d < .L8_START )
      stop("Landsat 8/9 is unavailable before ", format(.L8_START), ".",
           call. = FALSE)
    idx <- landsat_index_list[[index_name]]
    if( is.null(idx) )
      stop("Unknown index_name '", index_name, "'. Available: ",
           paste(names(landsat_index_list), collapse = ", "), call. = FALSE)
    return(list(source = "landsat-8", cfg = .stacCfg("landsat-8"), idx = idx))
  }
  idx <- landsat_index_list[[index_name]]
  if( is.null(idx) )
    stop("Unknown index_name '", index_name, "'.", call. = FALSE)
  return(list(source = "landsat-5", cfg = .stacCfg("landsat-5"), idx = idx))
}

#' Merge adjacent raster tiles, reprojecting only when their CRSs differ
#' @noRd
.mergeRasterTiles <- function(rasters, method = "bilinear"){
  rasters <- Filter(Negate(is.null), rasters)
  if( length(rasters) == 0L ) return(NULL)
  if( length(rasters) == 1L ) return(rasters[[1L]])
  ref_crs <- terra::crs(rasters[[1L]])
  aligned <- lapply(rasters, function(r){
    if( identical(terra::crs(r), ref_crs) ) r
    else terra::project(r, ref_crs, method = method)
  })
  merged <- do.call(terra::merge,
                    c(aligned, list(algo = 3, method = method,
                                    na.rm = TRUE)))
  ## algo=3 commonly returns a temporary VRT. Materialise the cropped AOI so
  ## an rs object does not retain a fragile reference to that temporary file.
  merged * 1
}

#' Mosaic same-date fine-source items into one acquisition scene
#' @noRd
.mosaicFineScenesByDate <- function(scenes, lookup_entry){
  if( length(scenes) < 2L ) return(scenes)
  keys <- vapply(scenes, function(sc) as.character(sc$date), character(1L))
  groups <- split(scenes, keys)
  out <- lapply(groups, function(group){
    bands <- lapply(lookup_entry$assets, function(an)
      .mergeRasterTiles(lapply(group, function(sc) sc$bands[[an]])))
    names(bands) <- lookup_entry$assets
    if( any(vapply(bands, is.null, logical(1L))) ) return(NULL)

    masks <- lapply(group, `[[`, "pixel_mask")
    masks <- Filter(Negate(is.null), masks)
    pixel_mask <- if( length(masks) == 0L ) NULL
                  else .mergeRasterTiles(masks, method = "near")

    asset_names <- unique(unlist(lapply(group, function(sc)
      names(sc$item_assets)), use.names = FALSE))
    item_assets <- lapply(asset_names, function(an){
      unique(unlist(lapply(group, function(sc) sc$item_assets[[an]]),
                    use.names = FALSE))
    })
    names(item_assets) <- asset_names

    list(date = group[[1L]]$date,
         bands = bands,
         index = .computeIndexFromBands(bands, lookup_entry),
         item_assets = item_assets,
         pixel_mask = pixel_mask)
  })
  Filter(Negate(is.null), out)
}

#' STAC items fetch via CQL2 + retry-with-backoff (internal)
#'
#' Builds a CQL2 filter via `rstac::ext_filter`, executes via
#' `post_request() |> items_fetch()`, and signs assets. The whole chain is
#' wrapped in retry-with-backoff per MPC's official guidance for handling
#' 503/504 under load. The cloud-cover filter is added server-side for sources
#' with a `cloud_key`; MCD43A4 has no scene-level `eo:cloud_cover` field.
#'
#' Closure-return pattern (`list(items, err)`) avoids `<<-` for capturing
#' the inner-tryCatch error to surface in the final message.
#'
#' @noRd
.fetchStacItems <- function(stac_url, collection_id, bbox_vec,
                            start_date, end_date,
                            max_cloud_cover, has_cloud,
                            retry_max     = 10L,
                            initial_pause = 2,
                            max_pause     = 60){

  bbox_geo <- rstac::cql2_bbox_as_geojson(bbox_vec)
  dt_int   <- rstac::cql2_interval(start_date, end_date)
  signer   <- rstac::sign_planetary_computer()

  attempt <- 0L
  repeat {
    res <- tryCatch({
      q <- if( has_cloud ){
        rstac::stac(stac_url, force_version = "1.0.0") |>
          rstac::ext_filter(
            collection == {{ collection_id }} &&
              t_intersects(datetime, {{ dt_int }}) &&
              s_intersects(geometry, {{ bbox_geo }}) &&
              `eo:cloud_cover` <= {{ max_cloud_cover }})
      } else {
        rstac::stac(stac_url, force_version = "1.0.0") |>
          rstac::ext_filter(
            collection == {{ collection_id }} &&
              t_intersects(datetime, {{ dt_int }}) &&
              s_intersects(geometry, {{ bbox_geo }}))
      }
      items <- rstac::post_request(q) |>
        rstac::items_fetch(progress = FALSE) |>
        rstac::items_sign(signer)
      list(items = items, err = NULL)
    }, error = function(e) list(items = NULL, err = e))

    if( !is.null(res$items) ) return(res$items)
    attempt <- attempt + 1L
    if( attempt >= retry_max ) break
    Sys.sleep(min(max_pause, initial_pause * 2^(attempt - 1L)))
  }
  message("    STAC query failed after ", retry_max, " retries: ",
          conditionMessage(res$err))
  return(NULL)
}

#' Consolidated STAC fetch (internal)
#'
#' Absorbs source config, query construction, query execution, scene reading.
#' Returns scenes harmonised to the first scene's grid.
#'
#' @noRd
.fetchStac <- function(bbox,
                       start_date,
                       end_date,
                       index_name      = "EVI2",
                       max_cloud_cover = 50,
                       scl_classes     = NULL,
                       source          = NULL){

  if( !.isScalarNumber(max_cloud_cover) ||
      max_cloud_cover < 0 || max_cloud_cover > 100 )
    stop(".fetchStac(): max_cloud_cover must be in [0, 100].",
         call. = FALSE)
  if( !inherits(bbox, "bbox") || anyNA(bbox) ||
      any(!is.finite(unclass(bbox))) )
    stop(".fetchStac(): bbox must be a finite sf bounding box.",
         call. = FALSE)
  if( !requireNamespace("rstac", quietly = TRUE) )
    stop("rstac package required.", call. = FALSE)

  res <- if( is.null(source) ) .resolveSource(start_date, index_name)
  else list(source = source, cfg = .stacCfg(source),
            idx = if( source == "mcd43a4" )        mcd43a4_index_list[[index_name]]
            else if( grepl("^landsat", source) )   landsat_index_list[[index_name]]
            else                                   s2_index_list[[index_name]])
  cfg <- res$cfg
  idx <- res$idx
  if( is.null(idx) )
    stop("Unknown index_name '", index_name, "' for source '", res$source, "'.",
         call. = FALSE)

  ## Grid behaviour is determined by cfg$grid$tiling from .stacCfg.
  is_sinusoidal <- identical(cfg$grid$tiling, "sinusoidal")

  ## MCD43A4 uses MODIS Sinusoidal tiles, which are 10 degrees square at the
  ## equator. STAC bbox search returns items whose WGS-84 footprint intersects
  ## the search bbox, but small AOIs near a tile seam may match only one of the
  ## two adjacent tiles based on footprint geometry alone -- even when the AOI's
  ## Sinusoidal coordinates lie in the neighbouring tile. Padding the search
  ## bbox by 5 km includes neighbouring coverage near a seam; the
  ## mosaic-per-date step downstream merges returned tiles into a single
  ## Sinusoidal raster covering the seam. Padding applies only to the search;
  ## the padded box is not used to enlarge the AOI crop.
  bbox_search <- bbox
  if( is_sinusoidal ){
    pad_deg <- 5000 / 111320
    bbox_search <- sf::st_bbox(c(xmin = max(-180, bbox[["xmin"]] - pad_deg),
                                 ymin = max(-90,  bbox[["ymin"]] - pad_deg),
                                 xmax = min(180, bbox[["xmax"]] + pad_deg),
                                 ymax = min(90,  bbox[["ymax"]] + pad_deg)),
                               crs = sf::st_crs(4326))
  }
  bbox_vec <- c(bbox_search[["xmin"]], bbox_search[["ymin"]],
                bbox_search[["xmax"]], bbox_search[["ymax"]])

  items <- .fetchStacItems(
    stac_url        = cfg$endpoint,
    collection_id   = cfg$collection,
    bbox_vec        = bbox_vec,
    start_date      = start_date,
    end_date        = end_date,
    max_cloud_cover = max_cloud_cover,
    has_cloud       = !is.null(cfg$cloud_key))
  if( is.null(items) || length(items[["features"]]) == 0L ) return(NULL)

  ## Cloud filter is handled server-side via CQL2 in .fetchStacItems.
  ## Platform filter is retained client-side: pushing it to CQL2 would
  ## use additional ext_filter syntax for marginal payload reduction.
  features <- items[["features"]]
  if( !is.null(cfg$platforms) ){
    features <- Filter(function(f){
      plat <- f[["properties"]][["platform"]]
      !is.null(plat) && plat %in% cfg$platforms
    }, features)
  }
  if( length(features) == 0L ) return(NULL)

  aoi_sfc <- sf::st_as_sfc(bbox)

  read_band <- function(asset_href, scale){
    if( is.null(asset_href) ) return(NULL)
    b <- terra::rast(paste0("/vsicurl/", asset_href))
    terra::scoff(b) <- NULL
    ## Crop each read to the AOI window. An MCD43A4 Sinusoidal tile is 2400 x 2400
    ## cells, while an AOI uses a small subset. A STAC bbox search can also
    ## return a tile whose data extent does not cover the AOI, so test overlap
    ## first and drop non-overlapping tiles. Multi-tile AOIs are reconstructed
    ## by the per-asset merge downstream.
    aoi_in_crs <- sf::st_transform(aoi_sfc, terra::crs(b))
    aoi_ext    <- terra::ext(sf::st_bbox(aoi_in_crs))
    te         <- terra::ext(b)
    overlaps   <- terra::xmin(aoi_ext) < terra::xmax(te) &&
      terra::xmax(aoi_ext) > terra::xmin(te) &&
      terra::ymin(aoi_ext) < terra::ymax(te) &&
      terra::ymax(aoi_ext) > terra::ymin(te)
    if( !overlaps ) return(NULL)
    b <- terra::crop(b, aoi_ext, snap = "out")
    if( terra::ncell(b) == 0 ) return(NULL)
    if( scale ){
      if( !is.null(cfg$valid_range) )
        b <- terra::ifel(b < cfg$valid_range[1] | b > cfg$valid_range[2], NA, b)
      if( !is.null(cfg$scale_factor) )
        b <- b * cfg$scale_factor + (cfg$scale_offset %||% 0)
    }
    b
  }

  ## ---- Sinusoidal-tiled path: group by date ->
  ## mosaic each asset in Sinusoidal coordinates, then compute the index.
  ## Stored scenes carry raw bands alongside the computed index so later callers
  ## (different search_index) can reuse the bands without re-fetching.
  if( is_sinusoidal ){
    feat_dates <- vapply(features, function(f){
      ds <- f[["properties"]][["start_datetime"]] %||%
        f[["properties"]][["datetime"]] %||% ""
      substr(ds, 1, 10)
    }, character(1L))

    scenes <- lapply(unique(feat_dates), function(date_str){
      group <- features[feat_dates == date_str]
      dt    <- tryCatch(as.Date(date_str), error = function(e) NULL)
      if( is.null(dt) || is.na(dt) ) return(NULL)

      ## Read bands per tile, organised by asset name.
      per_tile_bands <- lapply(group, function(f){
        tryCatch({
          tb <- lapply(idx$assets, function(an){
            read_band(f[["assets"]][[an]][["href"]], scale = TRUE)
          })
          names(tb) <- idx$assets
          Filter(Negate(is.null), tb)
        }, error = function(e) NULL)
      })
      per_tile_bands <- Filter(function(x) !is.null(x) && length(x) > 0L,
                               per_tile_bands)
      if( length(per_tile_bands) == 0L ) return(NULL)

      ## Merge per-asset across tiles -- gives one raster per asset
      ## covering the AOI's Sinusoidal extent.
      bands_merged <- lapply(idx$assets, function(an){
        asset_tiles <- lapply(per_tile_bands, function(tb) tb[[an]])
        asset_tiles <- Filter(Negate(is.null), asset_tiles)
        if( length(asset_tiles) == 0L ) return(NULL)
        if( length(asset_tiles) == 1L ) return(asset_tiles[[1L]])
        tryCatch(.mergeRasterTiles(asset_tiles, method = "bilinear"),
                 error = function(e){
                   warning(sprintf(
                     "Sinusoidal asset merge failed for asset %s on date %s; using first tile only. Reason: %s",
                     an, date_str, conditionMessage(e)),
                     call. = FALSE)
                   asset_tiles[[1L]]
                 })
      })
      names(bands_merged) <- idx$assets
      if( any(vapply(bands_merged, is.null, logical(1L))) ) return(NULL)

      idx_rast <- tryCatch({
        stk <- do.call(c, bands_merged)
        .finiteIndex(idx$fun(stk))
      }, error = function(e) NULL)
      if( is.null(idx_rast) ) return(NULL)

      ## Asset hrefs -- preserved from the first tile's STAC item so
      ## downstream can re-fetch via direct asset URL without searching.
      item_assets <- lapply(group[[1L]][["assets"]], function(a) a[["href"]])

      list(date = dt, bands = bands_merged, index = idx_rast,
           item_assets = item_assets)
    })
    scenes <- Filter(Negate(is.null), scenes)
    if( length(scenes) == 0L ) return(NULL)
    ## Deterministic order: production-time grouping is server-order; sort
    ## by date so rs$scenes and downstream anchors are reproducible.
    scenes <- scenes[order(vapply(scenes, function(s) as.numeric(s$date), numeric(1)))]
    return(scenes)
  }

  scenes <- lapply(features, function(f){
    ## datetime is empty-string for time-range items. %||% only fires on
    ## NULL, so coerce empty-string to NULL before falling back.
    dt_str <- f[["properties"]][["datetime"]]
    if( is.null(dt_str) || !nzchar(dt_str) )
      dt_str <- f[["properties"]][["start_datetime"]]
    dt <- tryCatch(as.Date(substr(dt_str, 1, 10)), error = function(e) NULL)
    if( is.null(dt) || is.na(dt) ) return(NULL)

    out <- tryCatch({
      pixel_mask <- NULL
      bands <- lapply(idx$assets, function(an){
        read_band(f[["assets"]][[an]][["href"]], scale = TRUE)
      })
      names(bands) <- idx$assets
      bands <- Filter(Negate(is.null), bands)
      if( length(bands) == 0L ) return(NULL)
      bands <- .harmoniseBands(bands)      # S2 mixes 10 m and 20 m assets
      stk <- do.call(c, bands)

      if( !is.null(cfg$mask_asset) && !is.null(scl_classes) ){
        mb <- read_band(f[["assets"]][[cfg$mask_asset]][["href"]], scale = FALSE)
        if( !is.null(mb) ){
          mb <- terra::resample(mb, stk, method = "near")
          cloud_mask <- if( identical(cfg$mask_kind, "scl") ){
            Reduce(`|`, lapply(scl_classes, function(cls) mb == cls))
          }else{
            bit_sum <- sum(vapply(scl_classes,
                                  function(b) bitwShiftL(1L, as.integer(b)),
                                  integer(1)))
            terra::app(mb, function(v) bitwAnd(as.integer(v), bit_sum) != 0L)
          }
          pixel_mask <- cloud_mask
          stk <- terra::mask(stk, cloud_mask, maskvalue = TRUE)
          ## Mask individual bands so reused-asset paths stay consistent.
          bands <- lapply(bands, function(b)
            terra::mask(b, cloud_mask, maskvalue = TRUE))
        }
      }
      idx_rast <- .finiteIndex(idx$fun(stk))
      item_assets <- lapply(f[["assets"]], function(a) a[["href"]])
      list(bands = bands, index = idx_rast, item_assets = item_assets,
           pixel_mask = pixel_mask)
    }, error = function(e) NULL)

    if( is.null(out) ) return(NULL)
    list(date        = dt,
         bands       = out$bands,
         index       = out$index,
         item_assets = out$item_assets,
         pixel_mask  = out$pixel_mask)
  })
  scenes <- Filter(Negate(is.null), scenes)
  if( length(scenes) == 0L ) return(NULL)

  ## A bbox may span adjacent source tiles. Combine every same-date acquisition
  ## before selecting a cross-date grid so coverage outside the first item is
  ## neither discarded nor represented as a duplicate date.
  scenes <- .mosaicFineScenesByDate(scenes, idx)
  if( length(scenes) == 0L ) return(NULL)

  ## Deterministic order before target-grid construction.
  scenes <- scenes[order(vapply(scenes, function(s) as.numeric(s$date), numeric(1)))]

  ## Build an AOI-wide union template rather than cropping every date to the
  ## first scene's extent. Bands stay in source CRS; only index rasters are
  ## aligned here.
  if( length(scenes) > 1L ){
    template <- terra::rast(.mergeRasterTiles(
      lapply(scenes, `[[`, "index"), method = "bilinear"))
    for( i in seq_along(scenes) ){
      if( !terra::compareGeom(scenes[[i]]$index, template,
                              stopOnError = FALSE) )
        scenes[[i]]$index <- terra::project(scenes[[i]]$index, template,
                                            method = "bilinear")
    }
  }
  return(scenes)
}
