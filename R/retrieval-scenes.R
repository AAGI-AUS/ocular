# =========================================================================
# Retrieval: scenes and public accessors
# =========================================================================

#' Retrieve Landsat or Sentinel-2 scenes
#'
#' Retrieves Landsat or Sentinel-2 scenes through the Microsoft Planetary
#' Computer SpatioTemporal Asset Catalog (STAC) and returns an \code{ocular}
#' object. The supplied point defines the centre of the analysis area. A
#' multi-temporal spectral feature stack is prepared when field boundary
#' delineation requires it.
#'
#' The \code{params} list is stored on the returned object and used by
#' subsequent processing and output functions. Named arguments supplied to
#' those functions override the stored setting for that call only; they do not
#' modify \code{rs$params}.
#'
#' @section Supported indices:
#' Sentinel-2 supports \code{NDVI}, \code{DVI}, \code{kNDVI}, \code{GNDVI},
#' \code{BNDVI}, \code{GBNDVI}, \code{GRNDVI}, \code{NIRv}, \code{NDWI},
#' \code{NDMI}, \code{EVI}, \code{SAVI}, \code{MSAVI}, \code{NDRE},
#' \code{MTCI}, \code{CIre_gitelson}, \code{CIRE}, \code{PSRI}, \code{BSI},
#' \code{SR}, \code{EVI2}, \code{OSAVI}, \code{ARVI}, \code{GARI},
#' \code{CIG}, \code{MNDWI}, \code{NBR}, \code{NBR2}, \code{BAI},
#' \code{MCARI}, and \code{IRECI}.
#'
#' Landsat supports \code{NDVI}, \code{DVI}, \code{kNDVI}, \code{GNDVI},
#' \code{BNDVI}, \code{GBNDVI}, \code{GRNDVI}, \code{NIRv}, \code{NDWI},
#' \code{NDMI}, \code{EVI}, \code{SAVI}, \code{MSAVI}, \code{BSI},
#' \code{SR}, \code{EVI2}, \code{OSAVI}, \code{ARVI}, \code{GARI},
#' \code{CIG}, \code{MNDWI}, \code{NBR}, \code{NBR2}, and \code{BAI}.
#'
#' @param longitude,latitude WGS-84 coordinates at the centre of the analysis
#'   area.
#' @param start_date,end_date Start and end of the retrieval period.
#' @param index_name Index calculated for the output scenes and used by
#'   \code{as_raster()} and \code{as_time_series()}. Field boundary delineation
#'   uses \code{search_index} instead (\code{"EVI2"} by default).
#' @param source Image source. \code{"auto"} selects a source according to
#'   \code{start_date}; alternatives are \code{"sentinel-2"},
#'   \code{"landsat-8"} (Landsat 8/9), and \code{"landsat-5"} (Landsat 5/7).
#' @param x_metres,y_metres Width and height of the analysis area in metres.
#'   At least one is required; if one is \code{NULL}, the other is used for
#'   both dimensions.
#' @param max_cloud_cover Maximum scene-level cloud cover, from 0 to 100,
#'   permitted by the STAC search.
#' @param scl_classes Optional unique non-negative integers to exclude from
#'   pixels: Sentinel-2 SCL class values from 0 to 11, or Landsat
#'   \code{QA_PIXEL} bit positions from 0 to 15. \code{NULL} applies no
#'   per-pixel class or bit mask; \code{max_cloud_cover} still filters scenes
#'   using STAC metadata.
#' @param params Optional complete or partial list of settings accepted by
#'   \code{rs_params()}. \code{NULL} uses the defaults.
#' @param search_index,search_start_date,search_end_date Inline overrides
#'   for the index and date range used for field boundary delineation.
#'   Search dates must lie within \code{start_date} and \code{end_date}; they
#'   do not restrict the dates available to output functions. \code{NULL}
#'   defers to \code{params}.
#' @returns An \code{ocular} object containing the query, retrieved scenes,
#'   and stored parameters.
#' @examples
#' \dontrun{
#' rs <- get_rs(
#'   longitude = 117.8, latitude = -32.0,
#'   start_date = "2021-05-01", end_date = "2021-11-30",
#'   index_name = "NDVI", x_metres = 1000
#' )
#' field <- boundary_delineation(rs)
#' as_time_series(field)
#' }
#' @export
get_rs <- function(longitude,
                   latitude,
                   start_date,
                   end_date,
                   index_name        = "EVI2",
                   source            = "auto",
                   x_metres          = NULL,
                   y_metres          = NULL,
                   max_cloud_cover   = 50,
                   scl_classes       = NULL,
                   params            = NULL,
                   search_index      = NULL,
                   search_start_date = NULL,
                   search_end_date   = NULL){

  ## Normalize params to a fully-populated rs_params list.
  if( is.null(params) ){
    params <- rs_params()
    param_user_set <- attr(params, "user_set", exact = TRUE)
  }else{
    if( !is.list(params) || is.null(names(params)) ||
        anyNA(names(params)) || any(!nzchar(names(params))) ||
        anyDuplicated(names(params)) )
      stop("get_rs(): params must be a uniquely named list (usually from rs_params()).",
           call. = FALSE)
    param_user_set <- .paramUserSet(params, "get_rs()")
    params <- do.call(rs_params, params)
  }
  ## Apply inline formal-arg overrides for fetch-time areas.
  if( !is.null(search_index) ){
    params$search_index <- search_index
    param_user_set <- union(param_user_set, "search_index")
  }
  if( !is.null(search_start_date) ){
    params$search_start_date <- search_start_date
    param_user_set <- union(param_user_set, "search_start_date")
  }
  if( !is.null(search_end_date) ){
    params$search_end_date <- search_end_date
    param_user_set <- union(param_user_set, "search_end_date")
  }
  ## Re-validate after inline overrides.
  params <- do.call(rs_params, params)
  attr(params, "user_set") <- param_user_set

  if( !.isScalarNumber(longitude) || longitude < -180 || longitude > 180 )
    stop("longitude must be a finite numeric in [-180, 180].", call. = FALSE)
  if( !.isScalarNumber(latitude) || latitude < -80 || latitude > 84 )
    stop("latitude must lie in the UTM-supported band [-80, 84].",
         call. = FALSE)
  if( !is.character(index_name) || length(index_name) != 1L ||
      is.na(index_name) || !nzchar(index_name) )
    stop("index_name must be a single non-empty string.", call. = FALSE)
  if( !is.character(source) || length(source) != 1L || is.na(source) )
    stop("source must be one of auto, sentinel-2, landsat-8, or landsat-5.",
         call. = FALSE)
  source <- match.arg(source,
                      c("auto", "sentinel-2", "landsat-8", "landsat-5"))
  if( !.isScalarNumber(max_cloud_cover) ||
      max_cloud_cover < 0 || max_cloud_cover > 100 )
    stop("max_cloud_cover must be a finite numeric in [0, 100].",
         call. = FALSE)
  if( !is.null(scl_classes) ){
    valid_scl <- is.numeric(scl_classes) && length(scl_classes) > 0L &&
      !anyNA(scl_classes) && all(is.finite(scl_classes)) &&
      all(scl_classes >= 0) && all(scl_classes <= .Machine$integer.max) &&
      all(scl_classes == trunc(scl_classes)) &&
      !anyDuplicated(scl_classes)
    if( !valid_scl )
      stop("scl_classes must be NULL or a unique vector of non-negative integers.",
           call. = FALSE)
    scl_classes <- as.integer(scl_classes)
  }

  if( is.null(x_metres) && is.null(y_metres) )
    stop("get_rs() requires at least one of x_metres or y_metres.", call. = FALSE)
  if( is.null(x_metres) ) x_metres <- y_metres
  if( is.null(y_metres) ) y_metres <- x_metres
  if( !.isScalarNumber(x_metres) || !.isScalarNumber(y_metres) ||
      x_metres <= 0 || y_metres <= 0 )
    stop("x_metres and y_metres must be finite positive numeric scalars.",
         call. = FALSE)

  start_date <- .normaliseOptionalDate(start_date, "start_date")
  end_date <- .normaliseOptionalDate(end_date, "end_date")
  if( is.null(start_date) || is.null(end_date) )
    stop("start_date and end_date are required.", call. = FALSE)

  search_index <- params$search_index
  res <- .resolveSource(start_date, search_index, source = source)
  source <- res$source

  ## Validate that index_name (used for output scenes) is also a known
  ## index for the resolved source.
  output_idx <- if( source == "sentinel-2" ) s2_index_list[[index_name]]
  else                         landsat_index_list[[index_name]]
  if( is.null(output_idx) )
    stop("Unknown index_name '", index_name, "' for source '", source, "'.",
         call. = FALSE)
  if( !is.null(scl_classes) ){
    max_class <- if( identical(source, "sentinel-2") ) 11L else 15L
    if( any(scl_classes > max_class) )
      stop("scl_classes must be in [0, ", max_class, "] for source '",
           source, "'.", call. = FALSE)
  }

  ## Output scenes span the full requested date range. Search dates restrict
  ## field boundary delineation, not output dates.
  out_start_d <- suppressWarnings(as.Date(start_date))
  out_end_d   <- suppressWarnings(as.Date(end_date))
  if( is.na(out_start_d) || is.na(out_end_d) )
    stop("Invalid start_date or end_date.", call. = FALSE)
  if( as.numeric(out_end_d - out_start_d) < 2 )
    stop("Date range too short.", call. = FALSE)

  ## Detection range = restricted search window when supplied, else the
  ## full output range. Must lie within the output range.
  det_start_d <- suppressWarnings(as.Date(params$search_start_date %||% start_date))
  det_end_d   <- suppressWarnings(as.Date(params$search_end_date   %||% end_date))
  if( is.na(det_start_d) || is.na(det_end_d) )
    stop("Invalid search_start_date or search_end_date.", call. = FALSE)
  if( det_start_d < out_start_d || det_end_d > out_end_d )
    stop("search_start_date/search_end_date must fall within [start_date, end_date].",
         call. = FALSE)
  if( as.numeric(det_end_d - det_start_d) < 2 )
    stop("Detection date range too short.", call. = FALSE)

  ## Construct the bounding box used by the STAC query.
  fetch_bbox <- point_to_bbox(longitude, latitude, x_metres, y_metres)

  ## Fetch output scenes over the full requested range.
  scenes <- .fetchStac(bbox            = fetch_bbox,
                       start_date      = as.character(out_start_d),
                       end_date        = as.character(out_end_d),
                       index_name      = index_name,
                       max_cloud_cover = max_cloud_cover,
                       scl_classes     = scl_classes,
                       source          = source)
  if( is.null(scenes) || length(scenes) == 0L )
    stop("No usable scenes were found for the requested point, dates, index, ",
         "and cloud filter.", call. = FALSE)

  ## Detection scenes: needed when detection differs from output on the
  ## index dimension (search_index != index_name) OR the date dimension
  ## (a restricted search window was supplied). Build them from the
  ## already-fetched items without a second STAC search. For the index
  ## dimension, .composeSceneFromBands reuses bands already present and
  ## fetches missing assets via .fetchAssets; bands are shared by
  ## reference between $scenes and $internals$detection_scenes --
  ## copy-on-modify keeps memory at one set of bands per date.
  need_detection <- !identical(search_index, index_name) ||
    !is.null(params$search_start_date) ||
    !is.null(params$search_end_date)
  detection_scenes <- if( need_detection ) list() else NULL
  if( need_detection && !is.null(scenes) && length(scenes) > 0L ){
    ## (a) date-restrict to the detection window.
    sc_dates <- do.call(c, lapply(scenes, `[[`, "date"))
    det_src  <- scenes[which(sc_dates >= det_start_d & sc_dates <= det_end_d)]
    ## (b) recompose onto search_index if it differs from index_name.
    if( !identical(search_index, index_name) ){
      search_lookup <- if( source == "sentinel-2" ) s2_index_list[[search_index]]
      else                         landsat_index_list[[search_index]]
      if( is.null(search_lookup) )
        stop("search_index '", search_index, "' is not known for source '",
             source, "'.", call. = FALSE)
      det_src <- lapply(det_src, function(sc){
        tryCatch(.composeSceneFromBands(sc, search_lookup, fetch_bbox, source),
                 error = function(e){
                   warning(sprintf("Could not compose %s for date %s: %s",
                                   search_index, format(sc$date),
                                   conditionMessage(e)), call. = FALSE)
                   NULL
                 })
      })
      det_src <- Filter(Negate(is.null), det_src)
    }
    detection_scenes <- det_src
  }
  if( need_detection && length(detection_scenes) == 0L )
    stop("No usable detection scenes remain for search_index '",
         search_index, "' in the requested detection date range.",
         call. = FALSE)

  rs <- .newOcular()
  rs$spec$point$longitude  <- longitude
  rs$spec$point$latitude   <- latitude
  rs$spec$start_date       <- as.character(out_start_d)
  rs$spec$end_date         <- as.character(out_end_d)
  rs$spec$x_metres         <- x_metres
  rs$spec$y_metres         <- y_metres
  rs$spec$index_name       <- index_name
  rs$spec$max_cloud_cover  <- max_cloud_cover
  rs$spec$scl_classes      <- scl_classes

  rs$geom$bbox             <- fetch_bbox
  rs$geom$source           <- source

  rs$scenes                       <- scenes
  rs$internals$detection_scenes   <- detection_scenes

  rs$params <- params

  ## Decide eager vs lazy feature_stack fetch by inspecting the parent
  ## call. Lazy only when piped directly into a scenes-only consumer
  ## (these never read feature_stack). Anything else -- terminal call,
  ## unknown parent, feature-stack consumer, standalone-shape wrapper --
  ## fetches eagerly. Conservative: when in doubt, eager.
  parent_call <- if( length(sys.calls()) >= 2L ) sys.call(-1L) else NULL
  parent_fn <- if( is.null(parent_call) ) NA_character_
  else if( is.symbol(parent_call[[1L]]) ) as.character(parent_call[[1L]])
  else NA_character_
  scenes_only_consumers <- c("as_raster", "as_time_series", "add_modis")
  need_features <- !(parent_fn %in% scenes_only_consumers)

  if( need_features ) rs <- .fetchFeatureStack(rs)
  return(rs)
}

#' Apply an index lookup's compute function to a named band list (internal)
#'
#' Builds a multi-layer SpatRaster from \code{bands} in the order required
#' by \code{lookup_entry$assets}, then applies \code{lookup_entry$fun}.
#' Errors if any required band is missing from \code{bands}.
#'
#' @noRd
.computeIndexFromBands <- function(bands, lookup_entry){

  needed <- lookup_entry$assets
  missing <- setdiff(needed, names(bands))
  if( length(missing) )
    stop(".computeIndexFromBands(): missing band(s) ",
         paste(missing, collapse = ", "), call. = FALSE)
  ordered <- .harmoniseBands(bands[needed])   # reused/refetched bands may differ
  stk <- do.call(c, ordered)
  return(.finiteIndex(lookup_entry$fun(stk)))
}

#' Bring a band list onto one grid before stacking (internal)
#'
#' Sentinel-2 indices mix 10 m assets (B02-B04, B08) with 20 m ones (B05-B07,
#' B8A, B11, B12): NDMI is B08 + B11, NDRE is B08 + B05, MNDWI is B03 + B11.
#' \code{terra::c()} requires a common geometry, so an unharmonised stack throws
#' and -- because the callers wrap the stack in \code{tryCatch(..., NULL)} --
#' every scene is dropped and the fetch reports "No scenes available
#' for feature_stack fetch". Landsat and MODIS are single-resolution per source
#' and are returned untouched.
#'
#' Resamples onto the finest band's grid (bilinear), preserving the spatial
#' detail delineation depends on; it adds no information to the coarse bands.
#'
#' @param bands Named list of SpatRasters.
#' @returns The list, all entries sharing the finest band's geometry.
#' @noRd
.harmoniseBands <- function(bands){

  if( length(bands) < 2L ) return(bands)
  cell_area <- vapply(bands, function(b){
    rr <- abs(terra::res(b))
    if( length(rr) != 2L || any(!is.finite(rr)) ) Inf else prod(rr)
  }, numeric(1L))
  if( all(!is.finite(cell_area)) ) return(bands)
  ref <- bands[[which.min(cell_area)]]
  return(lapply(bands, function(b){
    if( isTRUE(terra::compareGeom(b, ref, stopOnError = FALSE)) ) b
    else terra::project(b, ref, method = "bilinear")
  }))
}

#' Build a scene entry with a different index from existing bands (internal)
#'
#' Given a scene that already carries bands and a desired
#' \code{lookup_entry}, returns the scene with the index recomputed for the
#' new lookup, reusing bands already present and fetching only what is
#' missing via \code{.fetchAssets}. Used by \code{get_rs()} when
#' \code{search_index} differs from \code{index_name}, so detection runs on
#' one index while the output scenes carry another -- without a second STAC
#' search.
#'
#' @param scene A scene list (date, bands, index, item_assets, ...).
#' @param lookup_entry Index entry from \code{s2_index_list} / \code{landsat_index_list}.
#' @param bbox Search bbox (sf bbox object).
#' @param source One of "sentinel-2", "landsat-8", or "landsat-5".
#' @returns The scene with \code{$index} recomputed and any newly fetched
#'   bands appended to \code{$bands}.
#' @noRd
.composeSceneFromBands <- function(scene, lookup_entry, bbox, source){

  needed  <- lookup_entry$assets
  have    <- (scene$bands %||% list())[intersect(needed,
                                                 names(scene$bands %||% list()))]
  missing <- setdiff(needed, names(have))
  fetched <- list()
  if( length(missing) > 0L ){
    fetched <- .fetchAssets(item_assets        = scene$item_assets,
                            needed_asset_names = missing,
                            bbox               = bbox,
                            source             = source)
    still <- setdiff(missing, names(fetched))
    if( length(still) )
      stop("could not fetch band(s) ", paste(still, collapse = ", "),
           " for scene ", format(scene$date),
           " (asset unavailable, no AOI overlap, or expired signature).",
           call. = FALSE)
    ## Bands fetched for an alternate detection index must receive the same
    ## SCL/QA mask as the original output-index bands. Otherwise cloudy pixels
    ## can reappear simply because search_index uses a disjoint asset set.
    if( !is.null(scene$pixel_mask) ){
      fetched <- lapply(fetched, function(b){
        pm <- scene$pixel_mask
        if( !terra::compareGeom(pm, b, stopOnError = FALSE) )
          pm <- terra::project(pm, b, method = "near")
        terra::mask(b, pm, maskvalue = TRUE)
      })
    }
  }
  bands       <- c(have, fetched)[needed]
  scene$index <- .computeIndexFromBands(bands, lookup_entry)
  new_names   <- setdiff(names(fetched), names(scene$bands %||% list()))
  scene$bands <- c(scene$bands %||% list(), fetched[new_names])
  return(scene)
}

#' Direct asset fetch via stored item URIs (internal)
#'
#' Retrieves specific band rasters for a single STAC item without
#' re-issuing a STAC search. Used when an rs already carries
#' \code{item_assets} on each scene from the initial \code{.fetchStac}
#' call but a downstream consumer needs additional bands not currently
#' fetched. Honours the source-specific \code{cfg} (scaling, valid range)
#' and crops each asset read to the AOI window.
#'
#' @param item_assets Named list mapping asset name to href (STAC URL).
#' @param needed_asset_names Character vector of asset names to fetch.
#' @param bbox Search bbox (sf bbox object); each read is cropped to it.
#' @param source A source key supported by \code{.stacCfg()}:
#'   \code{"sentinel-2"}, \code{"landsat-8"}, \code{"landsat-5"}, or
#'   \code{"mcd43a4"}.
#' @returns Named list of SpatRasters keyed by asset name.
#' @noRd
.fetchAssets <- function(item_assets, needed_asset_names, bbox, source){

  cfg <- .stacCfg(source)
  aoi_sfc  <- sf::st_as_sfc(bbox)

  read_one <- function(href){
    tryCatch({
      b <- terra::rast(paste0("/vsicurl/", href))
      terra::scoff(b) <- NULL
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
      if( !is.null(cfg$valid_range) )
        b <- terra::ifel(b < cfg$valid_range[1] | b > cfg$valid_range[2], NA, b)
      if( !is.null(cfg$scale_factor) )
        b <- b * cfg$scale_factor + (cfg$scale_offset %||% 0)
      b
    }, error = function(e) NULL)
  }
  out <- lapply(needed_asset_names, function(an){
    hrefs <- unique(as.character(unlist(item_assets[[an]], use.names = FALSE)))
    hrefs <- hrefs[!is.na(hrefs) & nzchar(hrefs)]
    if( length(hrefs) == 0L ) return(NULL)
    .mergeRasterTiles(lapply(hrefs, read_one), method = "bilinear")
  })
  names(out) <- needed_asset_names
  return(Filter(Negate(is.null), out))
}

#' Fetch and assemble the detection feature_stack onto an rs (internal)
#'
#' Reads all required fetch state from \code{rs} (bbox, dates,
#' \code{params$search_index}, source, utm_epsg, centre_v, x_metres,
#' y_metres, max_cloud_cover, scl_classes). Populates
#' \code{rs$internals$feature_stack}, \code{rs$geom$nr}, \code{rs$geom$nc},
#' \code{rs$geom$centre_rc}, \code{rs$geom$pixel_size_m}, \code{rs$state$decided_mat}.
#' Idempotent: returns \code{rs} unchanged if feature_stack is already
#' populated.
#'
#' @noRd
.fetchFeatureStack <- function(rs){
  if( !is.null(rs$internals$feature_stack) ) return(rs)

  search_index <- rs$params$search_index %||% rs_params()$search_index
  res <- .resolveSource(rs$spec$start_date, search_index,
                        source = rs$geom$source)
  pixel_size_m <- res$idx$res

  ## Window breaks span the detection range (search window when set), not
  ## the full output range: detection_scenes only populate that sub-range,
  ## so breaks over the full range would leave empty windows and trip the
  ## no-scenes guard below.
  start_d <- as.Date(rs$params$search_start_date %||% rs$spec$start_date)
  end_d   <- as.Date(rs$params$search_end_date   %||% rs$spec$end_date)
  ## Window count uses one window per (2 * pixel_size_m) days, bounded to one
  ## through three windows. This engineering rule produces longer composite
  ## windows for coarser sensors: about 40 days at 10 m and 60 days at 30 m.
  num_windows <- min(3L, trunc(as.numeric(end_d - start_d) / (pixel_size_m * 2) + 1))
  if( num_windows < 1L ) num_windows <- 1L
  breaks <- seq(start_d, end_d, length.out = num_windows + 1L)

  ## Source for window-binning: detection_scenes if get_rs() pre-fetched
  ## search_index material (when search_index != index_name); otherwise
  ## the user-facing scenes (when search_index == index_name and the
  ## index area already carries the right raster). Reuse stored item_assets
  ## rather than issuing another STAC search.
  detection_scenes <- rs$internals$detection_scenes %||% rs$scenes
  if( length(detection_scenes) == 0L )
    stop("No scenes available for feature_stack fetch.", call. = FALSE)
  scene_dates <- do.call(c, lapply(detection_scenes, `[[`, "date"))

  composites <- lapply(seq_len(num_windows), function(w){
    ## Half-open membership [break_w, break_{w+1}) for all but the last
    ## window (closed on the right): a scene landing exactly on an interior
    ## break joins one window, not both adjacent composites.
    keep <- if( w < num_windows )
      which(scene_dates >= breaks[w] & scene_dates <  breaks[w + 1L])
    else
      which(scene_dates >= breaks[w] & scene_dates <= breaks[w + 1L])
    if( length(keep) == 0L ) return(NULL)
    rasters  <- lapply(detection_scenes[keep], `[[`, "index")
    template <- rasters[[1L]]
    rasters  <- lapply(rasters, function(r){
      if( !terra::compareGeom(template, r, stopOnError = FALSE) )
        terra::resample(r, template, method = "bilinear") else r
    })
    stk <- terra::rast(rasters)
    composite <- if( terra::nlyr(stk) > 1L )
      terra::app(stk, max, na.rm = TRUE) else stk
    ## max(..., na.rm = TRUE) is -Inf when every scene is missing at a cell.
    ## Restore missingness before bilinear projection can spread non-finite
    ## values into neighbouring cells.
    .finiteIndex(composite)
  })
  if( any(vapply(composites, is.null, logical(1))) )
    stop("Fetch failed: one or more windows produced no scenes. Try widening x_metres/y_metres, the date range, or max_cloud_cover.",
         call. = FALSE)

  comps_utm <- lapply(composites, function(r){
    tryCatch(terra::project(r, rs_utm(rs), method = "bilinear"),
             error = function(e) NULL)
  })
  if( any(vapply(comps_utm, is.null, logical(1))) )
    stop("UTM projection failed.", call. = FALSE)

  ref_r <- comps_utm[[1L]]
  for( i in seq_along(comps_utm) ){
    if( !terra::compareGeom(comps_utm[[i]], ref_r, stopOnError = FALSE) )
      comps_utm[[i]] <- terra::resample(comps_utm[[i]], ref_r, method = "bilinear")
  }
  feature_stack <- Reduce(c, comps_utm)
  names(feature_stack) <- paste0("w", 1:num_windows)

  ## Re-crop feature_stack to a tight box around the analysis centre in
  ## target UTM. After the WGS->source-UTM transform, snap-out crop,
  ## target UTM projection, and multi-window resample, feature_stack
  ## tends to be a bit larger than x_metres x y_metres and the centre
  ## isn't necessarily at the geometric centre. Both shift the reserve
  ## geometry's position relative to the seed and change BFS dynamics,
  ## so we trim back here.
  centre_xy <- .centreXY(rs)
  crop_extent <- terra::ext(
    centre_xy[["x"]] - rs$spec$x_metres / 2,
    centre_xy[["x"]] + rs$spec$x_metres / 2,
    centre_xy[["y"]] - rs$spec$y_metres / 2,
    centre_xy[["y"]] + rs$spec$y_metres / 2)
  feature_stack <- terra::crop(feature_stack, crop_extent, snap = "out")

  cell <- terra::cellFromXY(feature_stack[[1L]],
                            matrix(centre_xy, nrow = 1L))
  if( is.na(cell) )
    stop("Centre falls outside fetched grid.", call. = FALSE)
  centre_rc <- c(as.integer(terra::rowFromCell(feature_stack[[1L]], cell)),
                 as.integer(terra::colFromCell(feature_stack[[1L]], cell)))

  nr <- terra::nrow(feature_stack)
  nc <- terra::ncol(feature_stack)

  rs$internals$feature_stack <- feature_stack
  rs$geom$centre_rc          <- centre_rc
  rs$geom$pixel_size_m       <- as.numeric(terra::res(feature_stack)[1L])
  rs$geom$nr                 <- nr
  rs$geom$nc                 <- nc
  rs$state$decided_mat       <- matrix(FALSE, nr, nc)
  ## Run calibration once the detection feature stack and geometry exist.
  ## The eager/lazy need_features decision and the idempotence guard above
  ## bound it to one attempt. Failures degrade to a warning-bearing record,
  ## so retrieval itself never fails solely because calibration did.
  if( is.null(rs$internals$calibration) ){
    rs$internals$calibration <- tryCatch(
      .calibrateScene(rs),
      error = function(e) list(
        source   = "self",
        error    = conditionMessage(e),
        warnings = paste0(".calibrateScene() failed: ", conditionMessage(e))))
  }

  return(rs)
}

#' Analysed-area extent in target UTM (internal)
#'
#' Centred on \code{rs$spec$point}, sized exactly \code{x_metres} x
#' \code{y_metres}. Used by \code{as_raster} / \code{as_time_series} to
#' crop the output composite/scenes to the user's requested area --
#' independent of the detection feature_stack so output anchors to the
#' source-scene grid via \code{terra::crop(..., snap = "out")}.
#'
#' @noRd
.analysedExt <- function(rs){

  centre_xy <- .centreXY(rs)
  return(terra::ext(centre_xy[["x"]] - rs$spec$x_metres / 2,
                    centre_xy[["x"]] + rs$spec$x_metres / 2,
                    centre_xy[["y"]] - rs$spec$y_metres / 2,
                    centre_xy[["y"]] + rs$spec$y_metres / 2))
}

#' Resolve effective params at a call site (internal)
#'
#' Returns the fully-populated params list a stage should run with.
#' Merge order: \code{rs$params} (auto-applied base) -> \code{params}
#' (call-site overlay, overriding each area it sets) -> \code{overrides}
#' (inline named args from the stage, override anything they set).
#' \code{rs$params} is not mutated.
#'
#' Inline overrides are validated against both the rs_params name set and
#' the value constraints enforced by \code{rs_params()}; unknown or invalid
#' values fail before a stage mutates the object.
#'
#' @param use_calibration Whether to apply the current calibration baseline.
#'   A caller preparing calibration inputs sets this to \code{FALSE}, then
#'   resolves once more after recalibration.
#' @noRd
.resolveParams <- function(rs, params, overrides, use_calibration = TRUE){

  defaults <- rs_params()
  base <- defaults
  resolved_user_set <- character(0L)
  stored <- rs$params
  if( !is.null(stored) && length(stored) > 0L ){
    if( !is.list(stored) || is.null(names(stored)) ||
        anyNA(names(stored)) || any(!nzchar(names(stored))) ||
        anyDuplicated(names(stored)) )
      stop(".resolveParams(): rs$params must be a uniquely named list.",
           call. = FALSE)
    unknown_stored <- setdiff(names(stored), names(base))
    if( length(unknown_stored) )
      stop("Unknown stored params: ", paste(unknown_stored, collapse = ", "),
           call. = FALSE)
    stored_user_set <- attr(stored, "user_set", exact = TRUE)
    if( is.null(stored_user_set) )
      stop(".resolveParams(): rs$params is missing explicit parameter ",
           "ownership metadata. Recreate the ocular object with get_rs().",
           call. = FALSE)
    stored_user_set <- .paramUserSet(stored, ".resolveParams()")
    stored_user_set <- intersect(stored_user_set, names(stored))
    resolved_user_set <- union(resolved_user_set, stored_user_set)
    base[names(stored)] <- stored
  }
  retrieval_params <- base

  ## Calibration baselines are the lowest-priority overlay. Explicitly owned
  ## values, direct non-default mutations, params, and inline overrides take
  ## precedence. Baselines are not written back into rs$params.
  cal_record <- rs$internals$calibration
  cal <- if( isTRUE(use_calibration) &&
             identical(cal_record$schema_version,
                       .CALIBRATION_SCHEMA_VERSION) )
    cal_record$baselines
  else NULL
  if( !is.null(cal) && length(cal) > 0L ){
    for( nm in intersect(names(cal), names(base)) ){
      if( is.null(cal[[nm]]) ) next   ## no baseline was derived
      if( .paramIsAtDefault(rs, nm) ) base[[nm]] <- cal[[nm]]
    }
  }

  if( !is.null(params) ){
    if( !is.list(params) || is.null(names(params)) ||
        anyNA(names(params)) || any(!nzchar(names(params))) ||
        anyDuplicated(names(params)) )
      stop(".resolveParams(): params must be a named list (from rs_params()).",
           call. = FALSE)
    unknown <- setdiff(names(params), names(base))
    if( length(unknown) )
      stop("Unknown params: ", paste(unknown, collapse = ", "), call. = FALSE)
    params_user_set <- .paramUserSet(params, ".resolveParams()")
    resolved_user_set <- union(
      resolved_user_set, intersect(params_user_set, names(params)))
    base[names(params)] <- params
  }
  if( length(overrides) > 0L ){
    if( is.null(names(overrides)) || anyNA(names(overrides)) ||
        any(!nzchar(names(overrides))) || anyDuplicated(names(overrides)) )
      stop(".resolveParams(): inline overrides must be uniquely named.",
           call. = FALSE)
    unknown <- setdiff(names(overrides), names(base))
    if( length(unknown) )
      stop("Unknown params: ", paste(unknown, collapse = ", "), call. = FALSE)
    non_null <- overrides[!vapply(overrides, is.null, logical(1L))]
    if( length(non_null) > 0L ){
      resolved_user_set <- union(resolved_user_set, names(non_null))
      base[names(non_null)] <- non_null
    }
  }
  ## Centralised value validation is important because R function formals do
  ## not enforce types. It also normalises integer-valued and date parameters.
  resolved <- do.call(rs_params, base)
  attr(resolved, "user_set") <- resolved_user_set
  fetch_only <- c("search_index", "search_start_date", "search_end_date")
  changed_fetch <- fetch_only[!vapply(fetch_only, function(nm)
    identical(resolved[[nm]], retrieval_params[[nm]]), logical(1L))]
  if( length(changed_fetch) )
    stop("The fetch-time parameter(s) ", paste(changed_fetch, collapse = ", "),
         " cannot be changed after get_rs(). Re-run get_rs() with those ",
         "settings so the detection scenes and feature stack match.",
         call. = FALSE)
  return(resolved)
}

#' Apply call-local calibration inputs to an rs (internal)
#'
#' Recomputes calibration against a temporary copy carrying the effective
#' call-local params when calibration-relevant inputs differ. The returned rs
#' receives only the new calibration record: its stored \code{rs$params} are
#' never mutated. A failed recalibration warns and preserves the prior record.
#'
#' @noRd
.applyScenePriors <- function(rs, effective_params){
  if( is.null(rs$internals$feature_stack) ) return(rs)
  if( !is.list(effective_params) )
    stop(".applyScenePriors(): effective_params must be a params list.",
         call. = FALSE)
  cal <- rs$internals$calibration
  changed <- is.null(cal) ||
    !identical(cal$schema_version, .CALIBRATION_SCHEMA_VERSION) ||
    !identical(cal$strictness, effective_params$strictness) ||
    !identical(cal$area_threshold, effective_params$area_threshold) ||
    !identical(cal$index_name, effective_params$search_index)
  if( !changed ) return(rs)

  rs_for_cal <- rs
  rs_for_cal$params <- effective_params
  updated <- tryCatch(.calibrateScene(rs_for_cal), error = identity)
  if( inherits(updated, "error") ||
      (is.list(updated) && !is.null(updated$error)) ){
    msg <- if( inherits(updated, "error") ) conditionMessage(updated)
           else updated$error
    warning("Call-local recalibration failed; preserving the prior calibration: ",
            msg, call. = FALSE)
    return(rs)
  }
  rs$internals$calibration <- updated
  rs
}
