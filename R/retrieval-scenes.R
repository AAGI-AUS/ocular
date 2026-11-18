# =========================================================================
# Retrieval: scenes and public accessors
# =========================================================================

#' Retrieve Landsat or Sentinel-2 scenes
#'
#' Retrieves Landsat or Sentinel-2 scenes through the Microsoft Planetary
#' Computer SpatioTemporal Asset Catalog (STAC) and returns an \code{ocular}
#' object. The supplied point defines the centre of the analysis area. Field
#' boundary delineation uses spectral features from one or more periods.
#'
#' Settings supplied through \code{params} are stored on the returned object.
#' See \code{\link{rs_params}} for changing them in later calls.
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
#' @param start_date,end_date Start and end of the retrieval period, both
#'   inclusive.
#' @param index_name Index calculated for the output scenes and used by
#'   \code{as_raster()} and \code{as_time_series()}. Field boundary delineation
#'   uses \code{search_index} instead (\code{"EVI2"} by default).
#' @param source Image source. \code{"auto"} selects a source according to
#'   \code{start_date}. Alternatives are \code{"sentinel-2"},
#'   \code{"landsat-8"} (Landsat 8/9), and \code{"landsat-5"} (Landsat 5/7).
#' @param x_metres,y_metres Width and height of the analysis area in metres.
#'   At least one is required. If one is \code{NULL}, the other is used for
#'   both dimensions.
#' @param max_cloud_cover Maximum scene-level cloud cover, from 0 to 100,
#'   permitted by the STAC search.
#' @param scl_classes Optional unique non-negative integers to exclude from
#'   pixels: Sentinel-2 SCL class values from 0 to 11, or Landsat
#'   \code{QA_PIXEL} bit positions from 0 to 15. \code{NULL} applies no
#'   per-pixel class or bit mask. \code{max_cloud_cover} still filters scenes
#'   using STAC metadata.
#' @param params Optional complete or partial list of settings accepted by
#'   \code{rs_params()}. \code{NULL} uses the defaults.
#' @param search_index,search_start_date,search_end_date Inline overrides
#'   for the index and date range used for field boundary delineation.
#'   Search dates must lie within \code{start_date} and \code{end_date}. They
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
                   index_name = "EVI2",
                   source = "auto",
                   x_metres = NULL,
                   y_metres = NULL,
                   max_cloud_cover = 50,
                   scl_classes = NULL,
                   params = NULL,
                   search_index = NULL,
                   search_start_date = NULL,
                   search_end_date = NULL){

  # Normalise params to a complete rs_params list.
  if( is.null(params) ){
    params <- rs_params()
    param_user_set <- attr(params, "user_set", exact = TRUE)
  }else{
    if( !is.list(params) || is.null(names(params)) ||
        anyNA(names(params)) || any(!nzchar(names(params))) ||
        anyDuplicated(names(params)) ){
      stop("get_rs(): params must be a uniquely named list.",
           call. = FALSE)
    }
    param_user_set <- .paramUserSet(params, "get_rs()")
    params <- do.call(rs_params, params)
  }
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
  params <- do.call(rs_params, params)
  attr(params, "user_set") <- param_user_set

  if( !.isScalarNumber(longitude) || longitude < -180 || longitude > 180 ){
    stop("longitude must be a finite numeric in [-180, 180].", call. = FALSE)
  }
  if( !.isScalarNumber(latitude) || latitude < -80 || latitude > 84 ){
    stop("latitude must lie in the UTM-supported band [-80, 84].",
         call. = FALSE)
  }
  if( !is.character(index_name) || length(index_name) != 1L ||
      is.na(index_name) || !nzchar(index_name) ){
    stop("index_name must be a single non-empty string.", call. = FALSE)
  }
  if( !is.character(source) || length(source) != 1L || is.na(source) ){
    stop("source must be one of auto, sentinel-2, landsat-8, or landsat-5.",
         call. = FALSE)
  }
  source <- match.arg(source,
                      c("auto", "sentinel-2", "landsat-8", "landsat-5"))
  if( !.isScalarNumber(max_cloud_cover) ||
      max_cloud_cover < 0 || max_cloud_cover > 100 ){
    stop("max_cloud_cover must be a finite numeric in [0, 100].",
         call. = FALSE)
  }
  if( !is.null(scl_classes) ){
    valid_scl <- is.numeric(scl_classes) && length(scl_classes) > 0L &&
      !anyNA(scl_classes) && all(is.finite(scl_classes)) &&
      all(scl_classes >= 0) && all(scl_classes <= .Machine$integer.max) &&
      all(scl_classes == trunc(scl_classes)) &&
      !anyDuplicated(scl_classes)
    if( !valid_scl ){
      stop("scl_classes must be NULL or a unique vector of non-negative integers.",
           call. = FALSE)
    }
    scl_classes <- as.integer(scl_classes)
  }

  if( is.null(x_metres) && is.null(y_metres) ){
    stop("get_rs() requires at least one of x_metres or y_metres.", call. = FALSE)
  }
  if( is.null(x_metres) ){
    x_metres <- y_metres
  }
  if( is.null(y_metres) ){
    y_metres <- x_metres
  }
  if( !.isScalarNumber(x_metres) || !.isScalarNumber(y_metres) ||
      x_metres <= 0 || y_metres <= 0 ){
    stop("x_metres and y_metres must be finite positive numeric scalars.",
         call. = FALSE)
  }

  start_date <- .normaliseOptionalDate(start_date, "start_date")
  end_date <- .normaliseOptionalDate(end_date, "end_date")
  if( is.null(start_date) || is.null(end_date) ){
    stop("start_date and end_date are required.", call. = FALSE)
  }

  search_index <- params$search_index
  # EVI2 is available from every supported sensor.
  res <- .resolveSource(start_date, "EVI2", source = source)
  source <- res$source

  output_idx <- if( identical(source, "sentinel-2") ){
    s2_index_list[[index_name]]
  }else{
    landsat_index_list[[index_name]]
  }
  if( is.null(output_idx) ){
    stop("Unknown index_name '", index_name, "' for source '", source, "'.",
         call. = FALSE)
  }
  search_lookup <- if( identical(source, "sentinel-2") ){
    s2_index_list[[search_index]]
  }else{
    landsat_index_list[[search_index]]
  }
  if( is.null(search_lookup) ){
    stop("search_index '", search_index, "' is not known for source '",
         source, "'.", call. = FALSE)
  }
  if( !is.null(scl_classes) ){
    max_class <- if( identical(source, "sentinel-2") ) { 11L } else { 15L }
    if( any(scl_classes > max_class) ){
      stop("scl_classes must be in [0, ", max_class, "] for source '",
           source, "'.", call. = FALSE)
    }
  }

  # Search dates limit delineation. Output scenes still span the requested period.
  out_start_d <- suppressWarnings(as.Date(start_date))
  out_end_d <- suppressWarnings(as.Date(end_date))
  if( is.na(out_start_d) || is.na(out_end_d) ){
    stop("Invalid start_date or end_date.", call. = FALSE)
  }
  if( as.numeric(out_end_d - out_start_d) < 2 ){
    stop("Date range too short.", call. = FALSE)
  }

  det_start_d <- suppressWarnings(as.Date(params$search_start_date %||% start_date))
  det_end_d <- suppressWarnings(as.Date(params$search_end_date   %||% end_date))
  if( is.na(det_start_d) || is.na(det_end_d) ){
    stop("Invalid search_start_date or search_end_date.", call. = FALSE)
  }
  if( det_start_d < out_start_d || det_end_d > out_end_d ){
    stop("search_start_date/search_end_date must fall within [start_date, end_date].",
         call. = FALSE)
  }
  if( as.numeric(det_end_d - det_start_d) < 2 ){
    stop("Detection date range too short.", call. = FALSE)
  }

  fetch_bbox <- point_to_bbox(longitude, latitude, x_metres, y_metres)

  scenes <- .fetchStac(bbox = fetch_bbox,
                       start_date = as.character(out_start_d),
                       end_date = as.character(out_end_d),
                       index_name = index_name,
                       max_cloud_cover = max_cloud_cover,
                       scl_classes = scl_classes,
                       source = source,
                       s2_processing = params$s2_processing)
  if( is.null(scenes) || length(scenes) == 0L ){
    stop("No usable scenes were found for the requested point, dates, index, ",
         "and cloud filter.", call. = FALSE)
  }

  # Defer features only for recognised output calls that do not need them.
  # Other calls prepare features now. add_modis() is recognised by its bare name.
  parent_call <- if( length(sys.calls()) >= 2L ) { sys.call(-1L) } else { NULL }
  parent_fn <- if( is.null(parent_call) ){
    NA_character_
  }else if( is.symbol(parent_call[[1L]]) ){
    as.character(parent_call[[1L]])
  }else{
    NA_character_
  }
  # Compare the calling function so qualified and indirect output calls also work.
  parent_fun <- if( is.null(parent_call) ) { NULL } else { sys.function(-1L) }
  output_caller <- identical(parent_fun, as_raster) ||
    identical(parent_fun, as_time_series)
  scenes_only_caller <- output_caller ||
    identical(parent_fn, "add_modis")
  need_features <- !scenes_only_caller

  # A different search index or period needs a separate detection series.
  # Reuse the fetched items and existing bands, fetching only missing bands.
  need_detection <- !output_caller &&
    (!identical(search_index, index_name) ||
       !is.null(params$search_start_date) ||
       !is.null(params$search_end_date))
  detection_scenes <- if( need_detection ) { list() } else { NULL }
  if( need_detection && !is.null(scenes) && length(scenes) > 0L ){
    sc_dates <- do.call(c, lapply(scenes, `[[`, "date"))
    det_src <- scenes[which(sc_dates >= det_start_d & sc_dates <= det_end_d)]
    if( !identical(search_index, index_name) ){
      det_src <- lapply(det_src, function(sc){
        return(tryCatch(.composeSceneFromBands(sc, search_lookup, fetch_bbox, source),
                        error = function(e){
                          warning(sprintf("Could not compose %s for date %s: %s",
                                          search_index, format(sc$date),
                                          conditionMessage(e)), call. = FALSE)
                          return(NULL)
                        }))
      })
      det_src <- Filter(Negate(is.null), det_src)
    }
    detection_scenes <- det_src
  }
  if( need_detection && length(detection_scenes) == 0L ){
    stop("No usable detection scenes remain for search_index '",
         search_index, "' in the requested detection date range.",
         call. = FALSE)
  }

  rs <- .newOcular()
  rs$spec$point$longitude <- longitude
  rs$spec$point$latitude <- latitude
  rs$spec$start_date <- as.character(out_start_d)
  rs$spec$end_date <- as.character(out_end_d)
  rs$spec$x_metres <- x_metres
  rs$spec$y_metres <- y_metres
  rs$spec$index_name <- index_name
  rs$spec$max_cloud_cover <- max_cloud_cover
  rs$spec$scl_classes <- scl_classes

  rs$geom$bbox <- fetch_bbox
  rs$geom$source <- source

  rs$scenes <- scenes
  rs$internals$detection_scenes <- detection_scenes

  rs$params <- params

  if( need_features ){
    rs <- .fetchFeatureStack(rs)
  }
  return(rs)
}

#' Calculate an index from its required bands
#'
#' Order the bands using \code{lookup_entry$assets}, then apply its function.
#' Stop if a required band is missing.
#' @noRd
.computeIndexFromBands <- function(bands, lookup_entry){

  needed <- lookup_entry$assets
  missing <- setdiff(needed, names(bands))
  if( length(missing) ){
    stop(".computeIndexFromBands(): missing band(s) ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  ordered <- .harmoniseBands(bands[needed])   # May include additional bands.
  stk <- do.call(c, unname(ordered))
  names(stk) <- names(ordered)
  return(.finiteIndex(lookup_entry$fun(stk)))
}

#' Align bands before stacking them
#'
#' Sentinel-2 indices can combine bands at 10 m and 20 m resolution.
#' \code{terra::c()} requires a common geometry, so resample mismatched bands
#' to the finest band's grid using bilinear interpolation. Matching bands stay
#' unchanged. Resampling adds no detail to the coarser bands.
#'
#' @param bands Named list of \code{SpatRaster} objects.
#' @returns The bands on a common grid.
#' @noRd
.harmoniseBands <- function(bands){

  if( length(bands) < 2L ){
    return(bands)
  }
  cell_area <- vapply(bands, function(b){
    rr <- abs(terra::res(b))
    if( length(rr) != 2L || any(!is.finite(rr)) ) { Inf } else { prod(rr) }
  }, numeric(1L))
  if( all(!is.finite(cell_area)) ){
    return(bands)
  }
  ref <- bands[[which.min(cell_area)]]
  return(lapply(bands, function(b){
    if( isTRUE(terra::compareGeom(b, ref, stopOnError = FALSE)) ){
      b
    }else{
      terra::project(b, ref, method = "bilinear")
    }
  }))
}

#' Recalculate a scene with a different index
#'
#' Reuse its bands and fetch missing ones from the stored asset URLs, avoiding
#' a second STAC search. This separates the detection index from output indices.
#'
#' @param scene A scene with bands and stored asset URLs.
#' @param lookup_entry The required index's entry in a sensor lookup.
#' @param bbox Bounding box for asset reads.
#' @param source One of "sentinel-2", "landsat-8" or "landsat-5".
#' @returns The scene with the new index and any additional bands.
#' @noRd
.composeSceneFromBands <- function(scene, lookup_entry, bbox, source){

  needed <- lookup_entry$assets
  have <- (scene$bands %||% list())[intersect(needed,
                                              names(scene$bands %||% list()))]
  missing <- setdiff(needed, names(have))
  fetched <- list()
  if( length(missing) > 0L ){
    fetched <- .fetchAssets(item_assets = scene$item_assets,
                            needed_asset_names = missing,
                            bbox = bbox,
                            source = source,
                            dn_offset = scene$dn_offset)
    still <- setdiff(missing, names(fetched))
    if( length(still) ){
      stop("could not fetch band(s) ", paste(still, collapse = ", "),
           " for scene ", format(scene$date),
           " (asset unavailable, no AOI overlap, or expired signature).",
           call. = FALSE)
    }
    # Apply the same cloud mask to newly fetched bands, so changing the search
    # index cannot reintroduce masked pixels.
    if( !is.null(scene$pixel_mask) ){
      fetched <- lapply(fetched, function(b){
        pm <- scene$pixel_mask
        if( !terra::compareGeom(pm, b, stopOnError = FALSE) ){
          pm <- terra::project(pm, b, method = "near")
        }
        return(terra::mask(b, pm, maskvalue = TRUE))
      })
    }
  }
  bands <- c(have, fetched)[needed]
  scene$index <- .computeIndexFromBands(bands, lookup_entry)
  new_names <- setdiff(names(fetched), names(scene$bands %||% list()))
  scene$bands <- c(scene$bands %||% list(), fetched[new_names])
  return(scene)
}

#' Read additional bands from stored asset URLs
#'
#' Use the sensor's scaling and valid range, and crop each read to the query box.
#' No new STAC search is needed.
#'
#' @param item_assets Named list of asset URLs from a STAC item.
#' @param needed_asset_names Names of the required bands.
#' @param bbox Bounding box to crop each read.
#' @param source A source supported by \code{.stacCfg()}.
#' @param dn_offset Optional offsets, named by asset and ordered by URL.
#'   Subtract these before range checks and scaling. Unlisted assets have no offset.
#' @returns The band rasters, named by asset.
#' @noRd
.fetchAssets <- function(item_assets, needed_asset_names, bbox, source,
                         dn_offset = NULL){

  cfg <- .stacCfg(source)
  aoi_sfc <- sf::st_as_sfc(bbox)

  read_one <- function(href, shift){

    return(tryCatch({
      b <- terra::rast(paste0("/vsicurl/", href))
      terra::scoff(b) <- NULL
      aoi_in_crs <- sf::st_transform(aoi_sfc, terra::crs(b))
      aoi_ext <- terra::ext(sf::st_bbox(aoi_in_crs))
      te <- terra::ext(b)
      overlaps <- terra::xmin(aoi_ext) < terra::xmax(te) &&
        terra::xmax(aoi_ext) > terra::xmin(te) &&
        terra::ymin(aoi_ext) < terra::ymax(te) &&
        terra::ymax(aoi_ext) > terra::ymin(te)
      if( !overlaps ){
        return(NULL)
      }
      b <- terra::crop(b, aoi_ext, snap = "out")
      if( terra::ncell(b) == 0 ){
        return(NULL)
      }
      if( shift != 0 ){
        b <- b - shift
      }
      if( !is.null(cfg$valid_range) ){
        b <- terra::ifel(b < cfg$valid_range[1] | b > cfg$valid_range[2], NA, b)
      }
      if( !is.null(cfg$scale_factor) ){
        b <- b * cfg$scale_factor + (cfg$scale_offset %||% 0)
      }
      b
    }, error = function(e){
      return(NULL)
    }))
  }
  out <- lapply(needed_asset_names, function(an){
    hrefs <- as.character(unlist(item_assets[[an]], use.names = FALSE))
    shifts <- dn_offset[[an]] %||% rep(0, length(hrefs))
    keep <- !duplicated(hrefs) & !is.na(hrefs) & nzchar(hrefs)
    hrefs <- hrefs[keep]
    shifts <- shifts[keep]
    if( length(hrefs) == 0L ){
      return(NULL)
    }
    return(.mergeRasterTiles(mapply(read_one, hrefs, shifts, SIMPLIFY = FALSE,
                                    USE.NAMES = FALSE),
                             method = "bilinear"))
  })
  names(out) <- needed_asset_names
  return(Filter(Negate(is.null), out))
}

#' Prepare detection features and their grid
#'
#' Use the stored query and scenes to build temporal composites, the analysis
#' grid and its calibration. Return \code{rs} unchanged if features already exist.
#' @noRd
.fetchFeatureStack <- function(rs){

  if( !is.null(rs$internals$feature_stack) ){
    return(rs)
  }

  search_index <- rs$params$search_index %||% rs_params()$search_index
  res <- .resolveSource(rs$spec$start_date, search_index,
                        source = rs$geom$source)
  pixel_size_m <- res$idx$res

  # Split the detection period, not the full output period, to avoid empty windows.
  start_d <- as.Date(rs$params$search_start_date %||% rs$spec$start_date)
  end_d <- as.Date(rs$params$search_end_date   %||% rs$spec$end_date)
  # Use up to three equal windows, with more windows for longer periods and
  # fewer for coarser pixels.
  num_windows <- min(3L, trunc(as.numeric(end_d - start_d) / (pixel_size_m * 2) + 1))
  if( num_windows < 1L ){
    num_windows <- 1L
  }
  breaks <- seq(start_d, end_d, length.out = num_windows + 1L)

  # Reuse the separate detection series when present, otherwise use output scenes.
  detection_scenes <- rs$internals$detection_scenes %||% rs$scenes
  if( length(detection_scenes) == 0L ){
    stop("Cannot prepare ", search_index, " detection features: no scenes available from ",
         start_d, " to ", end_d, " (0 of ", num_windows, " windows have scenes).",
         call. = FALSE)
  }
  scene_dates <- do.call(c, lapply(detection_scenes, `[[`, "date"))

  # Put all dates on a common grid covering their combined extent before
  # compositing, including scenes recomputed in different source grids.
  first_index <- detection_scenes[[1L]]$index
  same_geom <- vapply(detection_scenes, function(sc){
    return(terra::compareGeom(first_index, sc$index, stopOnError = FALSE))
  }, logical(1L))
  if( !all(same_geom) ){
    template <- terra::rast(.mergeRasterTiles(
      lapply(detection_scenes, `[[`, "index"), method = "bilinear"))
    detection_scenes <- lapply(detection_scenes, function(sc){
      if( !terra::compareGeom(sc$index, template, stopOnError = FALSE) ){
        sc$index <- terra::project(sc$index, template, method = "bilinear")
      }
      return(sc)
    })
  }

  composites <- lapply(seq_len(num_windows), function(w){
    # A scene on a window boundary belongs to the later window. Include the
    # last date in the final window.
    keep <- if( w < num_windows ){
      which(scene_dates >= breaks[w] & scene_dates <  breaks[w + 1L])
    }else{
      which(scene_dates >= breaks[w] & scene_dates <= breaks[w + 1L])
    }
    if( length(keep) == 0L ){
      return(NULL)
    }
    rasters <- lapply(detection_scenes[keep], `[[`, "index")
    template <- rasters[[1L]]
    rasters <- lapply(rasters, function(r){
      if( !terra::compareGeom(template, r, stopOnError = FALSE) ){
        terra::resample(r, template, method = "bilinear")
      }else{
        r
      }
    })
    stk <- terra::rast(rasters)
    composite <- if( terra::nlyr(stk) > 1L ) { terra::app(stk, max, na.rm = TRUE) } else { stk }
    # max() returns -Inf for an all-missing cell. Restore NA before interpolation
    # can spread it to neighbouring cells.
    return(.finiteIndex(composite))
  })
  empty_windows <- which(vapply(composites, is.null, logical(1)))
  if( length(empty_windows) ){
    intervals <- vapply(empty_windows, function(w){
      # Scene dates are whole days. Report the inclusive dates tested by this window.
      first_day <- ceiling(as.numeric(breaks[w]))
      last_day <- if( w < num_windows ){
        ceiling(as.numeric(breaks[w + 1L])) - 1
      }else{
        floor(as.numeric(breaks[w + 1L]))
      }
      return(paste0("window ", w, " (", as.Date(first_day, origin = "1970-01-01"),
                    " to ", as.Date(last_day, origin = "1970-01-01"), ")"))
    }, character(1L))
    stop("Cannot prepare ", search_index, " detection features: no scenes in ",
         paste(intervals, collapse = ", "), ". ", num_windows - length(empty_windows),
         " of ", num_windows, " windows have scenes.", call. = FALSE)
  }

  comps_utm <- lapply(composites, function(r){
    return(tryCatch(terra::project(r, rs_utm(rs), method = "bilinear"),
                    error = function(e){
                      return(NULL)
                    }))
  })
  if( any(vapply(comps_utm, is.null, logical(1))) ){
    stop("UTM projection failed.", call. = FALSE)
  }

  ref_r <- comps_utm[[1L]]
  for( i in seq_along(comps_utm) ){
    if( !terra::compareGeom(comps_utm[[i]], ref_r, stopOnError = FALSE) ){
      comps_utm[[i]] <- terra::resample(comps_utm[[i]], ref_r, method = "bilinear")
    }
  }
  feature_stack <- Reduce(c, comps_utm)
  names(feature_stack) <- paste0("w", 1:num_windows)

  # Projection and resampling can enlarge or shift the feature grid.
  # Crop it back around the query point before flood fill.
  centre_xy <- .centreXY(rs)
  crop_extent <- terra::ext(
    centre_xy[["x"]] - rs$spec$x_metres / 2,
    centre_xy[["x"]] + rs$spec$x_metres / 2,
    centre_xy[["y"]] - rs$spec$y_metres / 2,
    centre_xy[["y"]] + rs$spec$y_metres / 2)
  feature_stack <- terra::crop(feature_stack, crop_extent, snap = "out")

  cell <- terra::cellFromXY(feature_stack[[1L]],
                            matrix(centre_xy, nrow = 1L))
  if( is.na(cell) ){
    stop("Centre falls outside fetched grid.", call. = FALSE)
  }
  centre_rc <- c(as.integer(terra::rowFromCell(feature_stack[[1L]], cell)),
                 as.integer(terra::colFromCell(feature_stack[[1L]], cell)))

  nr <- terra::nrow(feature_stack)
  nc <- terra::ncol(feature_stack)

  rs$internals$feature_stack <- feature_stack
  rs$geom$centre_rc <- centre_rc
  rs$geom$pixel_size_m <- as.numeric(terra::res(feature_stack)[1L])
  rs$geom$nr <- nr
  rs$geom$nc <- nc
  rs$state$decided_mat <- matrix(FALSE, nr, nc)
  # Calibration needs the completed feature grid. Record a failure without
  # failing the scene retrieval.
  if( is.null(rs$internals$calibration) ){
    rs$internals$calibration <- tryCatch(
      .calibrateScene(rs),
      error = function(e){
        return(list(
          source = "self",
          error = conditionMessage(e),
          warnings = paste0(".calibrateScene() failed: ", conditionMessage(e))))
      })
  }

  return(rs)
}

#' Construct the requested analysis extent in UTM
#'
#' Output functions crop scenes to this box independently of the detection grid.
#' @noRd
.analysedExt <- function(rs){

  centre_xy <- .centreXY(rs)
  return(terra::ext(centre_xy[["x"]] - rs$spec$x_metres / 2,
                    centre_xy[["x"]] + rs$spec$x_metres / 2,
                    centre_xy[["y"]] - rs$spec$y_metres / 2,
                    centre_xy[["y"]] + rs$spec$y_metres / 2))
}

#' Resolve settings for one stage call
#'
#' Start with stored settings and eligible calibration defaults. Apply explicit
#' settings from \code{params}, then named \code{overrides}. Validate the result
#' without changing \code{rs$params}.
#'
#' @param use_calibration Whether to use the current calibration. Set this to
#'   \code{FALSE} when preparing inputs for recalibration, then resolve again.
#' @param caller Function name to use in parameter errors.
#' @noRd
.resolveParams <- function(rs, params, overrides, use_calibration = TRUE,
                            caller = ".resolveParams()"){

  defaults <- rs_params()
  base <- defaults
  resolved_user_set <- character(0L)
  stored <- rs$params
  if( !is.null(stored) && length(stored) > 0L ){
    if( !is.list(stored) || is.null(names(stored)) ||
        anyNA(names(stored)) || any(!nzchar(names(stored))) ||
        anyDuplicated(names(stored)) ){
      stop(caller, ": rs$params must be a uniquely named list.",
           call. = FALSE)
    }
    unknown_stored <- setdiff(names(stored), names(base))
    if( length(unknown_stored) ){
      stop(caller, ": unknown stored params: ", paste(unknown_stored, collapse = ", "),
           call. = FALSE)
    }
    stored_user_set <- attr(stored, "user_set", exact = TRUE)
    if( is.null(stored_user_set) ){
      stop(caller, ": rs$params is missing explicit parameter ",
           "ownership metadata. Recreate the ocular object with get_rs().",
           call. = FALSE)
    }
    stored_user_set <- .paramUserSet(stored, caller, "rs$params")
    stored_user_set <- intersect(stored_user_set, names(stored))
    resolved_user_set <- union(resolved_user_set, stored_user_set)
    base[names(stored)] <- stored
  }
  retrieval_params <- base

  # Use calibration only where the user has not chosen or changed a setting.
  cal_record <- rs$internals$calibration
  cal <- if( isTRUE(use_calibration) &&
             identical(cal_record$schema_version,
                       .CALIBRATION_SCHEMA_VERSION) ){
    cal_record$baselines
  }else{
    NULL
  }
  if( !is.null(cal) && length(cal) > 0L ){
    for( nm in intersect(names(cal), names(base)) ){
      if( is.null(cal[[nm]]) ){
        next
      }
      if( .paramIsAtDefault(rs, nm) ){
        base[[nm]] <- cal[[nm]]
      }
    }
  }

  fetch_only <- c("search_index", "search_start_date", "search_end_date", "s2_processing")
  if( !is.null(params) ){
    if( !is.list(params) || is.null(names(params)) ||
        anyNA(names(params)) || any(!nzchar(names(params))) ||
        anyDuplicated(names(params)) ){
      stop(caller, ": params must be a uniquely named list.",
           call. = FALSE)
    }
    unknown <- setdiff(names(params), names(base))
    if( length(unknown) ){
      stop(caller, ": unknown params: ", paste(unknown, collapse = ", "), ".",
           call. = FALSE)
    }
    params_user_set <- .paramUserSet(params, caller)
    resolved_user_set <- union(
      resolved_user_set, intersect(params_user_set, names(params)))
    # Defaults filled in by rs_params() do not override stored or calibrated values.
    overlay <- names(params)
    unchosen <- vapply(overlay, function(nm){
      return(!(nm %in% params_user_set) && identical(params[[nm]], defaults[[nm]]))
    }, logical(1L))
    overlay <- overlay[!unchosen]
    base[overlay] <- params[overlay]
  }
  if( length(overrides) > 0L ){
    if( is.null(names(overrides)) || anyNA(names(overrides)) ||
        any(!nzchar(names(overrides))) || anyDuplicated(names(overrides)) ){
      stop(caller, ": inline overrides must be uniquely named.",
           call. = FALSE)
    }
    unknown <- setdiff(names(overrides), names(base))
    if( length(unknown) ){
      stop(caller, ": unknown params: ", paste(unknown, collapse = ", "), ".",
           call. = FALSE)
    }
    non_null <- overrides[!vapply(overrides, is.null, logical(1L))]
    if( length(non_null) > 0L ){
      resolved_user_set <- union(resolved_user_set, names(non_null))
      base[names(non_null)] <- non_null
    }
  }
  # Validate and normalise the merged settings in one place.
  resolved <- do.call(rs_params, base)
  attr(resolved, "user_set") <- resolved_user_set
  changed_fetch <- fetch_only[!vapply(fetch_only, function(nm){
    return(identical(resolved[[nm]], retrieval_params[[nm]]))
  }, logical(1L))]
  if( length(changed_fetch) ){
    stop(caller, ": the fetch-time parameter(s) ", paste(changed_fetch, collapse = ", "),
         " cannot be changed after get_rs(). Re-run get_rs() with those ",
         "settings so the detection scenes and feature stack match.",
         call. = FALSE)
  }
  return(resolved)
}

#' Start a shared record of failed calibration attempts
#' @noRd
.beginCalibrationInvocation <- function(rs){

  context <- attr(rs, ".calibration_context", exact = TRUE)
  owner <- !is.environment(context) || !isTRUE(context$active)
  if( owner ){
    context <- new.env(parent = emptyenv())
    context$active <- TRUE
    context$failures <- list()
    attr(rs, ".calibration_context") <- context
  }
  return(list(rs = rs, context = context, owner = owner))
}

#' Remove the temporary retry record from the result
#' @noRd
.endCalibrationInvocation <- function(result, invocation){

  if( invocation$owner && is_rs(result) ){
    attr(result, ".calibration_context") <- NULL
  }
  return(result)
}

#' Copy the values that determine calibration
#' @noRd
.calibrationInputSnapshot <- function(rs, effective_params){

  fstack <- rs$internals$feature_stack
  polygon <- rs$geom$ftw_prior$polygon
  if( inherits(polygon, "SpatVector") ){
    polygon <- list(type = terra::geomtype(polygon),
                    geometry = terra::geom(polygon), crs = terra::crs(polygon))
  }else if( inherits(polygon, c("sf", "sfc")) ){
    polygon <- sf::st_geometry(polygon)
  }
  return(list(
    schema = .CALIBRATION_SCHEMA_VERSION,
    params = effective_params[c("strictness", "area_threshold", "search_index")],
    index_name = rs$spec$index_name,
    point = rs$spec$point,
    centre = rs$geom$centre_rc,
    grid = c(rs$geom$nr, rs$geom$nc),
    feature_dim = dim(fstack),
    feature_extent = as.vector(terra::ext(fstack)),
    feature_crs = terra::crs(fstack),
    feature_values = terra::values(fstack),
    polygon = polygon,
    calibration = rs$internals$calibration
  ))
}

#' Refresh calibration for the settings used by this call
#'
#' Recalculate when the schema or settings differ, without changing stored
#' parameters. On failure, warn and keep the existing calibration.
#'
#' Within one workflow, remember failures so later stages do not repeat the
#' same attempt. Changed inputs or a new call allow another attempt. If the
#' inputs cannot be copied for comparison, each stage retries.
#' @noRd
.applyScenePriors <- function(rs, effective_params){

  if( is.null(rs$internals$feature_stack) ){
    return(rs)
  }
  if( !is.list(effective_params) ){
    stop(".applyScenePriors(): effective_params must be a params list.",
         call. = FALSE)
  }
  cal <- rs$internals$calibration
  changed <- is.null(cal) ||
    !identical(cal$schema_version, .CALIBRATION_SCHEMA_VERSION) ||
    !identical(cal$strictness, effective_params$strictness) ||
    !identical(cal$area_threshold, effective_params$area_threshold) ||
    !identical(cal$index_name, effective_params$search_index)
  if( !changed ){
    return(rs)
  }

  context <- attr(rs, ".calibration_context", exact = TRUE)
  inputs <- NULL
  if( is.environment(context) && isTRUE(context$active) ){
    inputs <- tryCatch(suppressWarnings(.calibrationInputSnapshot(rs, effective_params)),
                       error = function(e){
                         return(NULL)
                       })
    if( !is.null(inputs) && any(vapply(context$failures, function(previous){
      return(identical(previous, inputs))
    }, logical(1L))) ){
      return(rs)
    }
  }

  rs_for_cal <- rs
  rs_for_cal$params <- effective_params
  updated <- tryCatch(.calibrateScene(rs_for_cal), error = identity)
  if( inherits(updated, "error") ||
      (is.list(updated) && !is.null(updated$error)) ){
    msg <- if( inherits(updated, "error") ) { conditionMessage(updated) } else { updated$error }
    if( !is.null(inputs) ){
      context$failures[[length(context$failures) + 1L]] <- inputs
    }
    warning("Recalibration for this call failed. The prior calibration was preserved: ",
            msg, call. = FALSE)
    return(rs)
  }
  rs$internals$calibration <- updated
  return(rs)
}
