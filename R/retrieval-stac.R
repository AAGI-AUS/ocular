# =========================================================================
# Retrieval: source resolution and STAC fetch
# =========================================================================

#' Look up sensor metadata and reading settings
#' @noRd
.stacCfg <- function(source){

  # Tiling controls search padding and mosaics. Mask type controls how pixel
  # quality values are read.
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
                # MCD43A4 provides daily reflectance at 500 m. Fusion combines its bands
                # with Landsat bands before calculating an index.
                "mcd43a4" = list(
                  endpoint = .STAC_ENDPOINT, collection = "modis-43A4-061",
                  cloud_key = NULL, platforms = NULL,
                  scale_factor = 0.0001, scale_offset = 0,
                  valid_range = c(0, 32766), mask_asset = NULL,
                  grid = list(tiling = "sinusoidal"), mask_kind = NULL),
                stop("Unknown source: ", source, call. = FALSE)))
}

#' Digital-number offset of a Sentinel-2 item
#'
#' From processing baseline 04.00 the digital numbers of a Sentinel-2 item
#' carry an offset of 1000, so that reflectance is (DN - 1000) x 0.0001.
#'
#' @param feature A STAC item.
#' @returns 1000 or 0, or \code{NA} when the baseline cannot be read.
#' @noRd
.s2Offset <- function(feature){

  baseline <- feature[["properties"]][["s2:processing_baseline"]]
  if( is.null(baseline) || length(baseline) != 1L ){
    return(NA_real_)
  }
  baseline <- suppressWarnings(as.numeric(baseline))
  if( is.na(baseline) ){
    return(NA_real_)
  }
  offset <- if( baseline >= 4 ) { 1000 } else { 0 }
  return(offset)
}

#' Normalise acquisition times without rounding fractional seconds
#' @noRd
.s2AcquisitionTime <- function(x){

  if( is.na(x) ){
    return(NA_character_)
  }
  pattern <- paste0("^([0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2})",
                     "([.][0-9]+)?(Z|[+-][0-9]{2}:[0-9]{2})$")
  parts <- regmatches(x, regexec(pattern, x))[[1L]]
  if( length(parts) == 0L ){
    return(NA_character_)
  }
  whole <- as.POSIXct(parts[[2L]], format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
  if( is.na(whole) ||
      !identical(format(whole, "%Y-%m-%dT%H:%M:%S", tz = "UTC"), parts[[2L]]) ){
    return(NA_character_)
  }
  offset <- 0
  zone <- parts[[4L]]
  if( zone != "Z" ){
    hours <- as.numeric(substr(zone, 2L, 3L))
    minutes <- as.numeric(substr(zone, 5L, 6L))
    if( hours > 23 || minutes > 59 ){
      return(NA_character_)
    }
    offset <- (hours * 60 + minutes) * 60
    if( substr(zone, 1L, 1L) == "-" ){
      offset <- -offset
    }
  }
  fraction <- sub("0+$", "", substring(parts[[3L]], 2L))
  fraction <- if( nzchar(fraction) ) { paste0(".", fraction) } else { "" }
  return(paste0(format(whole - offset, "%Y-%m-%dT%H:%M:%S", tz = "UTC"),
                fraction, "Z"))
}

#' Select a processing version for each Sentinel-2 acquisition
#'
#' Datastrip sensing times distinguish complementary deliveries of one tile.
#' Keep all items for mosaicking when identity or baseline metadata is incomplete.
#' @noRd
.selectS2Processing <- function(features, s2_processing){

  text_one <- function(x){

    if( length(x) != 1L || is.na(x) || !nzchar(as.character(x)) ){
      return(NA_character_)
    }
    return(as.character(x))
  }
  metadata <- lapply(features, function(f){
    p <- f[["properties"]]
    identity <- vapply(c("platform", "datetime", "s2:mgrs_tile",
                         "sat:relative_orbit", "s2:datastrip_id"), function(nm){
      return(text_one(p[[nm]]))
    }, character(1L))
    identity[[2L]] <- .s2AcquisitionTime(identity[[2L]])
    strip <- identity[[5L]]
    if( is.na(strip) || !grepl("_S[0-9]{8}T[0-9]{6}_N", strip) ){
      identity[[5L]] <- NA_character_
    }else{
      identity[[5L]] <- sub(".*_S([0-9]{8}T[0-9]{6})_N.*", "\\1", strip)
    }
    baseline <- text_one(p[["s2:processing_baseline"]])
    version <- c(NA_real_, NA_real_)
    if( !is.na(baseline) && grepl("^[0-9]+[.][0-9]+$", baseline) ){
      version <- as.numeric(strsplit(baseline, ".", fixed = TRUE)[[1L]])
    }
    generated <- .s2AcquisitionTime(text_one(p[["s2:generation_time"]]))
    generation <- sub("Z$", "", generated)
    key <- if( anyNA(identity) ) { NA_character_ } else { paste(identity, collapse = "|") }
    return(list(key = key, baseline = baseline, version = version, generation = generation,
                id = text_one(f[["id"]])))
  })
  keys <- vapply(metadata, `[[`, character(1L), "key")
  groups <- split(which(!is.na(keys)), keys[!is.na(keys)])
  groups <- Filter(function(g){
    return(length(g) > 1L)
  }, groups)
  unidentified <- sum(is.na(keys))
  if( length(groups) == 0L && unidentified == 0L ){
    return(features)
  }
  keep <- rep(TRUE, length(features))
  baselines <- character(0L)
  unranked <- 0L
  direction <- if( identical(s2_processing, "highest_baseline") ) { -1 } else { 1 }
  for( group in groups ){
    m <- metadata[group]
    versions <- do.call(rbind, lapply(m, `[[`, "version"))
    if( any(!is.finite(versions)) ){
      unranked <- unranked + 1L
      next
    }
    generation <- vapply(m, `[[`, character(1L), "generation")
    generation_rank <- match(generation, sort(unique(generation), method = "radix"))
    ids <- vapply(m, `[[`, character(1L), "id")
    ranked <- order(direction * versions[, 1L], direction * versions[, 2L],
                    is.na(generation), direction * generation_rank, ids,
                    na.last = TRUE, method = "radix")
    selected <- group[ranked[[1L]]]
    keep[group] <- FALSE
    keep[selected] <- TRUE
    baselines <- c(baselines, metadata[[selected]]$baseline)
  }
  selected_count <- length(groups) - unranked
  label <- if( direction < 0 ) { "highest" } else { "lowest" }
  alternative <- if( direction < 0 ) { "lowest_baseline" } else { "highest_baseline" }
  notice <- character(0L)
  if( selected_count > 0L ){
    notice <- sprintf(paste0("%d Sentinel-2 acquisition(s) have multiple processings. ",
                            "The %s available processing baseline was kept ",
                            "(baseline %s). Set s2_processing = \"%s\" in rs_params() ",
                            "to select the other baseline rule."),
                      selected_count, label, paste(sort(unique(baselines)), collapse = ", "),
                      alternative)
  }
  if( unranked > 0L ){
    notice <- c(notice, sprintf(paste0("%d duplicate Sentinel-2 acquisition group(s) were ",
                                      "retained without processing selection because a ",
                                      "baseline could not be read."), unranked))
  }
  if( unidentified > 0L ){
    notice <- c(notice, sprintf(paste0("%d Sentinel-2 item(s) were retained without ",
                                      "processing selection because acquisition identity ",
                                      "was incomplete."), unidentified))
  }
  warning(paste(notice, collapse = " "), call. = FALSE)
  return(features[keep])
}

#' Resolve the sensor, index and pixel size for a query
#' @noRd
.resolveSource <- function(start_date, index_name, source = "auto"){

  d <- as.Date(start_date)
  if( is.na(d) ){
    stop("start_date could not be parsed: '", start_date, "'", call. = FALSE)
  }
  source <- match.arg(source,
                      c("auto", "sentinel-2", "landsat-8", "landsat-5"))
  if( identical(source, "auto") ){
    source <- if( d >= .S2_START ){
      "sentinel-2"
    }else if( d >= .L8_START ) { "landsat-8" } else { "landsat-5" }
    if( identical(source, "landsat-8") ){
      message("Sentinel-2 not available before ",
              format(.S2_START, "%Y-%m-%d"), ". Landsat 8/9 used.")
    }
    if( identical(source, "landsat-5") ){
      message("Landsat 8 unavailable before ",
              format(.L8_START, "%Y-%m-%d"), ". Landsat 5/7 used.")
    }
  }
  if( identical(source, "sentinel-2") ){
    if( d < .S2_START ){
      stop("Sentinel-2 is unavailable before ", format(.S2_START), ".",
           call. = FALSE)
    }
    idx <- s2_index_list[[index_name]]
    if( is.null(idx) ){
      stop("Unknown index_name '", index_name, "'. Available: ",
           paste(names(s2_index_list), collapse = ", "), call. = FALSE)
    }
    return(list(source = "sentinel-2", cfg = .stacCfg("sentinel-2"), idx = idx))
  }
  if( identical(source, "landsat-8") ){
    if( d < .L8_START ){
      stop("Landsat 8/9 is unavailable before ", format(.L8_START), ".",
           call. = FALSE)
    }
    idx <- landsat_index_list[[index_name]]
    if( is.null(idx) ){
      stop("Unknown index_name '", index_name, "'. Available: ",
           paste(names(landsat_index_list), collapse = ", "), call. = FALSE)
    }
    return(list(source = "landsat-8", cfg = .stacCfg("landsat-8"), idx = idx))
  }
  idx <- landsat_index_list[[index_name]]
  if( is.null(idx) ){
    stop("Unknown index_name '", index_name, "'.", call. = FALSE)
  }
  return(list(source = "landsat-5", cfg = .stacCfg("landsat-5"), idx = idx))
}

#' Merge adjacent raster tiles, reprojecting only when their CRSs differ
#' @noRd
.mergeRasterTiles <- function(rasters, method = "bilinear"){

  rasters <- Filter(Negate(is.null), rasters)
  if( length(rasters) == 0L ){
    return(NULL)
  }
  if( length(rasters) == 1L ){
    return(rasters[[1L]])
  }
  ref_crs <- terra::crs(rasters[[1L]])
  aligned <- lapply(rasters, function(r){
    if( identical(terra::crs(r), ref_crs) ){
      r
    }else{
      terra::project(r, ref_crs, method = method)
    }
  })
  merged <- do.call(terra::merge,
                    c(unname(aligned), list(algo = 3, method = method,
                                    na.rm = TRUE)))
  # Read the cropped values now so the object does not depend on a temporary VRT.
  return(merged * 1)
}

#' Combine fine sensor tiles from the same date
#' @noRd
.mosaicFineScenesByDate <- function(scenes, lookup_entry){

  if( length(scenes) < 2L ){
    return(scenes)
  }
  keys <- vapply(scenes, function(sc){
    return(as.character(sc$date))
  }, character(1L))
  groups <- split(scenes, keys)
  out <- lapply(groups, function(group){
    bands <- lapply(lookup_entry$assets, function(an){
      return(.mergeRasterTiles(lapply(group, function(sc){
        return(sc$bands[[an]])
      })))
    })
    names(bands) <- lookup_entry$assets
    if( any(vapply(bands, is.null, logical(1L))) ){
      return(NULL)
    }

    masks <- lapply(group, `[[`, "pixel_mask")
    masks <- Filter(Negate(is.null), masks)
    pixel_mask <- if( length(masks) == 0L ){
      NULL
    }else{
      .mergeRasterTiles(masks, method = "near")
    }

    asset_names <- unique(unlist(lapply(group, function(sc){
      return(names(sc$item_assets))
    }), use.names = FALSE))
    item_assets <- lapply(asset_names, function(an){
      return(unique(unlist(lapply(group, function(sc){
        return(sc$item_assets[[an]])
      }), use.names = FALSE)))
    })
    names(item_assets) <- asset_names

    scene <- list(date = group[[1L]]$date,
                  bands = bands,
                  index = .computeIndexFromBands(bands, lookup_entry),
                  item_assets = item_assets,
                  pixel_mask = pixel_mask)

    # Keep each URL's offset alongside it when duplicate URLs are removed.
    has_offset <- vapply(group, function(sc){
      return(!is.null(sc$dn_offset))
    }, logical(1L))
    if( any(has_offset) ){
      offset_names <- unique(unlist(lapply(group[has_offset], function(sc){
        return(names(sc$dn_offset))
      }), use.names = FALSE))
      dn_offset <- lapply(offset_names, function(an){
        hrefs <- unlist(lapply(group, function(sc){
          return(sc$item_assets[[an]])
        }), use.names = FALSE)
        shifts <- unlist(lapply(group, function(sc){
          return(sc$dn_offset[[an]] %||% rep(0, length(sc$item_assets[[an]])))
        }), use.names = FALSE)
        return(shifts[!duplicated(hrefs)])
      })
      names(dn_offset) <- offset_names
      scene$dn_offset <- dn_offset
    }
    return(scene)
  })
  return(Filter(Negate(is.null), out))
}

#' Recognise transient STAC failures
#'
#' Use HTTP status or curl class when available. rstac 1.0.1 drops these fields,
#' so recognise its specific message formats too. Leave unknown errors final.
#' @noRd
.isTransientStacError <- function(e){

  status <- e$status_code
  if( is.null(status) && is.list(e$response) ){
    status <- e$response$status_code
  }
  if( is.null(status) ){
    http_class <- grep("^http_[0-9]{3}$", class(e), value = TRUE)
    if( length(http_class) > 0L ){
      status <- as.integer(sub("^http_", "", http_class[[1L]]))
    }
  }
  msg <- conditionMessage(e)
  if( is.null(status) && grepl("^HTTP status '[0-9]{3}'[.]", msg) ){
    status <- as.integer(substr(msg, 14L, 16L))
  }
  if( length(status) == 1L && !is.na(status) ){
    return(status %in% c(408L, 429L, 500L, 502L, 503L, 504L))
  }
  transient_curl <- c("curl_error_couldnt_resolve_proxy",
                      "curl_error_couldnt_resolve_host",
                      "curl_error_couldnt_connect",
                      "curl_error_operation_timedout",
                      "curl_error_partial_file", "curl_error_got_nothing",
                      "curl_error_send_error", "curl_error_recv_error")
  if( any(class(e) %in% transient_curl) ){
    return(TRUE)
  }
  wrapper <- "^Error while requesting '[^\n]*'[.] *\n"
  if( !grepl(wrapper, msg) ){
    return(FALSE)
  }
  detail <- sub(wrapper, "", msg)
  return(grepl(paste0("^(Timeout was reached|Could not resolve host(name)?|",
                      "Could not resolve proxy|Could(n't| not) connect to server|",
                      "Failure when receiving data from the peer|",
                      "Failed sending data to the peer|",
                      "Server returned nothing [(]no headers, no data[)]|",
                      "Transferred a partial file)( \\[|:)"), detail))
}

#' Query STAC items with retries and sign their assets
#'
#' Filter dates and bounds with CQL2, adding cloud cover where the source supports
#' it. Retry failed requests with increasing delays.
#' @noRd
.fetchStacItems <- function(stac_url, collection_id, bbox_vec,
                            start_date, end_date,
                            max_cloud_cover, has_cloud,
                            retry_max = 10L,
                            initial_pause = 2,
                            max_pause = 60){

  # batch_rs() handles its own retry schedule. Its option lasts only through retrieval.
  if( isTRUE(getOption("ocular.batch_single_attempt")) ){
    retry_max <- 1L
  }
  bbox_geo <- rstac::cql2_bbox_as_geojson(bbox_vec)
  # A bare end date means midnight to the service. Include the rest of that day.
  dt_int <- rstac::cql2_interval(
    start_date, paste0(format(as.Date(end_date)), "T23:59:59.999999Z"))
  signer <- rstac::sign_planetary_computer()

  attempt <- 0L
  repeat{
    res <- tryCatch({
      q <- if( has_cloud ){
        rstac::stac(stac_url, force_version = "1.0.0") |>
          rstac::ext_filter(
            collection == {{ collection_id }} &&
              t_intersects(datetime, {{ dt_int }}) &&
              s_intersects(geometry, {{ bbox_geo }}) &&
              `eo:cloud_cover` <= {{ max_cloud_cover }})
      }else{
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
    }, error = function(e){
      return(list(items = NULL, err = e))
    })

    if( !is.null(res$items) ){
      return(res$items)
    }
    attempt <- attempt + 1L
    if( attempt >= retry_max ){
      break
    }
    Sys.sleep(min(max_pause, initial_pause * 2^(attempt - 1L)))
  }
  msg <- paste0("STAC query failed after ", retry_max,
                if( retry_max == 1L ) { " attempt: " } else { " attempts: " },
                .stripUrlQuery(conditionMessage(res$err)))
  stop(errorCondition(msg, class = c("ocular_stac_error", "simpleError"), call = NULL,
                      retryable = .isTransientStacError(res$err)))
}

#' Remove URL query strings before displaying text
#'
#' Asset signatures must not appear in messages.
#' @noRd
.stripUrlQuery <- function(x){

  stripped <- gsub("([[:alpha:]][[:alnum:]+.-]*://[^?[:space:]]*)\\?[^[:space:]'\"]*",
                   "\\1", x)
  return(stripped)
}

#' Retrieve scenes from STAC
#'
#' Read and scale the required bands, then compute each scene's index.
#' Align Landsat and Sentinel-2 index rasters on a grid covering all scenes.
#' Keep their bands on the original grids.
#' @noRd
.fetchStac <- function(bbox,
                       start_date,
                       end_date,
                       index_name = "EVI2",
                       max_cloud_cover = 50,
                       scl_classes = NULL,
                       source = NULL,
                       s2_processing = "highest_baseline"){

  if( !.isScalarNumber(max_cloud_cover) ||
      max_cloud_cover < 0 || max_cloud_cover > 100 ){
    stop(".fetchStac(): max_cloud_cover must be in [0, 100].",
         call. = FALSE)
  }
  if( !inherits(bbox, "bbox") || anyNA(bbox) ||
      any(!is.finite(unclass(bbox))) ){
    stop(".fetchStac(): bbox must be a finite sf bounding box.",
         call. = FALSE)
  }
  if( !requireNamespace("rstac", quietly = TRUE) ){
    stop("rstac package required.", call. = FALSE)
  }

  res <- if( is.null(source) ){
    .resolveSource(start_date, index_name)
  }else{
    list(source = source, cfg = .stacCfg(source),
         idx = if( source == "mcd43a4" ){
           mcd43a4_index_list[[index_name]]
         }else if( grepl("^landsat", source) ){
           landsat_index_list[[index_name]]
         }else{
           s2_index_list[[index_name]]
         })
  }
  cfg <- res$cfg
  idx <- res$idx
  if( is.null(idx) ){
    stop("Unknown index_name '", index_name, "' for source '", res$source, "'.",
         call. = FALSE)
  }

  is_sinusoidal <- identical(cfg$grid$tiling, "sinusoidal")

  # A small query near a MODIS tile seam can miss an adjacent tile because
  # its geographic footprint differs from the sinusoidal grid. Pad the search
  # by 5 km, then mosaic the tiles. Keep the crop at the requested extent.
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
    stac_url = cfg$endpoint,
    collection_id = cfg$collection,
    bbox_vec = bbox_vec,
    start_date = start_date,
    end_date = end_date,
    max_cloud_cover = max_cloud_cover,
    has_cloud = !is.null(cfg$cloud_key))
  if( is.null(items) || length(items[["features"]]) == 0L ){
    return(NULL)
  }

  # The server filters cloud cover. Filter platforms on the returned items.
  features <- items[["features"]]
  if( !is.null(cfg$platforms) ){
    features <- Filter(function(f){
      plat <- f[["properties"]][["platform"]]
      return(!is.null(plat) && plat %in% cfg$platforms)
    }, features)
  }
  if( length(features) == 0L ){
    return(NULL)
  }

  if( identical(res$source, "sentinel-2") ){
    features <- .selectS2Processing(features, s2_processing)
  }

  aoi_sfc <- sf::st_as_sfc(bbox)

  read_band <- function(asset_href, scale, shift = 0){

    if( is.null(asset_href) ){
      return(NULL)
    }
    b <- terra::rast(paste0("/vsicurl/", asset_href))
    terra::scoff(b) <- NULL
    # Read only the requested part of each tile. A STAC footprint can intersect
    # the query even when its raster does not, so check the raster overlap first.
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
    if( scale ){
      if( shift != 0 ){
        b <- b - shift
      }
      if( !is.null(cfg$valid_range) ){
        b <- terra::ifel(b < cfg$valid_range[1] | b > cfg$valid_range[2], NA, b)
      }
      if( !is.null(cfg$scale_factor) ){
        b <- b * cfg$scale_factor + (cfg$scale_offset %||% 0)
      }
    }
    return(b)
  }

  # Mosaic bands for each date before computing the index. Keep the bands
  # so another index can reuse them.
  if( is_sinusoidal ){
    # datetime gives the product date. start_datetime begins its 16-day source
    # window, eight days earlier.
    feat_dates <- vapply(features, function(f){
      ds <- f[["properties"]][["datetime"]]
      if( is.null(ds) || !nzchar(ds) ){
        ds <- f[["properties"]][["start_datetime"]] %||% ""
      }
      return(substr(ds, 1, 10))
    }, character(1L))
    # Window overlap can return products dated outside the requested period.
    in_period <- feat_dates >= format(as.Date(start_date)) &
      feat_dates <= format(as.Date(end_date))
    features <- features[in_period]
    feat_dates <- feat_dates[in_period]
    if( length(features) == 0L ){
      return(NULL)
    }

    date_results <- lapply(unique(feat_dates), function(date_str){
      group <- features[feat_dates == date_str]
      dt <- tryCatch(as.Date(date_str), error = function(e){
        return(NULL)
      })
      if( is.null(dt) || is.na(dt) ){
        return(list(scene = NULL, err = list()))
      }

      tile_reads <- lapply(group, function(f){
        return(tryCatch({
          tb <- lapply(idx$assets, function(an){
            return(read_band(f[["assets"]][[an]][["href"]], scale = TRUE))
          })
          names(tb) <- idx$assets
          list(bands = Filter(Negate(is.null), tb), err = NULL)
        }, error = function(e){
          return(list(bands = NULL, err = e))
        }))
      })
      tile_errors <- Filter(Negate(is.null), lapply(tile_reads, `[[`, "err"))
      per_tile_bands <- lapply(tile_reads, `[[`, "bands")
      per_tile_bands <- Filter(function(x){
        return(!is.null(x) && length(x) > 0L)
      }, per_tile_bands)
      if( length(per_tile_bands) == 0L ){
        return(list(scene = NULL, err = tile_errors))
      }

      bands_merged <- lapply(idx$assets, function(an){
        asset_tiles <- lapply(per_tile_bands, function(tb){
          return(tb[[an]])
        })
        asset_tiles <- Filter(Negate(is.null), asset_tiles)
        if( length(asset_tiles) == 0L ){
          return(NULL)
        }
        if( length(asset_tiles) == 1L ){
          return(asset_tiles[[1L]])
        }
        return(tryCatch(.mergeRasterTiles(asset_tiles, method = "bilinear"),
                        error = function(e){
                          warning(sprintf(
                            "Sinusoidal asset merge failed for asset %s on date %s. Using the first tile only. Reason: %s",
                            an, date_str, .stripUrlQuery(conditionMessage(e))),
                            call. = FALSE)
                          return(asset_tiles[[1L]])
                        }))
      })
      names(bands_merged) <- idx$assets
      if( any(vapply(bands_merged, is.null, logical(1L))) ){
        return(list(scene = NULL, err = tile_errors))
      }

      idx_rast <- tryCatch({
        stk <- do.call(c, unname(bands_merged))
        names(stk) <- names(bands_merged)
        .finiteIndex(idx$fun(stk))
      }, error = function(e){
        return(e)
      })
      if( inherits(idx_rast, "error") ){
        return(list(scene = NULL, err = c(tile_errors, list(idx_rast))))
      }

      # Keep asset URLs for later reads without another search.
      item_assets <- lapply(group[[1L]][["assets"]], function(a){
        return(a[["href"]])
      })

      return(list(scene = list(date = dt, bands = bands_merged, index = idx_rast,
                               item_assets = item_assets),
                  err = tile_errors))
    })
    tile_errors <- do.call(c, lapply(date_results, `[[`, "err"))
    if( length(tile_errors) > 0L ){
      warning(sprintf("MCD43A4 processing errors: %d. First error: %s",
                      length(tile_errors),
                      .stripUrlQuery(conditionMessage(tile_errors[[1L]]))),
              call. = FALSE)
    }
    scenes <- lapply(date_results, `[[`, "scene")
    scenes <- Filter(Negate(is.null), scenes)
    if( length(scenes) == 0L ){
      return(NULL)
    }
    scenes <- scenes[order(vapply(scenes, function(s){
      return(as.numeric(s$date))
    }, numeric(1)))]
    return(scenes)
  }

  # Apply offsets per item: processing baselines can differ on the same date.
  # Only reflectance assets use the offset.
  offsets <- if( identical(res$source, "sentinel-2") ){
    vapply(features, .s2Offset, numeric(1L))
  }else{
    rep(0, length(features))
  }
  offset_assets <- unique(unlist(lapply(s2_index_list, `[[`, "assets"),
                                 use.names = FALSE))

  # Track read errors separately from tiles outside the query area.
  read_results <- mapply(function(f, offset){
    shift <- if( is.na(offset) ) { 0 } else { offset }
    # Treat an empty datetime like NULL before falling back to start_datetime.
    dt_str <- f[["properties"]][["datetime"]]
    if( is.null(dt_str) || !nzchar(dt_str) ){
      dt_str <- f[["properties"]][["start_datetime"]]
    }
    dt <- tryCatch(as.Date(substr(dt_str, 1, 10)), error = function(e){
      return(NULL)
    })
    if( is.null(dt) || is.na(dt) ){
      return(list(scene = NULL, err = NULL))
    }

    out <- tryCatch({
      pixel_mask <- NULL
      bands <- lapply(idx$assets, function(an){
        return(read_band(f[["assets"]][[an]][["href"]], scale = TRUE, shift = shift))
      })
      names(bands) <- idx$assets
      bands <- Filter(Negate(is.null), bands)
      # This return leaves the item reader, including the enclosing tryCatch().
      if( length(bands) == 0L ){
        return(list(scene = NULL, err = NULL))
      }
      bands <- .harmoniseBands(bands)
      stk <- do.call(c, unname(bands))
      names(stk) <- names(bands)

      if( !is.null(cfg$mask_asset) && !is.null(scl_classes) ){
        mb <- read_band(f[["assets"]][[cfg$mask_asset]][["href"]], scale = FALSE)
        if( !is.null(mb) ){
          mb <- terra::resample(mb, stk, method = "near")
          cloud_mask <- if( identical(cfg$mask_kind, "scl") ){
            Reduce(`|`, lapply(scl_classes, function(cls){
              return(mb == cls)
            }))
          }else{
            bit_sum <- sum(vapply(scl_classes,
                                  function(b){
                                    return(bitwShiftL(1L, as.integer(b)))
                                  },
                                  integer(1)))
            terra::app(mb, function(v){
              return(bitwAnd(as.integer(v), bit_sum) != 0L)
            })
          }
          pixel_mask <- cloud_mask
          stk <- terra::mask(stk, cloud_mask, maskvalue = TRUE)
          # Mask bands too, since another index may reuse them.
          bands <- lapply(bands, function(b){
            return(terra::mask(b, cloud_mask, maskvalue = TRUE))
          })
        }
      }
      idx_rast <- .finiteIndex(idx$fun(stk))
      item_assets <- lapply(f[["assets"]], function(a){
        return(a[["href"]])
      })
      # Retain offsets for later asset reads.
      dn_offset <- NULL
      if( shift != 0 ){
        shifted <- intersect(names(item_assets), offset_assets)
        dn_offset <- lapply(item_assets[shifted], function(h){
          return(rep(shift, length(h)))
        })
      }
      list(bands = bands, index = idx_rast, item_assets = item_assets,
           pixel_mask = pixel_mask, dn_offset = dn_offset)
    }, error = function(e){
      return(e)
    })

    if( inherits(out, "error") ){
      return(list(scene = NULL, err = out))
    }
    scene <- list(date = dt,
                  bands = out$bands,
                  index = out$index,
                  item_assets = out$item_assets,
                  pixel_mask = out$pixel_mask)
    if( !is.null(out$dn_offset) ){
      scene$dn_offset <- out$dn_offset
    }
    return(list(scene = scene, err = NULL))
  }, features, offsets, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  read_errors <- Filter(Negate(is.null), lapply(read_results, `[[`, "err"))
  if( length(read_errors) > 0L ){
    warning(sprintf("%d of %d acquisitions dropped after an error. First error: %s",
                    length(read_errors), length(features),
                    .stripUrlQuery(conditionMessage(read_errors[[1L]]))),
            call. = FALSE)
  }
  kept <- !vapply(read_results, function(r){
    return(is.null(r$scene))
  }, logical(1L))
  n_unread <- sum(is.na(offsets) & kept)
  if( n_unread > 0L ){
    warning(sprintf(paste0("No offset was removed from %d Sentinel-2 item(s), ",
                           "because the processing baseline could not be read."),
                    n_unread),
            call. = FALSE)
  }
  scenes <- lapply(read_results, `[[`, "scene")
  scenes <- Filter(Negate(is.null), scenes)
  if( length(scenes) == 0L ){
    return(NULL)
  }

  # Combine adjacent tiles before aligning dates, preserving coverage and
  # one scene per date.
  scenes <- .mosaicFineScenesByDate(scenes, idx)
  if( length(scenes) == 0L ){
    return(NULL)
  }

  scenes <- scenes[order(vapply(scenes, function(s){
    return(as.numeric(s$date))
  }, numeric(1)))]

  # The index grid must cover every date's extent. Bands keep their source grids.
  if( length(scenes) > 1L ){
    template <- terra::rast(.mergeRasterTiles(
      lapply(scenes, `[[`, "index"), method = "bilinear"))
    for( i in seq_along(scenes) ){
      if( !terra::compareGeom(scenes[[i]]$index, template,
                              stopOnError = FALSE) ){
        scenes[[i]]$index <- terra::project(scenes[[i]]$index, template,
                                            method = "bilinear")
      }
    }
  }
  return(scenes)
}
