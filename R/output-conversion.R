# =========================================================================
# Output conversion and mask construction (matrix-based)
# =========================================================================

#' Build per-area + dissolved sf masks from alive_mat (internal)
#'
#' Uses areas_info (which carries the strict-separation result if applied).
#'
#' @noRd
.buildMaskSf <- function(alive_mat, areas_info, active_rast,
                         multiple_areas){

  if( !any(alive_mat) || areas_info$n_areas == 0L )
    return(list(mask_sf_per_area = NULL,
                mask_sf_dissolved = NULL,
                n_areas          = 0L))

  ## Build a SpatRaster carrying area rank per cell.
  patch <- terra::rast(active_rast[[1L]])
  terra::values(patch) <- NA_real_

  if( isFALSE(multiple_areas) ){
    keep <- areas_info$label_mat == areas_info$area_order[1L]
    keep_cells <- which(keep, arr.ind = TRUE)
    if( nrow(keep_cells) == 0L )
      return(list(mask_sf_per_area = NULL,
                  mask_sf_dissolved = NULL, n_areas = 0L))
    cells <- terra::cellFromRowCol(patch, keep_cells[, "row"],
                                   keep_cells[, "col"])
    patch[cells] <- 1
    n_areas <- 1L
  }else{
    rank_mat <- matrix(0L, nrow(alive_mat), ncol(alive_mat))
    for( ai in seq_len(areas_info$n_areas) ){
      lab <- ai
      rnk <- areas_info$areas[[ai]]$rank
      rank_mat[areas_info$label_mat == lab] <- rnk
    }
    keep_cells <- which(rank_mat > 0L, arr.ind = TRUE)
    if( nrow(keep_cells) == 0L )
      return(list(mask_sf_per_area = NULL,
                  mask_sf_dissolved = NULL, n_areas = 0L))
    cells <- terra::cellFromRowCol(patch, keep_cells[, "row"],
                                   keep_cells[, "col"])
    patch[cells] <- rank_mat[keep_cells]
    n_areas <- areas_info$n_areas
  }
  names(patch) <- "area_num"

  trimmed <- terra::trim(patch)
  poly_utm <- terra::as.polygons(trimmed, dissolve = TRUE)
  poly_wgs <- terra::project(poly_utm, "EPSG:4326")
  per_area <- sf::st_as_sf(poly_wgs)
  attr_cols <- setdiff(names(per_area), attr(per_area, "sf_column"))
  if( length(attr_cols) >= 1L ){
    names(per_area)[names(per_area) == attr_cols[1L]] <- "area_num"
    per_area$area_num <- as.integer(per_area$area_num)
    per_area <- per_area[order(per_area$area_num), , drop = FALSE]
    rownames(per_area) <- NULL
  }
  dissolved <- sf::st_sf(geometry = sf::st_union(per_area))

  return(list(mask_sf_per_area = per_area,
              mask_sf_dissolved = dissolved,
              n_areas          = as.integer(n_areas)))
}

# =========================================================================
# Composite scenes (internal)
# =========================================================================

#' Validate an optional user-supplied output mask (internal)
#' @noRd
.validateOutputMask <- function(mask, caller){
  if( is.null(mask) ) return(NULL)
  if( !inherits(mask, c("sf", "sfc")) )
    stop(caller, ": mask must be NULL or an sf/sfc geometry.",
         call. = FALSE)
  out <- if( inherits(mask, "sfc") ) sf::st_sf(geometry = mask) else mask
  if( nrow(out) == 0L || all(sf::st_is_empty(out)) )
    stop(caller, ": mask must contain at least one non-empty geometry.",
         call. = FALSE)
  if( is.na(sf::st_crs(out)) )
    stop(caller, ": mask must have a coordinate reference system.",
         call. = FALSE)
  out
}

#' Composite scenes per pixel across time (internal)
#' @noRd
.compositeScenes <- function(scenes, composite_function = "median"){

  if( is.null(scenes) || length(scenes) == 0L ) return(NULL)
  composite_function <- match.arg(composite_function,
                                  c("median", "mean", "max", "min", "sum"))
  stk <- Reduce(c, lapply(scenes, `[[`, "index"))
  out <- switch(composite_function,
                "median" = terra::median(stk, na.rm = TRUE),
                "mean"   = terra::mean(stk,   na.rm = TRUE),
                "max"    = max(stk, na.rm = TRUE),
                "min"    = min(stk, na.rm = TRUE),
                "sum"    = sum(stk, na.rm = TRUE))
  ## In particular, sum(..., na.rm = TRUE) maps an all-missing cell to zero.
  ## Preserve missingness for every composite function.
  valid_n <- terra::app(stk, function(v) sum(is.finite(v)))
  return(terra::ifel(valid_n == 0L, NA, out))
}

# =========================================================================
# rs_meta attachment (used by validate_data_fusion)
# =========================================================================

#' Attach rs metadata to an output (internal)
#' @noRd
.attachRsMeta <- function(out, rs){

  attr(out, "rs_meta") <- list(
    longitude  = rs$spec$point$longitude,
    latitude   = rs$spec$point$latitude,
    start_date = rs$spec$start_date,
    end_date   = rs$spec$end_date,
    index_name = rs$spec$index_name,
    source     = rs$geom$source,
    max_cloud_cover = rs$spec$max_cloud_cover,
    scl_classes     = rs$spec$scl_classes,
    params          = rs$params,
    x_metres   = rs$spec$x_metres,
    y_metres   = rs$spec$y_metres)
  return(out)
}

# =========================================================================
# as_raster
# =========================================================================

#' Create a composite spectral index raster
#'
#' Combines the requested index values across retrieved scenes for each pixel.
#' Before field boundary delineation, the result covers the full analysis area;
#' after delineation, it is masked to the delineated area or areas. Supply
#' \code{mask} to use an sf/sfc geometry directly.
#'
#' Parameters stored in \code{rs$params} are used by default. Named arguments
#' supplied here override the corresponding setting for this call only.
#'
#' @param x An \code{ocular} object, or a longitude in standalone use.
#' @param composite_function Summary used to combine each pixel across scenes:
#'   \code{"median"} (the default), \code{"mean"}, \code{"max"},
#'   \code{"min"}, or \code{"sum"}.
#' @param params Optional complete or partial list of settings accepted by
#'   \code{rs_params()}. These settings override \code{rs$params} for this
#'   call; \code{NULL} retains the stored settings.
#' @param multiple_areas Override for the corresponding \code{rs_params()}
#'   setting, controlling whether separate delineated areas are returned.
#' @param mask Optional sf/sfc geometry used directly as the area mask,
#'   without first requiring field boundary delineation. It may use any
#'   coordinate reference system and is reprojected to the composite raster.
#'   \code{NULL} uses the analysis area or delineated area, as applicable.
#' @param ... In standalone use, additional arguments passed to \code{get_rs()}.
#'   With an \code{ocular} object, further named \code{rs_params()} overrides
#'   for this call.
#' @returns A \code{SpatRaster}. With multiple delineated areas, layers are named
#'   \code{area_1}, \code{area_2}, and so on, followed by \code{union}. If
#'   delineation ran but found no area, values are all \code{NA}; an
#'   undelineated object returns a composite for the full analysis area.
#' @seealso \code{\link{as_time_series}}, \code{\link{get_rs}}
#' @export
as_raster <- function(x, composite_function = "median",
                      params = NULL,
                      multiple_areas = NULL,
                      mask = NULL,
                      ...){

  composite_function <- match.arg(composite_function,
                                  c("median", "mean", "max", "min", "sum"))
  mask <- .validateOutputMask(mask, "as_raster")

  formal_overrides <- list(
    multiple_areas = multiple_areas
  )
  if( is_rs(x) ){
    rs <- x
    params <- .resolveParams(rs, params, c(formal_overrides, list(...)))
  }else{
    rs <- get_rs(longitude = x, params = params, ...)
    params <- .resolveParams(rs, NULL, formal_overrides)
  }

  if( is.null(rs$scenes) || length(rs$scenes) == 0L )
    stop("as_raster: rs has no scenes.", call. = FALSE)

  composite <- .compositeScenes(rs$scenes, composite_function)
  ## Crop composite to the user's analysed area in target UTM, then to
  ## scenes' CRS. Anchors output to the source-scene grid via snap="out"
  ## -- independent of the (possibly resampled) detection feature_stack.
  fs_poly <- terra::as.polygons(.analysedExt(rs), crs = rs_utm(rs))
  fs_in_crs <- terra::project(fs_poly, terra::crs(composite))
  composite <- terra::crop(composite, terra::ext(fs_in_crs), snap = "out")

  ## User-supplied mask wins; bypass delineation-derived geometry.
  if( !is.null(mask) ){
    mv <- terra::project(terra::vect(mask), terra::crs(composite))
    return(.attachRsMeta(terra::mask(composite, mv), rs))
  }

  ## NULL means delineation has not run; an all-FALSE matrix means it ran and
  ## found no area. Do not turn the latter into a plausible whole-AOI
  ## output.
  if( is.null(rs$state$alive_mat) ){
    return(.attachRsMeta(composite, rs))
  }
  if( !any(rs$state$alive_mat) )
    return(.attachRsMeta(composite * NA_real_, rs))

  areas_info <- .applyStrictSeparation(rs$state$alive_mat, rs$geom$centre_rc,
                                       params$area_separation_strict)
  mp <- .buildMaskSf(rs$state$alive_mat, areas_info, rs$internals$feature_stack,
                     params$multiple_areas)
  if( is.null(mp$mask_sf_dissolved) )
    return(.attachRsMeta(composite, rs))

  ## Single-area
  if( isFALSE(params$multiple_areas) || mp$n_areas <= 1L ){
    mv <- terra::project(terra::vect(mp$mask_sf_dissolved), terra::crs(composite))
    return(.attachRsMeta(terra::mask(composite, mv), rs))
  }
  ## Multi-area: per-area layers + union
  per <- mp$mask_sf_per_area
  layers <- lapply(seq_len(nrow(per)), function(i){
    mv <- terra::project(terra::vect(per[i, ]), terra::crs(composite))
    layer <- terra::mask(composite, mv)
    names(layer) <- paste0("area_", per$area_num[i])
    layer
  })
  union_layer <- terra::mask(composite,
                             terra::project(terra::vect(mp$mask_sf_dissolved),
                                            terra::crs(composite)))
  names(union_layer) <- "union"
  out <- Reduce(c, c(layers, list(union_layer)))
  return(.attachRsMeta(out, rs))
}

# =========================================================================
# as_time_series
# =========================================================================

#' Per-mask scene aggregation (internal)
#'
#' Projects each mask once to scene CRS (caches by CRS string).
#'
#' @noRd
.aggregateMasks <- function(scenes, bbox, masks_named,
                            time_aggregate, aggregate_function){

  if( is.null(scenes) || length(scenes) == 0L ){
    return(data.frame(date = character(0), value = numeric(0),
                      stringsAsFactors = FALSE))
  }
  mask_names <- names(masks_named)
  has_fusion_meta <- !is.null(scenes[[1L]]$fused)

  ## Group scenes by their CRS so masks are projected once per group.
  scene_crs <- vapply(scenes, function(sc) terra::crs(sc$index), character(1))
  uniq_crs <- unique(scene_crs)
  proj_cache <- list()
  for( cs in uniq_crs ){
    proj_cache[[cs]] <- lapply(masks_named, function(m){
      if( is.null(m) || nrow(m) == 0L ) NULL
      else terra::project(terra::vect(m), cs)
    })
  }

  plot_sfc <- sf::st_as_sfc(bbox)
  day   <- vapply(scenes, function(sc) as.character(sc$date), character(1))
  dates <- switch(time_aggregate,
                  "monthly" = substr(day, 1, 7),
                  "yearly"  = substr(day, 1, 4),
                  day)

  scene_stats <- lapply(seq_along(scenes), function(si){
    sc <- scenes[[si]]
    r           <- sc$index
    plot_in_crs <- sf::st_transform(plot_sfc, terra::crs(r))
    r_crop      <- terra::crop(r, terra::ext(sf::st_bbox(plot_in_crs)),
                               snap = "out")
    if( terra::ncell(r_crop) == 0L ){
      vals_named <- stats::setNames(rep(NA_real_, length(mask_names)), mask_names)
      return(list(values = vals_named, mae = NA_real_, nmae = NA_real_))
    }
    m_projs <- proj_cache[[scene_crs[si]]]
    vals_named <- vapply(seq_along(mask_names), function(i){
      r_m <- if( !is.null(m_projs[[i]]) ) terra::mask(r_crop, m_projs[[i]])
      else r_crop
      n_valid <- terra::global(is.finite(r_m), "sum", na.rm = TRUE)[[1L]]
      v <- if( !is.finite(n_valid) || n_valid == 0L ) NA_real_ else
        terra::global(r_m, aggregate_function, na.rm = TRUE)[[1L]]
      as.numeric(v)
    }, numeric(1))
    names(vals_named) <- mask_names

    if( isTRUE(sc$fused) && !is.null(sc$residual) ){
      diagnostic_i <- match("union", mask_names)
      if( is.na(diagnostic_i) ) diagnostic_i <- 1L
      diagnostic_mask <- m_projs[[diagnostic_i]]
      resid <- terra::crop(sc$residual,
                           terra::ext(sf::st_bbox(plot_in_crs)), snap = "out")
      lst1  <- terra::crop(sc$landsat_t1,
                           terra::ext(sf::st_bbox(plot_in_crs)), snap = "out")
      if( !is.null(diagnostic_mask) ){
        resid <- terra::mask(resid, diagnostic_mask)
        lst1  <- terra::mask(lst1,  diagnostic_mask)
      }
      mae     <- as.numeric(terra::global(resid, "mean", na.rm = TRUE)[[1L]])
      mean_ls <- as.numeric(terra::global(lst1,  "mean", na.rm = TRUE)[[1L]])
      nmae    <- if( is.finite(mean_ls) && abs(mean_ls) > 1e-6 )
        mae / abs(mean_ls) else NA_real_
      return(list(values = vals_named, mae = mae, nmae = nmae))
    }
    list(values = vals_named, mae = NA_real_, nmae = NA_real_)
  })

  values_matrix <- do.call(rbind, lapply(scene_stats, `[[`, "values"))
  values_df     <- as.data.frame(values_matrix, stringsAsFactors = FALSE)
  names(values_df) <- mask_names

  df <- if( has_fusion_meta ){
    cbind(
      data.frame(date = dates, day = day, stringsAsFactors = FALSE),
      values_df,
      data.frame(
        is_fused    = vapply(scenes, `[[`, logical(1), "fused"),
        anchor_days = vapply(scenes, function(sc) sc$anchor_days %||% NA_integer_, integer(1)),
        anchor_mae  = vapply(scene_stats, `[[`, numeric(1), "mae"),
        anchor_nmae = vapply(scene_stats, `[[`, numeric(1), "nmae"),
        stringsAsFactors = FALSE))
  }else{
    cbind(data.frame(date = dates, stringsAsFactors = FALSE), values_df)
  }

  if( length(mask_names) == 1L ){
    keep <- is.finite(df[[mask_names]])
  }else{
    val_mat <- as.matrix(df[, mask_names, drop = FALSE])
    keep    <- rowSums(!is.finite(val_mat)) < ncol(val_mat)
  }
  df <- df[keep, , drop = FALSE]
  if( nrow(df) == 0L ){
    df$day <- NULL  ## helper column (fusion path only); keep schema uniform
    return(df)
  }

  if( has_fusion_meta ){
    ## Remove a fused estimate that duplicates a real same-day acquisition
    ## before temporal aggregation. Ordering by day and is_fused
    ## keeps the real (is_fused=FALSE) scene over the fused one. Deduping on the
    ## bucketed date would collapse each month/year to a single scene.
    df <- df[order(df$day, df$is_fused), , drop = FALSE]
    df <- df[!duplicated(df$day), , drop = FALSE]
    if( !identical(time_aggregate, "daily") ){
      agg_cols <- c(mask_names, "anchor_days", "anchor_mae", "anchor_nmae")
      ## NB: no na.action arg -- aggregate.data.frame has no such parameter
      ## (only aggregate.formula does); it would fall through `...` into
      ## The aggregate function handles missing values through na.rm = TRUE.
      df <- stats::aggregate(
        df[, agg_cols, drop = FALSE],
        by  = list(date = df$date, is_fused = df$is_fused),
        FUN = function(x){
          x <- x[is.finite(x)]
          if( length(x) == 0L ) NA_real_
          else do.call(aggregate_function, list(x, na.rm = TRUE))
        })
    }else{
      df$day <- NULL
    }
  }else{
    df <- stats::aggregate(df[, mask_names, drop = FALSE],
                           by  = list(date = df$date),
                           FUN = function(x){
                             x <- x[is.finite(x)]
                             if( length(x) == 0L ) NA_real_
                             else do.call(aggregate_function,
                                          list(x, na.rm = TRUE))
                           })
  }
  return(df)
}

#' Create a spectral index time series
#'
#' Summarises the requested index values within the analysis area or delineated
#' field for each retrieved scene. Results may be grouped by day, month, or
#' year. Supply \code{mask} to use an sf/sfc geometry directly.
#'
#' Parameters stored in \code{rs$params} are used by default. Named arguments
#' supplied here override the corresponding setting for this call only.
#'
#' @param x An \code{ocular} object, or a longitude in standalone use.
#' @param time_aggregate Temporal grouping: \code{"daily"} (the default),
#'   \code{"monthly"}, or \code{"yearly"}.
#' @param aggregate_function A function, or its name, used to summarise values
#'   within each area and time group. It must accept numeric input and
#'   \code{na.rm}, and return one numeric value. The default is \code{"mean"}.
#' @param params Optional complete or partial list of settings accepted by
#'   \code{rs_params()}. These settings override \code{rs$params} for this
#'   call; \code{NULL} retains the stored settings.
#' @param multiple_areas Override for the corresponding \code{rs_params()}
#'   setting, controlling whether separate delineated areas are returned.
#' @param mask Optional sf/sfc geometry used directly as the area mask,
#'   without first requiring field boundary delineation. It may use any
#'   coordinate reference system and is reprojected for each scene. \code{NULL}
#'   uses the analysis area or delineated area, as applicable.
#' @param ... In standalone use, additional arguments passed to \code{get_rs()}.
#'   With an \code{ocular} object, further named \code{rs_params()} overrides
#'   for this call.
#' @returns A \code{data.frame}. Single-area or undelineated output contains
#'   \code{date} and \code{value}. Multiple-area output contains \code{date},
#'   one \code{area_*} column per component, and \code{union}. A series with
#'   fused estimates also contains \code{is_fused}, \code{anchor_days},
#'   \code{anchor_mae}, and \code{anchor_nmae}. The latter two measure mismatch
#'   between the fine and coarse anchors over the diagnostic mask (the union in
#'   multiple-area output); they are not held-out prediction errors.
#'   \code{is_fused} does not denote a raw MODIS observation. An empty
#'   delineation returns zero rows.
#' @seealso \code{\link{as_raster}}, \code{\link{get_rs}}
#' @export
as_time_series <- function(x, time_aggregate = "daily",
                           aggregate_function = "mean",
                           params = NULL,
                           multiple_areas = NULL,
                           mask = NULL,
                           ...){

  valid_time_aggregates <- c("daily", "monthly", "yearly")
  if( !is.character(time_aggregate) || length(time_aggregate) != 1L ||
      is.na(time_aggregate) || !nzchar(time_aggregate) ||
      !(time_aggregate %in% valid_time_aggregates) )
    stop("as_time_series: time_aggregate must be daily, monthly, or yearly.",
         call. = FALSE)
  if( !(is.function(aggregate_function) ||
        (is.character(aggregate_function) &&
         length(aggregate_function) == 1L && !is.na(aggregate_function) &&
         nzchar(aggregate_function))) )
    stop("as_time_series: aggregate_function must be a function or its name.",
         call. = FALSE)
  aggregate_function <- tryCatch(
    match.fun(aggregate_function),
    error = function(e)
      stop("as_time_series: unknown aggregate_function: ",
           conditionMessage(e), call. = FALSE))
  probe <- tryCatch(do.call(aggregate_function,
                            list(c(1, 2, NA_real_), na.rm = TRUE)),
                    error = identity)
  if( inherits(probe, "error") || !is.numeric(probe) || length(probe) != 1L )
    stop("as_time_series: aggregate_function must accept numeric input plus ",
         "na.rm and return one numeric value.", call. = FALSE)
  mask <- .validateOutputMask(mask, "as_time_series")

  formal_overrides <- list(
    multiple_areas = multiple_areas
  )
  if( is_rs(x) ){
    rs <- x
    params <- .resolveParams(rs, params, c(formal_overrides, list(...)))
  }else{
    rs <- get_rs(longitude = x, params = params, ...)
    params <- .resolveParams(rs, NULL, formal_overrides)
  }
  if( is.null(rs$scenes) || length(rs$scenes) == 0L ){
    return(.attachRsMeta(
      data.frame(date = character(0), value = numeric(0),
                 stringsAsFactors = FALSE), rs))
  }

  ## Delineation ran but found no field. Returning the whole analysed box here
  ## would silently convert a failed detection into plausible observations.
  if( is.null(mask) && !is.null(rs$state$alive_mat) &&
      !any(rs$state$alive_mat) )
    return(.attachRsMeta(
      data.frame(date = character(0), value = numeric(0),
                 stringsAsFactors = FALSE), rs))

  masks_named <- list(value = NULL)
  if( !is.null(mask) ){
    masks_named <- list(value = mask)
  }else if( !is.null(rs$state$alive_mat) && any(rs$state$alive_mat) ){
    areas_info <- .applyStrictSeparation(rs$state$alive_mat, rs$geom$centre_rc,
                                         params$area_separation_strict)
    mp <- .buildMaskSf(rs$state$alive_mat, areas_info, rs$internals$feature_stack,
                       params$multiple_areas)
    if( !is.null(mp$mask_sf_dissolved) ){
      if( isFALSE(params$multiple_areas) || mp$n_areas <= 1L ){
        masks_named <- list(value = mp$mask_sf_dissolved)
      }else{
        per <- mp$mask_sf_per_area
        parts <- lapply(seq_len(nrow(per)), function(i) per[i, ])
        names(parts) <- paste0("area_", per$area_num)
        masks_named <- c(parts, list(union = mp$mask_sf_dissolved))
      }
    }
  }

  ## Centred-area WGS bbox: derive from rs's analysed extent so the
  ## per-scene crop is anchored to the user's requested area, not to
  ## the (possibly resampled) detection feature_stack.
  fs_poly <- terra::as.polygons(.analysedExt(rs), crs = rs_utm(rs))
  fs_wgs  <- terra::project(fs_poly, "EPSG:4326")
  out_bbox <- sf::st_bbox(sf::st_as_sf(fs_wgs))

  df <- .aggregateMasks(scenes             = rs$scenes,
                        bbox               = out_bbox,
                        masks_named        = masks_named,
                        time_aggregate     = time_aggregate,
                        aggregate_function = aggregate_function)
  return(.attachRsMeta(df, rs))
}
