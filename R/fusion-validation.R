# =========================================================================
# Fusion agreement calculations
# =========================================================================

#' Estimate Landsat scenes by leaving each one out in turn
#'
#' Use the nearest remaining Landsat scene as the anchor and retrieve MCD43A4
#' near both dates. Return the estimates and observations for later comparison.
#' Return \code{NULL} if fewer than two Landsat scenes are available, the index
#' cannot be fused, or no omitted scene can be estimated.
#' @noRd
.validateCompute <- function(rs, params){

  if( !rs$geom$source %in% c("landsat-8", "landsat-5") ){
    stop("validate_data_fusion: Landsat scenes are required. The input source ",
         "is '", rs$geom$source, "'. Re-fetch with source = 'landsat-8' ",
         "or source = 'landsat-5'.", call. = FALSE)
  }
  ls_scenes <- rs$scenes
  if( length(ls_scenes) < 2L ){
    return(NULL)
  }

  if( rs$spec$index_name %in% .fusionUnsupportedIndices() ){
    warning("validate_data_fusion: MODIS fusion not supported for index '",
            rs$spec$index_name, "'. Validation skipped.", call. = FALSE)
    return(NULL)
  }

  l_entry <- landsat_index_list[[rs$spec$index_name]]
  m_entry <- mcd43a4_index_list[[rs$spec$index_name]]
  if( is.null(l_entry) || is.null(m_entry) ){
    return(NULL)
  }
  m_assets <- m_entry$assets

  # Validation needs coarse scenes near the Landsat dates, not every day
  # in the requested period.
  ls_date_vec <- vapply(ls_scenes, function(s){
    return(as.numeric(s$date))
  }, numeric(1))
  coarse_error <- NULL
  coarse_scenes <- tryCatch(
    .fetchCoarseNearDates(bbox = rs$geom$bbox,
                          fine_dates = ls_date_vec,
                          index_name = rs$spec$index_name,
                          max_cloud_cover = rs$spec$max_cloud_cover,
                          source = "mcd43a4"),
    error = function(e){
      coarse_error <<- conditionMessage(e)
      return(NULL)
    })
  if( !is.null(coarse_error) ){
    warning("validate_data_fusion: MCD43A4 retrieval failed: ", coarse_error,
            call. = FALSE)
  }
  if( is.null(coarse_scenes) || length(coarse_scenes) == 0L ){
    return(NULL)
  }

  ref_r <- ls_scenes[[1L]]$index
  coarse_scenes <- lapply(coarse_scenes, function(sc){
    b <- .orderedBands(sc$bands, m_assets)
    if( is.null(b) ){
      return(NULL)
    }
    return(list(date = sc$date, bands = .bandsToMemory(b)))
  })
  coarse_scenes <- Filter(Negate(is.null), coarse_scenes)
  if( length(coarse_scenes) == 0L ){
    return(NULL)
  }

  coarse_dates <- vapply(coarse_scenes, function(sc){
    return(as.numeric(sc$date))
  }, numeric(1))

  coarse_tol_d <- 1L
  fusion_errors <- character(0L)
  n_skipped <- 0L
  results <- lapply(seq_along(ls_scenes), function(holdout_i){
    target <- ls_scenes[[holdout_i]]
    remaining <- ls_scenes[-holdout_i]
    if( length(remaining) == 0L ){
      return(NULL)
    }
    rem_dates <- vapply(remaining, function(sc){
      return(as.numeric(sc$date))
    }, numeric(1))
    anchor <- remaining[[which.min(abs(rem_dates -
                                            as.numeric(target$date)))]]
    fine_raw <- .orderedBands(anchor$bands, l_entry$assets)
    if( is.null(fine_raw) ){
      return(NULL)
    }
    fine_bands <- lapply(fine_raw, function(b){
      return(terra::project(b, ref_r, method = "bilinear"))
    })
    names(fine_bands) <- l_entry$assets
    i_t1 <- which.min(abs(coarse_dates - as.numeric(anchor$date)))
    i_t2 <- which.min(abs(coarse_dates - as.numeric(target$date)))
    # Both dates need a nearby coarse scene. Reusing one would produce no change.
    if( abs(coarse_dates[i_t1] - as.numeric(anchor$date)) > coarse_tol_d ||
        abs(coarse_dates[i_t2] - as.numeric(target$date)) > coarse_tol_d ||
        i_t1 == i_t2 ){
      n_skipped <<- n_skipped + 1L
      return(NULL)
    }
    c_t1 <- coarse_scenes[[i_t1]]$bands
    c_t2 <- coarse_scenes[[i_t2]]$bands
    fused <- tryCatch(
      .fuseLandsatMODIS(fine_bands, c_t1, c_t2, l_entry$fun),
      error = function(e){
        fusion_errors <<- c(
          fusion_errors,
          paste0(as.character(target$date), ": ", conditionMessage(e)))
        return(NULL)
      })
    if( is.null(fused) ){
      return(NULL)
    }
    p <- fused$index
    if( !terra::compareGeom(p, target$index, stopOnError = FALSE) ){
      p <- terra::resample(p, target$index, method = "bilinear")
    }
    return(list(date = target$date,
                predicted_rast = p,
                target_rast = target$index))
  })
  if( length(fusion_errors) > 0L ){
    warning("validate_data_fusion: fusion failed for ",
            length(fusion_errors), " held-out scene(s). First error: ",
            fusion_errors[[1L]], call. = FALSE)
  }
  if( n_skipped > 0L ){
    warning("validate_data_fusion: ", n_skipped, " holdout(s) skipped, because ",
            "the anchor and the target date had no separate MCD43A4 scene ",
            "within 1 day.", call. = FALSE)
  }
  results <- Filter(Negate(is.null), results)
  if( length(results) == 0L ){
    return(NULL)
  }
  return(results)
}

#' Measure agreement between paired observations and estimates
#'
#' Return RMSE, MAE, bias (estimate minus observation) and Pearson correlation.
#' Empty inputs give \code{NA_real_}. Correlation also requires two pairs
#' and nonzero variance in both vectors.
#' @noRd
.agreementMetrics <- function(obs, pred){

  n <- length(obs)
  if( n < 1L ){
    return(list(rmse = NA_real_, mae = NA_real_,
                bias = NA_real_, r = NA_real_))
  }
  diff <- pred - obs
  rmse <- sqrt(mean(diff^2))
  mae <- mean(abs(diff))
  bias <- mean(diff)
  r <- if( n >= 2L ){
    sd_o <- stats::sd(obs)
    sd_p <- stats::sd(pred)
    if( !is.finite(sd_o) || !is.finite(sd_p) || sd_o == 0 || sd_p == 0 ){
      NA_real_
    }else{
      stats::cor(obs, pred)
    }
  }else{
    NA_real_
  }
  return(list(rmse = rmse, mae = mae, bias = bias, r = r))
}

#' Summarise fusion agreement for each scene
#'
#' Compare paired pixels after excluding missing values. Return missing metrics
#' when fewer than \code{min_valid_n} pairs remain.
#' @noRd
.collapseToScene <- function(loo_results, min_valid_n = 4L){

  rows <- lapply(loo_results, function(res){
    obs <- terra::values(res$target_rast,    na.rm = FALSE)
    pred <- terra::values(res$predicted_rast, na.rm = FALSE)
    valid <- is.finite(obs) & is.finite(pred)
    n_pix <- sum(valid)
    if( n_pix < min_valid_n ){
      data.frame(date = as.character(res$date),
                 rmse = NA_real_,
                 mae = NA_real_,
                 bias = NA_real_,
                 r = NA_real_,
                 n_pixels = n_pix,
                 stringsAsFactors = FALSE)
    }else{
      m <- .agreementMetrics(obs[valid], pred[valid])
      data.frame(date = as.character(res$date),
                 rmse = m$rmse,
                 mae = m$mae,
                 bias = m$bias,
                 r = m$r,
                 n_pixels = n_pix,
                 stringsAsFactors = FALSE)
    }
  })
  return(do.call(rbind, rows))
}

#' Summarise fusion agreement for each pixel
#'
#' Compare paired values across scenes after excluding missing values. Return
#' missing metrics when fewer than \code{min_valid_n} pairs remain.
#'
#' @returns A raster with \code{rmse}, \code{mae}, \code{bias}, \code{r} and
#'   \code{n_scenes} layers on the template grid.
#' @noRd
.collapseToPixel <- function(loo_results, template, min_valid_n = 4L){

  n_scenes <- length(loo_results)
  nr <- terra::nrow(template)
  nc <- terra::ncol(template)
  n_cells <- nr * nc

  pred_mat <- matrix(NA_real_, n_cells, n_scenes)
  obs_mat <- matrix(NA_real_, n_cells, n_scenes)
  for( i in seq_along(loo_results) ){
    p <- loo_results[[i]]$predicted_rast
    t <- loo_results[[i]]$target_rast
    if( !terra::compareGeom(p, template, stopOnError = FALSE) ){
      p <- terra::resample(p, template, method = "bilinear")
    }
    if( !terra::compareGeom(t, template, stopOnError = FALSE) ){
      t <- terra::resample(t, template, method = "bilinear")
    }
    pred_mat[, i] <- terra::values(p, na.rm = FALSE)
    obs_mat[, i] <- terra::values(t, na.rm = FALSE)
  }

  rmse_v <- numeric(n_cells)
  mae_v <- numeric(n_cells)
  bias_v <- numeric(n_cells)
  r_v <- numeric(n_cells)
  ns_v <- integer(n_cells)
  for( k in seq_len(n_cells) ){
    obs_k <- obs_mat[k, ]
    pred_k <- pred_mat[k, ]
    valid <- is.finite(obs_k) & is.finite(pred_k)
    n_k <- sum(valid)
    ns_v[k] <- n_k
    if( n_k < min_valid_n ){
      rmse_v[k] <- NA_real_
      mae_v[k] <- NA_real_
      bias_v[k] <- NA_real_
      r_v[k] <- NA_real_
    }else{
      m <- .agreementMetrics(obs_k[valid], pred_k[valid])
      rmse_v[k] <- m$rmse
      mae_v[k] <- m$mae
      bias_v[k] <- m$bias
      r_v[k] <- m$r
    }
  }

  rmse_r <- terra::rast(template)
  terra::values(rmse_r) <- rmse_v
  mae_r <- terra::rast(template)
  terra::values(mae_r) <- mae_v
  bias_r <- terra::rast(template)
  terra::values(bias_r) <- bias_v
  r_r <- terra::rast(template)
  terra::values(r_r) <- r_v
  ns_r <- terra::rast(template)
  terra::values(ns_r) <- ns_v
  out <- c(rmse_r, mae_r, bias_r, r_r, ns_r)
  names(out) <- c("rmse", "mae", "bias", "r", "n_scenes")
  return(out)
}

#' Retrieve scenes for fusion validation
#'
#' Pass the selected index, cloud settings, parameters and dimensions to
#' \code{get_rs()}. Standalone calls require at least one dimension.
#' Piped calls reuse dimensions stored on their input.
#' @noRd
.validateBuildRs <- function(longitude, latitude, start_date, end_date,
                             x_metres, y_metres, params,
                             max_cloud_cover, scl_classes,
                             index_name = "EVI2", source = "auto"){

  if( is.null(params) ){
    params <- rs_params()
  }
  if( identical(source, "auto") ){
    source <- if( as.Date(start_date) >= .L8_START ) { "landsat-8" } else { "landsat-5" }
  }
  return(get_rs(longitude = longitude,
                latitude = latitude,
                start_date = start_date,
                end_date = end_date,
                index_name = index_name,
                source = source,
                x_metres = x_metres,
                y_metres = y_metres,
                max_cloud_cover = max_cloud_cover,
                scl_classes = scl_classes,
                params = params))
}

#' Assess Landsat-MODIS fusion using leave-one-out agreement
#'
#' Holds out each available Landsat scene in turn, estimates its spectral index
#' values from the nearest remaining Landsat anchor and MCD43A4 NBAR data, and
#' compares the estimate with the held-out Landsat observation. A held-out
#' scene is skipped when it or its anchor has no MCD43A4 scene of its own within
#' one day. The resulting metrics describe agreement for the selected data,
#' place, and period. They do not establish general accuracy or replace
#' independent validation.
#'
#' \strong{Modes:}
#' \itemize{
#'   \item Standalone \code{validate_data_fusion(longitude, latitude,
#'         start_date, end_date, ...)}: returns a \code{data.frame} with
#'         agreement metrics for each held-out scene (one row per scene).
#'   \item Piped data.frame (after \code{as_time_series()}): appends
#'         agreement metric columns to the input data.frame, matched
#'         on daily \code{date}. Monthly/yearly bucket keys cannot be matched
#'         and receive \code{NA} metrics with a warning.
#'   \item Piped SpatRaster (after \code{as_raster()}): returns a
#'         five-layer SpatRaster of per-pixel agreement metrics
#'         (\code{rmse}, \code{mae}, \code{bias}, \code{r},
#'         \code{n_scenes}) on the input grid.
#' }
#'
#' \strong{Metrics:} RMSE, MAE, bias (signed mean of \code{pred -
#' obs}), and Pearson \code{r}. Per-scene metrics aggregate over
#' pixels. Per-pixel metrics aggregate over held-out scenes. Pixels (or
#' scenes) with fewer than \code{min_valid_n} valid pairs receive
#' \code{NA} for all four metrics. The count column or layer
#' (\code{n_pixels} or \code{n_scenes}) is always populated so the
#' user can audit data sufficiency.
#'
#' @param x A \code{SpatRaster}, a \code{data.frame} with daily dates, or numeric
#'   longitude. Convert an \code{ocular} object with \code{as_raster()} or
#'   \code{as_time_series()} first.
#' @param latitude,start_date,end_date Required for standalone form. Otherwise,
#'   these arguments are ignored and values are read from \code{rs_meta}.
#' @param x_metres,y_metres Assessment-area dimensions. At least
#'   one is required for the standalone form. Piped modes use the
#'   dimensions recorded on the input.
#' @param params Optional validated parameter list. In standalone mode,
#'   \code{NULL} uses \code{rs_params()} defaults. Piped modes instead reuse
#'   the parameters recorded on the input unless this argument is supplied.
#' @param index_name Spectral index assessed in the standalone form
#'   (default \code{"EVI2"}). Fusion supports NDVI, EVI, DVI, kNDVI, BNDVI,
#'   NIRv, SAVI, MSAVI, SR, EVI2, OSAVI, ARVI, NBR, and BAI. Piped modes take
#'   the index from the input's \code{rs_meta} (recorded by
#'   \code{as_time_series()} or \code{as_raster()}), so agreement is calculated
#'   for the index that was produced. MCD43A4 product-specific quality layers
#'   are not yet applied.
#' @param source Landsat source for the standalone assessment: \code{"auto"}
#'   (default) selects Landsat 8/9 from 2013 onward and Landsat 5/7 earlier.
#'   Alternatives are \code{"landsat-8"} and \code{"landsat-5"}. Piped modes use the
#'   source recorded on their input and reject non-Landsat inputs.
#' @param min_valid_n Minimum number of valid (non-NA) paired
#'   values required for metric computation. The default is \code{4L}. Pearson
#'   correlation itself requires at least two pairs. Pixels or scenes below
#'   this threshold receive \code{NA} metrics. Raise it for stricter analysis.
#' @param max_cloud_cover STAC cloud-cover ceiling (default 50).
#' @param scl_classes Optional Landsat \code{QA_PIXEL} bit positions to
#'   exclude, passed to \code{get_rs()}. \code{NULL} applies no per-pixel bit
#'   mask.
#' @returns A \code{data.frame} for standalone or time-series input, or a
#'   \code{SpatRaster} for raster input. If no held-out scene can be estimated,
#'   standalone calls return \code{NULL} and other forms return the input
#'   unchanged. Input without \code{rs_meta} is returned unchanged with a warning.
#' @export
validate_data_fusion <- function(x,
                                 latitude = NULL,
                                 start_date = NULL,
                                 end_date = NULL,
                                 x_metres = NULL,
                                 y_metres = NULL,
                                 params = NULL,
                                 min_valid_n = 4L,
                                 max_cloud_cover = 50,
                                 scl_classes = NULL,
                                 index_name = "EVI2",
                                 source = "auto"){

  cloud_was_missing <- missing(max_cloud_cover)
  scl_was_missing <- missing(scl_classes)
  if( !.isWholeNumber(min_valid_n) || min_valid_n < 1L ){
    stop("validate_data_fusion: min_valid_n must be a positive integer.",
         call. = FALSE)
  }
  min_valid_n <- as.integer(min_valid_n)
  if( !is.character(source) || length(source) != 1L || is.na(source) ){
    stop("validate_data_fusion: source must be auto, landsat-8, or landsat-5.",
         call. = FALSE)
  }
  source <- match.arg(source, c("auto", "landsat-8", "landsat-5"))

  if( is_rs(x) ){
    stop("validate_data_fusion: ocular objects are not supported directly. ",
         "Use as_raster() or daily as_time_series() output as the input.", call. = FALSE)
  }

  if( inherits(x, "SpatRaster") ){
    meta <- attr(x, "rs_meta")
    if( is.null(meta) ){
      warning("validate_data_fusion: SpatRaster lacks rs_meta.", call. = FALSE)
      return(x)
    }
    if( !meta$source %in% c("landsat-8", "landsat-5") ){
      stop("validate_data_fusion: piped input must originate from Landsat. ",
           "Found '", meta$source, "'.", call. = FALSE)
    }
    params_use <- params %||% meta$params
    cloud_use <- if( cloud_was_missing ) { meta$max_cloud_cover } else { max_cloud_cover }
    scl_use <- if( scl_was_missing ) { meta$scl_classes } else { scl_classes }
    rs <- .validateBuildRs(longitude = meta$longitude,
                           latitude = meta$latitude,
                           start_date = meta$start_date,
                           end_date = meta$end_date,
                           x_metres = meta$x_metres,
                           y_metres = meta$y_metres,
                           params = params_use,
                           max_cloud_cover = cloud_use,
                           scl_classes = scl_use,
                           index_name = meta$index_name %||% "EVI2",
                           source = meta$source)
    loo <- .validateCompute(rs, params_use)
    if( is.null(loo) ){
      message("validate_data_fusion: insufficient scenes for validation.")
      return(x)
    }
    return(.collapseToPixel(loo, template = x[[1L]],
                            min_valid_n = min_valid_n))
  }

  if( is.data.frame(x) ){
    meta <- attr(x, "rs_meta")
    if( is.null(meta) ){
      warning("validate_data_fusion: data.frame lacks rs_meta.", call. = FALSE)
      return(x)
    }
    if( !"date" %in% names(x) ){
      stop("validate_data_fusion: piped data.frame must contain `date`.",
           call. = FALSE)
    }
    if( !meta$source %in% c("landsat-8", "landsat-5") ){
      stop("validate_data_fusion: piped input must originate from Landsat. ",
           "Found '", meta$source, "'.", call. = FALSE)
    }
    params_use <- params %||% meta$params
    cloud_use <- if( cloud_was_missing ) { meta$max_cloud_cover } else { max_cloud_cover }
    scl_use <- if( scl_was_missing ) { meta$scl_classes } else { scl_classes }
    rs <- .validateBuildRs(longitude = meta$longitude,
                           latitude = meta$latitude,
                           start_date = meta$start_date,
                           end_date = meta$end_date,
                           x_metres = meta$x_metres,
                           y_metres = meta$y_metres,
                           params = params_use,
                           max_cloud_cover = cloud_use,
                           scl_classes = scl_use,
                           index_name = meta$index_name %||% "EVI2",
                           source = meta$source)
    loo <- .validateCompute(rs, params_use)
    if( is.null(loo) ){
      return(x)
    }
    vres <- .collapseToScene(loo, min_valid_n = min_valid_n)
    # Metrics use daily dates. Warn when month or year labels prevent a match.
    midx <- match(as.character(x$date), vres$date)
    if( all(is.na(midx)) ){
      warning("validate_data_fusion: no input dates matched the per-scene ",
              "validation dates. This mode expects daily input. Frames aggregated ",
              "by month or year cannot be matched. Metric ",
              "columns will be NA.", call. = FALSE)
    }
    x$rmse <- vres$rmse[midx]
    x$mae <- vres$mae[midx]
    x$bias <- vres$bias[midx]
    x$r <- vres$r[midx]
    x$n_pixels <- vres$n_pixels[midx]
    return(x)
  }

  if( is.null(latitude) || is.null(start_date) || is.null(end_date) ){
    stop("validate_data_fusion standalone form requires latitude, ",
         "start_date, end_date.", call. = FALSE)
  }
  if( is.null(x_metres) && is.null(y_metres) ){
    stop("validate_data_fusion standalone form requires at least one ",
         "of x_metres or y_metres.", call. = FALSE)
  }
  rs <- .validateBuildRs(longitude = x,
                         latitude = latitude,
                         start_date = start_date,
                         end_date = end_date,
                         x_metres = x_metres,
                         y_metres = y_metres,
                         params = params,
                         max_cloud_cover = max_cloud_cover,
                         scl_classes = scl_classes,
                         index_name = index_name,
                         source = source)
  loo <- .validateCompute(rs, params)
  if( is.null(loo) ){
    return(NULL)
  }
  return(.collapseToScene(loo, min_valid_n = min_valid_n))
}
