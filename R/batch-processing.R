# =========================================================================
# Batch retrieval, caching and retries
# =========================================================================
# batch_rs() visits unfinished sites in rounds. Each visit retrieves once
# and processes a successful retrieval once.

.RS_CACHE_VERSION <- "v2"


# =========================================================================
# Batch helpers
# =========================================================================

#' Versioned scene cache for the current user
#' @noRd
.rsCacheDir <- function(create = TRUE){

  base <- tools::R_user_dir("ocular", which = "cache")
  d <- file.path(base, "rs", .RS_CACHE_VERSION)
  if( isTRUE(create) ){
    d <- .initialiseCacheDir(d, "rs", .RS_CACHE_VERSION, ".rsCacheDir()")
  }
  return(d)
}

#' Cache key for a complete batch request
#' @noRd
.rsCacheKey <- function(longitude, latitude, start_date, end_date, index_name,
                        source, x_metres, y_metres, max_cloud_cover, params,
                        retries, cloud_relax, pipeline_cache_key,
                        package_version = as.character(utils::packageVersion("ocular"))){

  # Use one numeric type so equivalent requests share a key. Include the package
  # version to keep results from different versions separate.
  query <- list(
    longitude = round(longitude, 7L),
    latitude = round(latitude, 7L),
    start_date = as.character(start_date),
    end_date = as.character(end_date),
    index_name = index_name,
    source = source,
    x_metres = as.numeric(x_metres),
    y_metres = as.numeric(y_metres),
    max_cloud_cover = as.numeric(max_cloud_cover),
    params = params,
    retries = retries,
    retry_policy = "rounds-transient-v1",
    relaxed_final_attempt = is.finite(retries) && retries > 1L,
    cloud_relax = as.numeric(cloud_relax),
    pipeline_cache_key = pipeline_cache_key,
    package_version = package_version)
  return(paste0(substr(.stableHash(query), 1L, 24L), ".rds"))
}

#' Validate and normalise one batch site without fetching
#' @noRd
.validateBatchSite <- function(site){

  id <- as.character(site$id[[1L]])
  if( length(id) != 1L || is.na(id) || !nzchar(id) ){
    stop("site id must be non-missing and non-empty.", call. = FALSE)
  }
  longitude <- site$longitude[[1L]]
  latitude <- site$latitude[[1L]]
  if( !.isScalarNumber(longitude) || longitude < -180 || longitude > 180 ){
    stop("longitude must be a finite numeric in [-180, 180].", call. = FALSE)
  }
  if( !.isScalarNumber(latitude) || latitude < -80 || latitude > 84 ){
    stop("latitude must be a finite numeric in [-80, 84].", call. = FALSE)
  }
  start_date <- .normaliseOptionalDate(site$start_date[[1L]], "start_date")
  end_date <- .normaliseOptionalDate(site$end_date[[1L]], "end_date")
  if( as.Date(start_date) > as.Date(end_date) ){
    stop("start_date must be on or before end_date.", call. = FALSE)
  }
  if( as.numeric(as.Date(end_date) - as.Date(start_date)) < 2 ){
    stop("date range must span at least two days.", call. = FALSE)
  }
  return(list(id = id, longitude = longitude, latitude = latitude,
              start_date = start_date, end_date = end_date))
}

#' Validate daily dates or year/month bucket keys returned by a batch pipeline
#' @noRd
.normaliseBatchDates <- function(x, n){

  if( length(x) != n || anyNA(x) ){
    return(NULL)
  }
  if( inherits(x, "Date") ){
    return(x)
  }
  value <- as.character(x)
  if( anyNA(value) || any(!nzchar(value)) ){
    return(NULL)
  }
  full <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", value)
  if( all(full) ){
    out <- suppressWarnings(as.Date(value))
    if( anyNA(out) ){
      return(NULL)
    }
    return(out)
  }
  month <- grepl("^[0-9]{4}-(0[1-9]|1[0-2])$", value)
  year <- grepl("^[0-9]{4}$", value)
  if( all(month | year) ){
    return(value)
  }
  return(NULL)
}

#' Read one batch cache entry, treating corruption as a miss
#' @noRd
.rsCacheRead <- function(path){

  obj <- tryCatch(readRDS(path), error = function(e){
    return(NULL)
  })
  if( !is.list(obj) || !is.data.frame(obj$data) ||
      !all(c("source", "date", "value") %in% names(obj$data)) ){
    return(NULL)
  }
  if( !is.character(obj$data$source) ||
      length(obj$data$source) != nrow(obj$data) ||
      !is.numeric(obj$data$value) || !is.null(dim(obj$data$value)) ||
      length(obj$data$value) != nrow(obj$data) ){
    return(NULL)
  }
  obj$data$date <- .normaliseBatchDates(obj$data$date, nrow(obj$data))
  if( is.null(obj$data$date) ){
    return(NULL)
  }
  return(obj)
}

#' Clear the on-disk scene-series cache
#'
#' Removes only ocular scene-series cache entries from the default cache or an
#' ownership-marked custom directory. It does not recursively remove the
#' directory or delete unrelated files. Cleared results are retrieved and
#' processed again when they are next requested.
#'
#' @param cache_dir Cache directory to clear. \code{NULL} selects the versioned
#'   per-user cache. A custom directory must previously have been initialised
#'   by a cache-enabled \code{batch_rs()} call. Unmarked custom directories and
#'   broad paths are refused.
#' @returns \code{TRUE}, invisibly.
#' @seealso \code{\link{batch_rs}}, \code{\link{clear_ftw_cache}}
#' @export
clear_rs_cache <- function(cache_dir = NULL){

  if( !is.null(cache_dir) &&
      (!is.character(cache_dir) || length(cache_dir) != 1L ||
       is.na(cache_dir) || !nzchar(cache_dir)) ){
    stop("clear_rs_cache(): `cache_dir` must be NULL or one non-empty path.",
         call. = FALSE)
  }
  default_dir <- .rsCacheDir(create = FALSE)
  if( is.null(cache_dir) ){
    cache_dir <- default_dir
  }
  cache_dir <- .validateCacheClearDir(
    cache_dir, "rs", .RS_CACHE_VERSION, "clear_rs_cache()", default_dir)
  entries <- .cacheEntryFiles(cache_dir, "^[0-9a-f]{24}[.]rds$")
  if( length(entries) > 0L ){
    status <- unlink(entries, recursive = FALSE, force = FALSE)
    if( !identical(status, 0L) ){
      warning("clear_rs_cache(): one or more cache entries could not be removed.",
              call. = FALSE)
    }
  }
  return(invisible(TRUE))
}

#' Reduce a delineated site to a time series
#'
#' Use the union of its areas so each site and date has one value.
#' @noRd
.defaultBatchPipeline <- function(rs){

  out <- as_time_series(boundary_delineation(rs))
  if( !"value" %in% names(out) && "union" %in% names(out) ){
    out$value <- out$union
  }
  return(out)
}


# =========================================================================
# Site retrieval helpers
# =========================================================================

#' Retrieve once within a batch visit
#' @noRd
.batchGetRs <- function(...){

  old <- options(ocular.batch_single_attempt = TRUE)
  on.exit(options(old), add = TRUE)
  return(get_rs(...))
}

#' Wait between rounds, capped at one minute
#' @noRd
.batchPause <- function(round){

  Sys.sleep(min(60, 2^min(round, 6)))
  return(invisible(NULL))
}

#' Retrieve and process one site
#'
#' With a finite retry budget, relax the cloud limit on the last attempt only
#' if more than one attempt was requested. Unlimited retries keep the requested
#' limit. Retry transient retrieval failures, but run the pipeline only once.
#' @noRd
.batchFetchOne <- function(longitude, latitude, start_date, end_date, index_name,
                           source, x_metres, y_metres, max_cloud_cover, params,
                           pipeline, retries, cloud_relax, attempt = 1){

  relax <- is.finite(retries) && retries > 1L && attempt == retries
  cc <- if( relax ) { max(max_cloud_cover, cloud_relax) } else { max_cloud_cover }
  rs <- tryCatch(
    .batchGetRs(longitude = longitude, latitude = latitude,
                start_date = as.character(start_date),
                end_date = as.character(end_date), index_name = index_name,
                source = source, x_metres = x_metres, y_metres = y_metres,
                max_cloud_cover = cc, params = params),
    error = identity)
  if( inherits(rs, "error") ){
    return(errorCondition(.stripUrlQuery(conditionMessage(rs)),
                           class = "ocular_batch_error", call = NULL,
                           retryable = inherits(rs, "ocular_stac_error") &&
                             isTRUE(rs$retryable)))
  }

  # A pipeline error is not a transient retrieval failure.
  d <- tryCatch(pipeline(rs), error = identity)
  if( inherits(d, "error") ){
    return(simpleError(.stripUrlQuery(conditionMessage(d))))
  }
  if( !is.data.frame(d) || !all(c("date", "value") %in% names(d)) ){
    return(simpleError("the pipeline must return a data.frame with `date` and `value`."))
  }
  dates <- .normaliseBatchDates(d$date, nrow(d))
  if( is.null(dates) ){
    return(simpleError("the pipeline returned invalid dates/bucket keys."))
  }
  if( !is.numeric(d$value) || !is.null(dim(d$value)) ||
      length(d$value) != nrow(d) ){
    return(simpleError("the pipeline must return `value` as one numeric vector."))
  }
  return(list(data = data.frame(source = rep(rs$geom$source, nrow(d)), date = dates,
                                value = d$value, stringsAsFactors = FALSE),
              n_scenes = length(rs$scenes), cloud = cc))
}


# =========================================================================
# Batch processing
# =========================================================================

#' Build spectral index time series for multiple sites
#'
#' For each row of \code{sites}, retrieves scenes with \code{get_rs()} and
#' passes the resulting \code{ocular} object to a post-retrieval function. The
#' default function runs field boundary delineation and creates a time series.
#' Per-site results are combined in one long \code{data.frame}. Completed
#' results can be cached on disk, and failed retrievals can be attempted again.
#'
#' @param sites A \code{data.frame} with columns \code{id}, \code{longitude},
#'   \code{latitude}, \code{start_date}, and \code{end_date}, with one row per
#'   site and date range.
#' @param index_name Index to calculate. See \code{get_rs()}.
#' @param source Image source passed to \code{get_rs()}. The default is
#'   \code{"auto"}.
#' @param x_metres,y_metres Width and height of the analysis area in metres.
#'   If \code{y_metres} is \code{NULL}, \code{x_metres} is used for both.
#' @param max_cloud_cover Maximum scene-level cloud cover for the ordinary
#'   retrieval attempts.
#' @param params Optional complete or partial list of settings accepted by
#'   \code{rs_params()}, passed to each \code{get_rs()} call.
#' @param pipeline Post-retrieval function that accepts an \code{ocular} object
#'   and returns a \code{data.frame} with \code{date} and \code{value}. By
#'   default, \code{boundary_delineation()} is followed by
#'   \code{as_time_series()}. If \code{multiple_areas = TRUE}, the default
#'   function uses the \code{union} column as \code{value}.
#' @param cache Whether to read and write the on-disk cache. See
#'   \code{\link{clear_rs_cache}}.
#' @param cache_dir Cache location. \code{NULL} uses the versioned per-user
#'   cache directory.
#' @param pipeline_cache_key Stable character identifier for a custom
#'   \code{pipeline}. A custom function is not cached unless this is supplied,
#'   because its body does not identify captured values reliably.
#' @param retries Maximum retrieval attempts per site, or \code{Inf} for no
#'   limit. Unfinished sites are visited in input order once per round. Only
#'   classified transient STAC query or transport failures are retried.
#'   Invalid requests, no scenes, unclassified errors, and pipeline failures
#'   are final. Asset-reading failures are not retried by the batch schedule.
#'   Rounds wait 2, 4, 8, 16, 32, then 60 seconds
#'   between retries. Each visit makes one STAC query attempt, including any
#'   pagination and signing. The pipeline runs once after retrieval succeeds.
#' @param cloud_relax Cloud-cover ceiling considered on the final attempt of a
#'   finite budget when \code{retries > 1}. The final ceiling is the greater of
#'   this value and \code{max_cloud_cover}. Unlimited retries keep
#'   \code{max_cloud_cover}. A no-scenes result is final, so it does not advance
#'   to a later attempt with a relaxed ceiling.
#' @param on_error \code{"skip"} (default) records a final failure and continues.
#'   \code{"stop"} aborts when a failure is final or its attempts are exhausted.
#'   An interrupt returns completed results for either setting.
#' @param quiet Whether to suppress per-site progress messages.
#' @returns A \code{data.frame} with \code{id}, \code{index}, \code{source},
#'   \code{date}, and \code{value}. Failed sites are absent and recorded in the
#'   \code{"failures"} attribute with their \code{id}, \code{index},
#'   \code{reason}, \code{attempts}, and \code{last_error}. Attempts count
#'   calls to \code{get_rs()}, including one interrupted during retrieval or
#'   its pipeline. If interrupted, completed rows are returned and unfinished
#'   sites are recorded as interrupted, retaining their last retrieval error.
#'   Results and failures follow the input site order.
#'   The \code{"outcomes"} attribute has one row per input site, in input
#'   order, including repeated IDs. Its columns are \code{id}, \code{index},
#'   \code{status}, \code{dates}, \code{attempts}, and \code{cached}. Status is
#'   \code{"success"}, \code{"empty"}, \code{"failure"}, or
#'   \code{"interrupted"}. Empty successes produce no result rows or failures.
#'   \code{dates} counts returned date rows, with zero for an empty success
#'   and \code{NA} for failure or interruption. \code{attempts} counts
#'   retrieval calls in this invocation, so a cache-only hit has zero attempts.
#'   \code{cached} is true only for a completed result read from cache.
#' @seealso \code{\link{get_rs}}, \code{\link{as_time_series}},
#'   \code{\link{clear_rs_cache}}
#' @examples
#' \dontrun{
#' sites <- data.frame(id = c("a", "b"),
#'                     longitude = c(-97.29971, -99.17693),
#'                     latitude  = c(30.60044, 29.40344),
#'                     start_date = c("2019-03-26", "2019-03-11"),
#'                     end_date   = c("2019-08-20", "2019-08-05"))
#' ts <- batch_rs(sites, index_name = "GNDVI", x_metres = 1200)
#' attr(ts, "outcomes")
#' attr(ts, "failures")
#' }
#' @export
batch_rs <- function(sites, index_name = "EVI2", source = "auto",
                     x_metres = 1000L, y_metres = NULL,
                     max_cloud_cover = 50L, params = NULL, pipeline = NULL,
                     cache = TRUE, cache_dir = NULL, pipeline_cache_key = NULL,
                     retries = 3L,
                     cloud_relax = 80L, on_error = c("skip", "stop"),
                     quiet = FALSE){

  on_error <- match.arg(on_error)
  need <- c("id", "longitude", "latitude", "start_date", "end_date")
  if( !is.data.frame(sites) || !all(need %in% names(sites)) ){
    stop("batch_rs(): `sites` must be a data.frame with columns ",
         paste(need, collapse = ", "), ". Please add them and try again.",
         call. = FALSE)
  }
  if( nrow(sites) == 0L ){
    stop("batch_rs(): `sites` has no rows. Please supply at least one site.",
         call. = FALSE)
  }
  if( !is.character(index_name) || length(index_name) != 1L ||
      is.na(index_name) || !nzchar(index_name) ||
      !index_name %in% unique(c(names(s2_index_list), names(landsat_index_list))) ){
    stop("batch_rs(): `index_name` is not a supported fine-source index.",
         call. = FALSE)
  }
  if( !is.character(source) || length(source) != 1L || is.na(source) ){
    stop("batch_rs(): `source` must be auto, sentinel-2, landsat-8, or landsat-5.",
         call. = FALSE)
  }
  source <- match.arg(source,
                      c("auto", "sentinel-2", "landsat-8", "landsat-5"))
  if( !.isScalarNumber(x_metres) || x_metres <= 0 ){
    stop("batch_rs(): `x_metres` must be a finite positive numeric.",
         call. = FALSE)
  }
  if( is.null(y_metres) ){
    y_metres <- x_metres
  }
  if( !.isScalarNumber(y_metres) || y_metres <= 0 ){
    stop("batch_rs(): `y_metres` must be NULL or a finite positive numeric.",
         call. = FALSE)
  }
  if( !.isScalarNumber(max_cloud_cover) || max_cloud_cover < 0 ||
      max_cloud_cover > 100 ){
    stop("batch_rs(): `max_cloud_cover` must be in [0, 100].", call. = FALSE)
  }
  if( !.isScalarNumber(cloud_relax) || cloud_relax < 0 || cloud_relax > 100 ){
    stop("batch_rs(): `cloud_relax` must be in [0, 100].", call. = FALSE)
  }
  unlimited <- is.numeric(retries) && !is.complex(retries) &&
    length(retries) == 1L && !is.na(retries) && retries == Inf
  if( is.complex(retries) ||
      (!unlimited && (!.isWholeNumber(retries) || retries < 1L)) ){
    stop("batch_rs(): `retries` must be a single positive integer or Inf.",
         call. = FALSE)
  }
  retries <- if( unlimited ) { Inf } else { as.integer(retries) }
  if( !is.logical(cache) || length(cache) != 1L || is.na(cache) ){
    stop("batch_rs(): `cache` must be TRUE or FALSE.", call. = FALSE)
  }
  if( !is.logical(quiet) || length(quiet) != 1L || is.na(quiet) ){
    stop("batch_rs(): `quiet` must be TRUE or FALSE.", call. = FALSE)
  }
  if( !is.null(cache_dir) &&
      (!is.character(cache_dir) || length(cache_dir) != 1L ||
       is.na(cache_dir) || !nzchar(cache_dir)) ){
    stop("batch_rs(): `cache_dir` must be NULL or one non-empty path.",
         call. = FALSE)
  }
  if( !is.null(pipeline_cache_key) &&
      (!is.character(pipeline_cache_key) || length(pipeline_cache_key) != 1L ||
       is.na(pipeline_cache_key) || !nzchar(pipeline_cache_key)) ){
    stop("batch_rs(): `pipeline_cache_key` must be NULL or one non-empty string.",
         call. = FALSE)
  }

  if( is.null(params) ){
    params <- rs_params()
  }else{
    if( !is.list(params) || is.null(names(params)) || anyNA(names(params)) ||
        any(!nzchar(names(params))) || anyDuplicated(names(params)) ){
      stop("batch_rs(): `params` must be a uniquely named list.", call. = FALSE)
    }
    param_user_set <- .paramUserSet(params, "batch_rs()")
    params <- do.call(rs_params, params)
    attr(params, "user_set") <- param_user_set
  }

  default_pipeline <- is.null(pipeline)
  if( default_pipeline ){
    pipeline <- .defaultBatchPipeline
    pipeline_cache_key <- "default-boundary-timeseries-v2"
  }else if( !is.function(pipeline) ){
    stop("batch_rs(): `pipeline` must be NULL or a function.", call. = FALSE)
  }else if( isTRUE(cache) && is.null(pipeline_cache_key) ){
    warning("batch_rs(): caching disabled for a custom pipeline without ",
            "`pipeline_cache_key`.", call. = FALSE)
    cache <- FALSE
  }
  n <- nrow(sites)
  rows <- vector("list", n)
  fails <- vector("list", n)
  checked <- vector("list", n)
  attempts <- numeric(n)
  cached <- logical(n)
  last_error <- rep("", n)
  pending <- rep(TRUE, n)
  ids <- vapply(seq_len(n), function(i){
    id <- tryCatch(as.character(sites$id[[i]]), error = function(e){
      return(NULL)
    })
    return(if( length(id) == 1L ) { id } else { paste0("row ", i) })
  }, character(1L))

  # Return completed results on interruption. With on_error = "stop",
  # failures propagate.
  interrupted <- tryCatch({
    dir_use <- if( !isTRUE(cache) ){
      NULL
    }else if( is.null(cache_dir) ){
      .rsCacheDir()
    }else{
      .initialiseCacheDir(cache_dir, "rs", .RS_CACHE_VERSION, "batch_rs()")
    }
    round <- 1
    while( any(pending) ){
      for( i in which(pending) ){
        if( is.null(checked[[i]]) ){
          checked[i] <- list(tryCatch(
            .validateBatchSite(sites[i, , drop = FALSE]), error = identity))
        }
        s <- checked[[i]]
        if( inherits(s, "error") ){
          got <- s
        }else{
          key <- .rsCacheKey(s$longitude, s$latitude, s$start_date, s$end_date,
                             index_name, source, x_metres, y_metres,
                             max_cloud_cover, params, retries, cloud_relax,
                             pipeline_cache_key)
          f <- if( is.null(dir_use) ) { NULL } else { file.path(dir_use, key) }
          hit <- NULL
          if( !is.null(f) && file.exists(f) ){
            hit <- .rsCacheRead(f)
            if( is.null(hit) && attempts[[i]] == 0 && !isTRUE(quiet) ){
              message("  batch_rs: ignoring unreadable cache entry for ", s$id)
            }
          }
          if( !is.null(hit) ){
            cached[[i]] <- TRUE
            rows[[i]] <- cbind(id = rep(s$id, nrow(hit$data)),
                               index = rep(index_name, nrow(hit$data)), hit$data,
                               stringsAsFactors = FALSE)
            pending[[i]] <- FALSE
            if( !isTRUE(quiet) ){
              message("  batch_rs: ", s$id, " / ", index_name, " (cached, ",
                      nrow(hit$data), " dates)")
            }
            next
          }
          attempts[[i]] <- attempts[[i]] + 1
          got <- tryCatch(
            .batchFetchOne(s$longitude, s$latitude, s$start_date, s$end_date,
                           index_name, source, x_metres, y_metres,
                           max_cloud_cover, params, pipeline, retries,
                           cloud_relax, attempt = attempts[[i]]),
            error = identity)
        }
        if( inherits(got, "error") ){
          msg <- .stripUrlQuery(conditionMessage(got))
          last_error[[i]] <- msg
          if( inherits(got, "ocular_batch_error") && isTRUE(got$retryable) &&
              attempts[[i]] < retries ){
            next
          }
          fails[[i]] <- data.frame(id = ids[[i]], index = index_name,
                                    reason = msg, attempts = attempts[[i]],
                                    last_error = msg, stringsAsFactors = FALSE)
          pending[[i]] <- FALSE
          if( identical(on_error, "stop") ){
            stop("batch_rs(): ", ids[[i]], " / ", index_name, ": ", msg,
                 call. = FALSE)
          }
          if( !isTRUE(quiet) ){
            message("  batch_rs: skipped ", ids[[i]], " -- ", msg)
          }
          next
        }

        # Keep the completed result if cache writing or progress output is interrupted.
        rows[[i]] <- cbind(id = rep(s$id, nrow(got$data)),
                           index = rep(index_name, nrow(got$data)), got$data,
                           stringsAsFactors = FALSE)
        pending[[i]] <- FALSE
        if( !is.null(f) ){
          cache_obj <- list(data = got$data,
                            meta = list(key = key, created_at = Sys.time(),
                                        cloud_used = got$cloud,
                                        attempts = attempts[[i]],
                                        retry_policy = "rounds-transient-v1"))
          tryCatch(saveRDS(cache_obj, f), error = function(e){
            warning("batch_rs(): cache write failed for ", s$id, ": ",
                    .stripUrlQuery(conditionMessage(e)), call. = FALSE)
          })
        }
        if( !isTRUE(quiet) ){
          message("  batch_rs: ", s$id, " / ", index_name, " -- ", got$n_scenes,
                  " scenes (", nrow(got$data), " dates)")
        }
      }
      if( any(pending) ){
        .batchPause(round)
        round <- round + 1
      }
    }
    FALSE
  }, interrupt = function(e){
    return(TRUE)
  })
  if( interrupted ){
    unfinished <- vapply(rows, is.null, logical(1L)) &
      vapply(fails, is.null, logical(1L))
    for( i in which(unfinished) ){
      fails[[i]] <- data.frame(id = ids[[i]], index = index_name,
                                reason = "Interrupted before completion.",
                                attempts = attempts[[i]],
                                last_error = last_error[[i]],
                                stringsAsFactors = FALSE)
    }
  }
  # An interrupted cache hit is complete only after its rows have been stored.
  dates <- vapply(rows, function(x){
    return(if( is.null(x) ) { NA_integer_ } else { nrow(x) })
  }, integer(1L))
  completed <- !is.na(dates)
  status <- rep("failure", n)
  status[completed] <- "success"
  status[completed & dates == 0L] <- "empty"
  if( interrupted ){
    status[unfinished] <- "interrupted"
  }
  outcomes <- data.frame(id = ids, index = rep(index_name, n), status = status,
                          dates = dates, attempts = attempts,
                          cached = completed & cached, stringsAsFactors = FALSE)
  rows <- Filter(Negate(is.null), rows)
  fails <- Filter(Negate(is.null), fails)

  out <- if( length(rows) > 0L ){
    do.call(rbind, rows)
  }else{
    data.frame(id = character(0L), index = character(0L), source = character(0L),
               date = as.Date(character(0L)), value = numeric(0L),
               stringsAsFactors = FALSE)
  }
  rownames(out) <- NULL
  attr(out, "failures") <- if( length(fails) > 0L ) { do.call(rbind, fails) } else { NULL }
  attr(out, "outcomes") <- outcomes
  return(out)
}
