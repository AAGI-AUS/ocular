# =========================================================================
# Batch processing, caching, and retry handling
#
# Batch requests can retry transient retrieval failures and cache completed
# results so equivalent requests do not repeat the same retrieval work.
#
# The cache key is derived from the complete result-affecting query -- point,
# dates, source, index, dimensions, cloud behaviour, params, and pipeline
# identity -- never from a caller's row label, so equivalent requests share an
# entry. Layout and lifecycle follow the FTW cache (per-user, versioned,
# cleared by a matching clear_* verb).
#
# .batchFetchOne applies the retry and cloud-relaxation policy for batch calls.
# =========================================================================

.RS_CACHE_VERSION <- "v2"


# ---- cache primitives ---------------------------------------------------

#' On-disk scene-series cache directory (per-user, versioned)
#' @noRd
.rsCacheDir <- function(create = TRUE){
  base <- tools::R_user_dir("ocular", which = "cache")
  d <- file.path(base, "rs", .RS_CACHE_VERSION)
  if( isTRUE(create) )
    d <- .initialiseCacheDir(d, "rs", .RS_CACHE_VERSION, ".rsCacheDir()")
  return(d)
}

#' Safe content-hash cache key for one complete batch query.
#' @noRd
.rsCacheKey <- function(longitude, latitude, start_date, end_date, index_name,
                        source, x_metres, y_metres, max_cloud_cover, params,
                        retries, cloud_relax, pipeline_cache_key){
  query <- list(
    longitude = round(longitude, 7L),
    latitude = round(latitude, 7L),
    start_date = as.character(start_date),
    end_date = as.character(end_date),
    index_name = index_name,
    source = source,
    x_metres = x_metres,
    y_metres = y_metres,
    max_cloud_cover = max_cloud_cover,
    params = params,
    retries = retries,
    relaxed_final_attempt = retries > 1L,
    cloud_relax = cloud_relax,
    pipeline_cache_key = pipeline_cache_key)
  paste0(substr(.stableHash(query), 1L, 24L), ".rds")
}

#' Validate and normalize one batch site without fetching (internal)
#' @noRd
.validateBatchSite <- function(site){
  id <- as.character(site$id[[1L]])
  if( length(id) != 1L || is.na(id) || !nzchar(id) )
    stop("site id must be non-missing and non-empty.", call. = FALSE)
  longitude <- site$longitude[[1L]]
  latitude <- site$latitude[[1L]]
  if( !.isScalarNumber(longitude) || longitude < -180 || longitude > 180 )
    stop("longitude must be a finite numeric in [-180, 180].", call. = FALSE)
  if( !.isScalarNumber(latitude) || latitude < -80 || latitude > 84 )
    stop("latitude must be a finite numeric in [-80, 84].", call. = FALSE)
  start_date <- .normaliseOptionalDate(site$start_date[[1L]], "start_date")
  end_date <- .normaliseOptionalDate(site$end_date[[1L]], "end_date")
  if( as.Date(start_date) > as.Date(end_date) )
    stop("start_date must be on or before end_date.", call. = FALSE)
  if( as.numeric(as.Date(end_date) - as.Date(start_date)) < 2 )
    stop("date range must span at least two days.", call. = FALSE)
  list(id = id, longitude = longitude, latitude = latitude,
       start_date = start_date, end_date = end_date)
}

#' Validate daily dates or year/month bucket keys returned by a batch pipeline
#' @noRd
.normaliseBatchDates <- function(x, n){
  if( length(x) != n || anyNA(x) ) return(NULL)
  if( inherits(x, "Date") ) return(x)
  value <- as.character(x)
  if( anyNA(value) || any(!nzchar(value)) ) return(NULL)
  full <- grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", value)
  if( all(full) ){
    out <- suppressWarnings(as.Date(value))
    if( anyNA(out) ) return(NULL)
    return(out)
  }
  month <- grepl("^[0-9]{4}-(0[1-9]|1[0-2])$", value)
  year <- grepl("^[0-9]{4}$", value)
  if( all(month | year) ) return(value)
  NULL
}

#' Read one batch cache entry, treating corruption as a miss (internal)
#' @noRd
.rsCacheRead <- function(path){
  obj <- tryCatch(readRDS(path), error = function(e) NULL)
  if( !is.list(obj) || !is.data.frame(obj$data) ||
      !all(c("source", "date", "value") %in% names(obj$data)) )
    return(NULL)
  if( !is.character(obj$data$source) ||
      length(obj$data$source) != nrow(obj$data) ||
      !is.numeric(obj$data$value) || !is.null(dim(obj$data$value)) ||
      length(obj$data$value) != nrow(obj$data) )
    return(NULL)
  obj$data$date <- .normaliseBatchDates(obj$data$date, nrow(obj$data))
  if( is.null(obj$data$date) ) return(NULL)
  obj
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
#'   by a cache-enabled \code{batch_rs()} call; unmarked custom directories and
#'   broad paths are refused.
#' @returns \code{TRUE}, invisibly.
#' @seealso \code{\link{batch_rs}}, \code{\link{clear_ftw_cache}}
#' @export
clear_rs_cache <- function(cache_dir = NULL){
  if( !is.null(cache_dir) &&
      (!is.character(cache_dir) || length(cache_dir) != 1L ||
       is.na(cache_dir) || !nzchar(cache_dir)) )
    stop("clear_rs_cache(): `cache_dir` must be NULL or one non-empty path.",
         call. = FALSE)
  default_dir <- .rsCacheDir(create = FALSE)
  if( is.null(cache_dir) ) cache_dir <- default_dir
  cache_dir <- .validateCacheClearDir(
    cache_dir, "rs", .RS_CACHE_VERSION, "clear_rs_cache()", default_dir)
  entries <- .cacheEntryFiles(cache_dir, "^[0-9a-f]{24}[.]rds$")
  if( length(entries) > 0L ){
    status <- unlink(entries, recursive = FALSE, force = FALSE)
    if( !identical(status, 0L) )
      warning("clear_rs_cache(): one or more cache entries could not be removed.",
              call. = FALSE)
  }
  return(invisible(TRUE))
}

#' Default batch reducer (internal)
#'
#' Multiple-area time series expose per-area columns plus their union. Batch
#' output has one value per site/date, so the union is its canonical reducer.
#' @noRd
.defaultBatchPipeline <- function(rs){
  out <- as_time_series(boundary_delineation(rs))
  if( !"value" %in% names(out) && "union" %in% names(out) )
    out$value <- out$union
  out
}


# ---- one site, with retries ---------------------------------------------

#' Fetch and reduce a single site, retrying transient failures (internal)
#'
#' Retries retrieval only, relaxing the cloud ceiling to \code{cloud_relax} on
#' the final attempt when more than one attempt was requested. Once retrieval
#' succeeds, the deterministic pipeline is evaluated exactly once. Returns the
#' reduced data.frame plus provenance, or the last condition object if every
#' retrieval attempt failed.
#' @noRd
.batchFetchOne <- function(longitude, latitude, start_date, end_date, index_name,
                           source, x_metres, y_metres, max_cloud_cover, params,
                           pipeline, retries, cloud_relax){
  last <- NULL
  rs <- NULL
  for( i in seq_len(retries) ){
    relax <- retries > 1L && i == retries
    cc <- if( relax ) max(max_cloud_cover, cloud_relax) else max_cloud_cover
    fetched <- tryCatch(
      get_rs(longitude = longitude, latitude = latitude,
             start_date = as.character(start_date),
             end_date = as.character(end_date), index_name = index_name,
             source = source, x_metres = x_metres, y_metres = y_metres,
             max_cloud_cover = cc, params = params),
      error = identity)
    if( !inherits(fetched, "error") ){
      rs <- fetched
      break
    }
    last <- fetched
  }
  if( is.null(rs) ) return(last)

  ## Pipeline errors are deterministic for the fetched object; do not repeat
  ## them as if they were transient endpoint failures.
  d <- tryCatch(pipeline(rs), error = identity)
  if( inherits(d, "error") ) return(d)
  if( !is.data.frame(d) || !all(c("date", "value") %in% names(d)) )
    return(simpleError("the pipeline must return a data.frame with `date` and `value`."))
  dates <- .normaliseBatchDates(d$date, nrow(d))
  if( is.null(dates) )
    return(simpleError("the pipeline returned invalid dates/bucket keys."))
  if( !is.numeric(d$value) || !is.null(dim(d$value)) ||
      length(d$value) != nrow(d) )
    return(simpleError("the pipeline must return `value` as one numeric vector."))
  list(data = data.frame(source = rep(rs$geom$source, nrow(d)), date = dates,
                         value = d$value, stringsAsFactors = FALSE),
       n_scenes = length(rs$scenes), cloud = cc)
}


# ---- the sweep ----------------------------------------------------------

#' Build spectral index time series for multiple sites
#'
#' For each row of \code{sites}, retrieves scenes with \code{get_rs()} and
#' passes the resulting \code{ocular} object to a post-retrieval function. The
#' default function runs field boundary delineation and creates a time series.
#' Per-site results are combined in one long \code{data.frame}. Completed
#' results can be cached on disk, and failed retrievals can be attempted again.
#'
#' @param sites A \code{data.frame} with columns \code{id}, \code{longitude},
#'   \code{latitude}, \code{start_date}, and \code{end_date}; one row per
#'   site and date range.
#' @param index_name Index to calculate; see \code{get_rs()}.
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
#' @param retries Number of retrieval attempts per site. The post-retrieval
#'   function is evaluated once, after retrieval succeeds.
#' @param cloud_relax Cloud-cover ceiling considered on the final attempt when
#'   \code{retries > 1}. The final ceiling is the greater of this value and
#'   \code{max_cloud_cover}.
#' @param on_error \code{"skip"} (default) records the failure and continues;
#'   \code{"stop"} aborts the sweep.
#' @param quiet Whether to suppress per-site progress messages.
#' @returns A \code{data.frame} with \code{id}, \code{index}, \code{source},
#'   \code{date}, and \code{value}. Failed sites are absent and recorded in the
#'   \code{"failures"} attribute with their \code{id}, \code{index}, and
#'   \code{reason}.
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
#' attr(ts, "failures")
#' }
#' @export
batch_rs <- function( sites, index_name = "EVI2", source = "auto",
                      x_metres = 1000L, y_metres = NULL,
                      max_cloud_cover = 50L, params = NULL, pipeline = NULL,
                      cache = TRUE, cache_dir = NULL, pipeline_cache_key = NULL,
                      retries = 3L,
                      cloud_relax = 80L, on_error = c("skip", "stop"),
                      quiet = FALSE ){

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
      !index_name %in% unique(c(names(s2_index_list), names(landsat_index_list))) )
    stop("batch_rs(): `index_name` is not a supported fine-source index.",
         call. = FALSE)
  if( !is.character(source) || length(source) != 1L || is.na(source) )
    stop("batch_rs(): `source` must be auto, sentinel-2, landsat-8, or landsat-5.",
         call. = FALSE)
  source <- match.arg(source,
                      c("auto", "sentinel-2", "landsat-8", "landsat-5"))
  if( !.isScalarNumber(x_metres) || x_metres <= 0 )
    stop("batch_rs(): `x_metres` must be a finite positive numeric.",
         call. = FALSE)
  if( is.null(y_metres) ) y_metres <- x_metres
  if( !.isScalarNumber(y_metres) || y_metres <= 0 )
    stop("batch_rs(): `y_metres` must be NULL or a finite positive numeric.",
         call. = FALSE)
  if( !.isScalarNumber(max_cloud_cover) || max_cloud_cover < 0 ||
      max_cloud_cover > 100 )
    stop("batch_rs(): `max_cloud_cover` must be in [0, 100].", call. = FALSE)
  if( !.isScalarNumber(cloud_relax) || cloud_relax < 0 || cloud_relax > 100 )
    stop("batch_rs(): `cloud_relax` must be in [0, 100].", call. = FALSE)
  if( !.isWholeNumber(retries) || retries < 1L ){
    stop("batch_rs(): `retries` must be a single positive integer.", call. = FALSE)
  }
  retries <- as.integer(retries)
  if( !is.logical(cache) || length(cache) != 1L || is.na(cache) )
    stop("batch_rs(): `cache` must be TRUE or FALSE.", call. = FALSE)
  if( !is.logical(quiet) || length(quiet) != 1L || is.na(quiet) )
    stop("batch_rs(): `quiet` must be TRUE or FALSE.", call. = FALSE)
  if( !is.null(cache_dir) &&
      (!is.character(cache_dir) || length(cache_dir) != 1L ||
       is.na(cache_dir) || !nzchar(cache_dir)) )
    stop("batch_rs(): `cache_dir` must be NULL or one non-empty path.",
         call. = FALSE)
  if( !is.null(pipeline_cache_key) &&
      (!is.character(pipeline_cache_key) || length(pipeline_cache_key) != 1L ||
       is.na(pipeline_cache_key) || !nzchar(pipeline_cache_key)) )
    stop("batch_rs(): `pipeline_cache_key` must be NULL or one non-empty string.",
         call. = FALSE)

  if( is.null(params) ) params <- rs_params()
  else {
    if( !is.list(params) || is.null(names(params)) || anyNA(names(params)) ||
        any(!nzchar(names(params))) || anyDuplicated(names(params)) )
      stop("batch_rs(): `params` must be a uniquely named list.", call. = FALSE)
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
  dir_use <- if( !isTRUE(cache) ) NULL else if( is.null(cache_dir) )
    .rsCacheDir() else .initialiseCacheDir(
      cache_dir, "rs", .RS_CACHE_VERSION, "batch_rs()")

  rows <- list(); fails <- list()
  for( i in seq_len(nrow(sites)) ){
    s_raw <- sites[i, , drop = FALSE]
    checked <- tryCatch(.validateBatchSite(s_raw), error = identity)
    if( inherits(checked, "error") ){
      id <- tryCatch(as.character(s_raw$id[[1L]]), error = function(e) paste0("row ", i))
      msg <- conditionMessage(checked)
      if( identical(on_error, "stop") )
        stop("batch_rs(): ", id, " / ", index_name, ": ", msg,
             call. = FALSE)
      if( !isTRUE(quiet) ) message("  batch_rs: skipped ", id, " -- ", msg)
      fails[[length(fails) + 1L]] <- data.frame(id = id, index = index_name,
                                                reason = msg,
                                                stringsAsFactors = FALSE)
      next
    }
    s <- checked
    key <- .rsCacheKey(s$longitude, s$latitude, s$start_date, s$end_date,
                       index_name, source, x_metres, y_metres, max_cloud_cover,
                       params, retries, cloud_relax, pipeline_cache_key)
    f   <- if( is.null(dir_use) ) NULL else file.path(dir_use, key)

    if( !is.null(f) && file.exists(f) ){
      hit <- .rsCacheRead(f)
      if( !is.null(hit) ){
        if( !isTRUE(quiet) )
          message("  batch_rs: ", s$id, " / ", index_name, " (cached)")
        if( nrow(hit$data) > 0L )
          rows[[length(rows) + 1L]] <- cbind(
            id = rep(s$id, nrow(hit$data)),
            index = rep(index_name, nrow(hit$data)), hit$data,
            stringsAsFactors = FALSE)
        next
      }
      if( !isTRUE(quiet) )
        message("  batch_rs: ignoring unreadable cache entry for ", s$id)
    }

    got <- tryCatch(
      .batchFetchOne(s$longitude, s$latitude, s$start_date, s$end_date,
                     index_name, source, x_metres, y_metres,
                     max_cloud_cover, params, pipeline, retries,
                     cloud_relax),
      error = identity)
    if( inherits(got, "error") ){
      msg <- conditionMessage(got)
      if( identical(on_error, "stop") ){
        stop("batch_rs(): ", s$id, " / ", index_name, ": ", msg, call. = FALSE)
      }
      if( !isTRUE(quiet) ) message("  batch_rs: skipped ", s$id, " -- ", msg)
      fails[[ length(fails) + 1L ]] <- data.frame(id = s$id, index = index_name,
                                                  reason = msg,
                                                  stringsAsFactors = FALSE)
      next
    }
    if( !is.null(f) ){
      cache_obj <- list(data = got$data,
                        meta = list(key = key, created_at = Sys.time(),
                                    cloud_used = got$cloud))
      tryCatch(saveRDS(cache_obj, f),
               error = function(e)
                 warning("batch_rs(): cache write failed for ", s$id, ": ",
                         conditionMessage(e), call. = FALSE))
    }
    if( !isTRUE(quiet) ){
      message("  batch_rs: ", s$id, " / ", index_name, " -- ", got$n_scenes,
              " scenes (", nrow(got$data), " dates)")
    }
    if( nrow(got$data) > 0L )
      rows[[ length(rows) + 1L ]] <- cbind(
        id = rep(s$id, nrow(got$data)),
        index = rep(index_name, nrow(got$data)), got$data,
        stringsAsFactors = FALSE)
  }

  out <- if( length(rows) > 0L ) do.call(rbind, rows) else
    data.frame(id = character(0L), index = character(0L), source = character(0L),
               date = as.Date(character(0L)), value = numeric(0L),
               stringsAsFactors = FALSE)
  rownames(out) <- NULL
  attr(out, "failures") <- if( length(fails) > 0L ) do.call(rbind, fails) else NULL
  return(out)
}
