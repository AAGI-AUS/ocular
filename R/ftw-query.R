# =========================================================================
# FTW query and cache
# =========================================================================
# Query field polygons from GeoParquet and cache them for add_ftw_prior().
# Default column names follow fiboa. The caller supplies the dataset's CRS.
# DuckDB extensions must already be installed: spatial reads geometry,
# httpfs reads remote files, and icu extracts source years from timestamps.

.FTW_CACHE_VERSION <- "v2"

#' Identify a field dataset without exposing credentials or local paths
#' @noRd
.safeFtwSourceLabel <- function(source){

  no_query <- sub("[?#].*$", "", source)
  if( grepl("^[A-Za-z][A-Za-z0-9+.-]*://", no_query) ){
    scheme <- tolower(sub(":.*$", "", no_query))
    path_part <- sub("^[A-Za-z][A-Za-z0-9+.-]*://[^/]*", "", no_query)
    leaf <- basename(sub("/+$", "", path_part))
    if( !nzchar(leaf) || leaf %in% c(".", "/") ){
      leaf <- "source"
    }
    return(paste0(scheme, "://.../", leaf))
  }
  leaf <- basename(no_query)
  if( !nzchar(leaf) || leaf %in% c(".", "/") ) { "local source" } else { leaf }
}

#' Redact a source string from a caught condition
#' @noRd
.redactFtwCondition <- function(condition, source){

  msg <- conditionMessage(condition)
  # DuckDB can include a truncated query after "LINE n:", beyond the source
  # string that the replacement below can match.
  msg <- gsub("\nLINE [0-9]+:[^\n]*\n[ ^]*", "", msg, perl = TRUE)
  if( nzchar(source) ){
    msg <- gsub(source, .safeFtwSourceLabel(source), msg, fixed = TRUE)
  }
  return(msg)
}


# =========================================================================
# FTW cache helpers
# =========================================================================

#' Versioned field polygon cache for the current user
#' @noRd
.ftwCacheDir <- function(create = TRUE){

  base <- tools::R_user_dir("ocular", which = "cache")
  d <- file.path(base, "ftw", .FTW_CACHE_VERSION)
  if( isTRUE(create) ){
    d <- .initialiseCacheDir(d, "ftw", .FTW_CACHE_VERSION, ".ftwCacheDir()")
  }
  return(d)
}

#' Round a query box outward to a shared cache grid
#'
#' Nearby requests can share an entry. Query the enlarged box, then select
#' the field containing the point.
#' @noRd
.ftwQuantiseBbox <- function(bbox_wgs84, grid_deg = 0.01){

  required <- c("xmin", "ymin", "xmax", "ymax")
  if( !is.numeric(bbox_wgs84) || !all(required %in% names(bbox_wgs84)) ||
      anyNA(bbox_wgs84[required]) ||
      any(!is.finite(bbox_wgs84[required])) ||
      bbox_wgs84[["xmin"]] >= bbox_wgs84[["xmax"]] ||
      bbox_wgs84[["ymin"]] >= bbox_wgs84[["ymax"]] ||
      bbox_wgs84[["xmin"]] < -180 || bbox_wgs84[["xmax"]] > 180 ||
      bbox_wgs84[["ymin"]] < -90 || bbox_wgs84[["ymax"]] > 90 ){
    stop("bbox_wgs84 must be a finite named c(xmin, ymin, xmax, ymax) ",
         "with increasing bounds.", call. = FALSE)
  }
  if( !.isScalarNumber(grid_deg) || grid_deg <= 0 ){
    stop("grid_deg must be a finite positive numeric scalar.", call. = FALSE)
  }
  xmin <- floor(bbox_wgs84[["xmin"]] / grid_deg) * grid_deg
  ymin <- floor(bbox_wgs84[["ymin"]] / grid_deg) * grid_deg
  xmax <- ceiling(bbox_wgs84[["xmax"]] / grid_deg) * grid_deg
  ymax <- ceiling(bbox_wgs84[["ymax"]] / grid_deg) * grid_deg
  return(c(xmin = unname(xmin), ymin = unname(ymin),
           xmax = unname(xmax), ymax = unname(ymax)))
}

#' Hash an R value for a cache key
#' @noRd
.ftwStableHash <- function(payload){

  return(.stableHash(payload))
}

#' Cache key for a field dataset, bounding box and query settings
#' @noRd
.ftwCacheKey <- function(source_id, bbox_quantised, query = list()){

  parts <- list(source_id = source_id,
                bbox_quantised = round(bbox_quantised, 6L),
                query = query)
  return(substr(.ftwStableHash(parts), 1L, 16L))
}

#' Identify a dataset for caching
#'
#' For a local file, include its full path, size and modification time so a
#' replacement file gets a new entry. Use a URL or glob string as supplied.
#' @noRd
.ftwSourceId <- function(source){

  is_url <- grepl("^[A-Za-z][A-Za-z0-9+.-]*://", source)
  is_glob <- grepl("[*?[]", source) && !file.exists(source)
  if( is_url || is_glob ){
    return(source)
  }
  info <- file.info(source, extra_cols = FALSE)
  return(list(path = normalizePath(source, mustWork = FALSE),
              size = as.numeric(info$size),
              mtime = as.numeric(info$mtime)))
}


# =========================================================================
# Polygon cache
# =========================================================================

#' Read cached polygons and update their access time
#'
#' Return \code{list(polygons, meta)} on a hit, or \code{NULL} on a miss.
#' @noRd
.ftwCacheRead <- function(key, cache_dir = .ftwCacheDir()){

  f <- file.path(cache_dir, paste0(key, ".rds"))
  if( !file.exists(f) ){
    return(NULL)
  }
  obj <- tryCatch(readRDS(f), error = function(e){
    return(NULL)
  })
  if( !is.list(obj) || !inherits(obj$polygons, "sf") ){
    return(NULL)
  }
  if( !is.list(obj$meta) ){
    obj$meta <- list()
  }
  obj$meta$last_accessed_at <- Sys.time()
  tryCatch(saveRDS(obj, f), error = function(e){
    return(NULL)
  })
  return(obj)
}

#' Cache polygons with their metadata
#' @noRd
.ftwCacheWrite <- function(key, polygons, meta = list(),
                           cache_dir = .ftwCacheDir()){

  f <- file.path(cache_dir, paste0(key, ".rds"))
  meta$created_at <- Sys.time()
  meta$last_accessed_at <- meta$created_at
  tryCatch(saveRDS(list(polygons = polygons, meta = meta), f),
           error = function(e){
             return(warning("FTW cache write failed: ",
                            conditionMessage(e), call. = FALSE))
           })
  return(invisible(f))
}

#' @noRd
.ftwCacheList <- function(cache_dir = .ftwCacheDir()){

  return(.cacheEntryFiles(cache_dir, "^[0-9a-f]{16}[.]rds$"))
}

#' Remove the least recently accessed entries above the cache limit
#' @noRd
.ftwCacheEvictLRU <- function(cache_dir = .ftwCacheDir(), max_entries = 200L){

  if( !.isWholeNumber(max_entries) || max_entries < 1L ){
    stop("max_entries must be a positive integer.", call. = FALSE)
  }
  max_entries <- as.integer(max_entries)
  files <- .ftwCacheList(cache_dir)
  if( length(files) <= max_entries ){
    return(invisible(0L))
  }
  atime <- vapply(files, function(f){
    obj <- tryCatch(readRDS(f), error = function(e){
      return(NULL)
    })
    stamp <- if( is.list(obj) && is.list(obj$meta) ) { obj$meta$last_accessed_at } else { NULL }
    value <- suppressWarnings(as.numeric(stamp %||% file.info(f)$mtime))
    if( length(value) != 1L || !is.finite(value) ){
      value <- suppressWarnings(as.numeric(file.info(f)$mtime))
    }
    return(value)
  }, numeric(1L))
  drop <- files[order(atime)][seq_len(length(files) - max_entries)]
  unlink(drop)
  return(invisible(length(drop)))
}

#' Clear the FTW on-disk cache
#'
#' Removes only ocular FTW cache-entry files from its versioned per-user
#' directory. It never recursively removes the directory or unrelated files.
#' @examples
#' \dontrun{
#'   clear_ftw_cache()   # remove cached field boundary query results
#' }
#' @returns \code{TRUE}, invisibly.
#' @export
clear_ftw_cache <- function(){

  default_dir <- .ftwCacheDir(create = FALSE)
  cache_dir <- .validateCacheClearDir(
    default_dir, "ftw", .FTW_CACHE_VERSION, "clear_ftw_cache()", default_dir)
  entries <- .cacheEntryFiles(cache_dir, "^[0-9a-f]{16}[.]rds$")
  if( length(entries) > 0L ){
    status <- unlink(entries, recursive = FALSE, force = FALSE)
    if( !identical(status, 0L) ){
      warning("clear_ftw_cache(): one or more cache entries could not be removed.",
              call. = FALSE)
    }
  }
  return(invisible(TRUE))
}


# =========================================================================
# Field polygon query
# =========================================================================

#' Read field polygons whose bounds overlap a query box
#'
#' Query GeoParquet through DuckDB spatial. Keep the dataset's CRS and derive
#' source years from the configured date column when available.
#'
#' @param bbox_wgs84 Bounds in EPSG:4326, ordered xmin, ymin, xmax, ymax.
#' @param source Local path, HTTPS URL or S3 key for the GeoParquet dataset.
#' @param crs Dataset CRS, as an EPSG code or \code{sf} CRS. Transform the query
#'   box to this CRS before filtering.
#' @param geometry_col,datetime_col Geometry and date columns, using fiboa
#'   names by default.
#' @returns An \code{sf} polygon layer with \code{source_year}, or \code{NULL}
#'   when the query returns no candidate polygons.
#' @noRd
.ftwQuery <- function(bbox_wgs84, source, crs = 4326L,
                      geometry_col = "geometry",
                      datetime_col = "determination_datetime"){

  for( pkg in c("duckdb", "DBI", "sf") ){
    if( !requireNamespace(pkg, quietly = TRUE) ){
      stop(".ftwQuery() needs the '", pkg, "' package.", call. = FALSE)
    }
  }

  .ftwQuantiseBbox(bbox_wgs84, grid_deg = 1)
  if( !is.character(source) || length(source) != 1L || is.na(source) ||
      !nzchar(trimws(source)) ){
    stop(".ftwQuery(): source must be one non-empty path or URL.",
         call. = FALSE)
  }
  if( !is.character(geometry_col) || length(geometry_col) != 1L ||
      is.na(geometry_col) || !nzchar(geometry_col) ){
    stop(".ftwQuery(): geometry_col must be one non-empty column name.",
         call. = FALSE)
  }
  if( !is.null(datetime_col) &&
      (!is.character(datetime_col) || length(datetime_col) != 1L ||
       is.na(datetime_col) || !nzchar(datetime_col)) ){
    stop(".ftwQuery(): datetime_col must be NULL or one non-empty column name.",
         call. = FALSE)
  }

  # Filter in the dataset's CRS. Delineation later projects the selected
  # polygon to the feature grid.
  out_crs <- tryCatch(sf::st_crs(crs), error = function(e){
    return(sf::st_crs(NA))
  })
  if( is.na(out_crs) ){
    stop(".ftwQuery(): `crs` is missing or invalid. Supply the source ",
         "dataset CRS explicitly.", call. = FALSE)
  }
  flt <- bbox_wgs84
  if( out_crs != sf::st_crs(4326) ){
    bb <- sf::st_bbox(c(xmin = bbox_wgs84[["xmin"]], ymin = bbox_wgs84[["ymin"]],
                        xmax = bbox_wgs84[["xmax"]], ymax = bbox_wgs84[["ymax"]]),
                      crs = 4326)
    bbT <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(bb), out_crs))
    flt <- c(xmin = bbT[["xmin"]], ymin = bbT[["ymin"]],
             xmax = bbT[["xmax"]], ymax = bbT[["ymax"]])
  }

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # Disable automatic downloads before loading installed extensions.
  for( setting in c("autoinstall_known_extensions",
                    "autoload_known_extensions") ){
    tryCatch(
      DBI::dbExecute(con, paste("SET", setting, "= false")),
      error = function(e){
        stop(".ftwQuery(): could not disable DuckDB extension auto-install ",
             "(`", setting, "`): ", conditionMessage(e), call. = FALSE)
      })
  }
  tryCatch(
    DBI::dbExecute(con, "LOAD spatial"),
    error = function(e){
      stop(".ftwQuery(): DuckDB's spatial extension is not available locally. ",
           "ocular will not download it automatically. ", conditionMessage(e),
           call. = FALSE)
    })
  have_icu <- !is.null(datetime_col) && tryCatch({
    DBI::dbExecute(con, "LOAD icu")
    TRUE
  }, error = function(e){
    return(FALSE)
  })
  if( grepl("^(https?|s3)://", source, ignore.case = TRUE) ){
    tryCatch(
      DBI::dbExecute(con, "LOAD httpfs"),
      error = function(e){
        stop(".ftwQuery(): a remote source needs DuckDB's httpfs extension, ",
             "which is not available locally. ocular will not download it ",
             "automatically. ", conditionMessage(e), call. = FALSE)
      })
  }

  # DuckDB cannot bind read_parquet()'s path argument. Quote it as a SQL
  # literal, while keeping the numeric bounds as bound parameters.
  src_lit <- DBI::dbQuoteString(con, source)
  # The bounding box filter finds candidates. The caller selects the field
  # containing the point.
  geom_id <- DBI::dbQuoteIdentifier(con, geometry_col)
  dt_sel <- if( have_icu && !is.null(datetime_col) ){
    paste0(", EXTRACT(epoch FROM ",
           DBI::dbQuoteIdentifier(con, datetime_col), ") AS dt_epoch")
  }else{
    ""
  }
  sql <- sprintf(
    "SELECT ST_AsWKB(%1$s) AS wkb%2$s
       FROM read_parquet(%3$s)
      WHERE ST_XMax(%1$s) >= ? AND ST_XMin(%1$s) <= ?
        AND ST_YMax(%1$s) >= ? AND ST_YMin(%1$s) <= ?",
    geom_id, dt_sel, src_lit)

  res <- tryCatch(
    DBI::dbGetQuery(con, sql, params = list(
      flt[["xmin"]], flt[["xmax"]],
      flt[["ymin"]], flt[["ymax"]])),
    error = function(e){
      stop(".ftwQuery() failed for ", .safeFtwSourceLabel(source), ": ",
           .redactFtwCondition(e, source),
           "\nCheck `geometry_col`, `datetime_col`, and `crs` against the ",
           "source schema.", call. = FALSE)
    })

  if( is.null(res) || nrow(res) == 0L ){
    return(NULL)
  }

  wkb <- res$wkb
  if( !is.list(wkb) ){
    wkb <- as.list(wkb)
  }
  geom <- sf::st_as_sfc(structure(lapply(wkb, as.raw), class = "WKB"),
                        EWKB = FALSE)
  sf::st_crs(geom) <- out_crs
  yr <- if( have_icu && !is.null(datetime_col) && !is.null(res$dt_epoch) ){
    suppressWarnings(as.integer(format(
      as.POSIXct(res$dt_epoch, origin = "1970-01-01", tz = "UTC"), "%Y")))
  }else{
    rep(NA_integer_, nrow(res))
  }
  if( !is.null(datetime_col) && !isTRUE(have_icu) ){
    warning("FTW: DuckDB's icu extension is unavailable locally. ",
            "The source_year column was set to NA.", call. = FALSE)
  }
  return(sf::st_sf(source_year = yr, geometry = geom))
}


# =========================================================================
# Attach a field boundary prior
# =========================================================================

#' Attach an optional field boundary prior to an ocular object
#'
#' Queries a user-supplied field boundary GeoParquet, then stores the polygon
#' containing the supplied point and its representative \code{source_year} on
#' the \code{ocular} object. The default column names follow the fiboa schema.
#' Fields of The World (FTW) is one field boundary dataset that can be supplied
#' in a compatible GeoParquet. fiboa and FTW are not interchangeable names.
#'
#' The selected polygon can provide calibration support for
#' \code{segment_area()} and \code{boundary_delineation()}. When a source year
#' is available, it can also act as an age-weighted soft prior during the
#' initial area search. Use \code{diagnose_against_ftw()} to compare the
#' resulting mask with a field boundary reference.
#'
#' This feature is experimental. A polygon used to guide delineation cannot
#' also provide independent evidence of its accuracy. Assess improvement
#' against a separate reference.
#'
#' Queries require the optional \code{duckdb} and \code{DBI} packages and a
#' locally installed DuckDB spatial extension. Remote sources also require the
#' DuckDB \code{httpfs} extension. The DuckDB \code{icu} extension is needed
#' for \code{source_year}, without which the soft prior is not applied.
#' \code{ocular} does not download or install DuckDB extensions
#' automatically.
#'
#' @param rs An \code{ocular} object returned by \code{get_rs()}.
#' @param source Path or URL to a compatible field boundary GeoParquet covering
#'   the area of interest. Local files, HTTPS objects, and S3 keys are accepted.
#' @param crs The source dataset's geometry CRS, supplied as an EPSG code or
#'   \code{sf} CRS. Set it when the source is not EPSG:4326.
#' @param refresh If \code{TRUE}, bypass the cache and repeat the query.
#' @param max_cache_entries Maximum number of entries retained in the on-disk
#'   cache. The default is 200. Least-recently accessed entries are removed
#'   when the limit is exceeded.
#' @param geometry_col,datetime_col Geometry and date column names. Their
#'   defaults follow the fiboa schema. Override them for other compatible
#'   exports, or set \code{datetime_col = NULL} when no observation date is
#'   available.
#' @returns The \code{ocular} object with the selected field boundary attached,
#'   or the unchanged object with a message when no field contains the supplied
#'   point.
#' @examples
#' \dontrun{
#'   rs <- get_rs(171.80, -43.90, "2021-10-01", "2022-04-30", x_metres = 2000)
#'   rs <- add_ftw_prior(rs,
#'     source = "/path/to/field-boundaries.parquet",
#'     crs = 4326)
#'   agreement <- diagnose_against_ftw(boundary_delineation(rs))
#' }
#' @export
add_ftw_prior <- function(rs, source, crs = 4326L, refresh = FALSE,
                          max_cache_entries = 200L,
                          geometry_col = "geometry",
                          datetime_col = "determination_datetime"){

  if( !is_rs(rs) ){
    stop("add_ftw_prior(): rs must be an ocular object.", call. = FALSE)
  }
  if( missing(source) || !is.character(source) || length(source) != 1L ||
      is.na(source) || !nzchar(trimws(source)) ){
    stop("add_ftw_prior(): `source` must be a single path or URL to a field ",
         "boundary GeoParquet.", call. = FALSE)
  }
  if( !is.logical(refresh) || length(refresh) != 1L || is.na(refresh) ){
    stop("add_ftw_prior(): `refresh` must be TRUE or FALSE.", call. = FALSE)
  }
  if( !.isWholeNumber(max_cache_entries) || max_cache_entries < 1L ){
    stop("add_ftw_prior(): `max_cache_entries` must be a positive integer.",
         call. = FALSE)
  }
  max_cache_entries <- as.integer(max_cache_entries)

  source_crs <- tryCatch(sf::st_crs(crs), error = function(e){
    return(sf::st_crs(NA))
  })
  if( is.na(source_crs) ){
    stop("add_ftw_prior(): `crs` is missing or invalid.", call. = FALSE)
  }
  if( !is.character(geometry_col) || length(geometry_col) != 1L ||
      is.na(geometry_col) || !nzchar(geometry_col) ){
    stop("add_ftw_prior(): `geometry_col` must be one non-empty name.",
         call. = FALSE)
  }
  if( !is.null(datetime_col) &&
      (!is.character(datetime_col) || length(datetime_col) != 1L ||
       is.na(datetime_col) || !nzchar(datetime_col)) ){
    stop("add_ftw_prior(): `datetime_col` must be NULL or one non-empty name.",
         call. = FALSE)
  }

  bbox <- rs$geom$bbox
  if( is.null(bbox) ){
    stop("add_ftw_prior(): rs$geom$bbox is missing.", call. = FALSE)
  }
  bbox <- c(xmin = unname(bbox[["xmin"]]), ymin = unname(bbox[["ymin"]]),
            xmax = unname(bbox[["xmax"]]), ymax = unname(bbox[["ymax"]]))

  qb <- .ftwQuantiseBbox(bbox)
  query <- list(crs_wkt = source_crs$wkt,
                geometry_col = geometry_col,
                datetime_col = datetime_col)
  key <- .ftwCacheKey(.ftwSourceId(source), qb, query = query)

  polygons <- NULL
  if( !refresh ){
    hit <- .ftwCacheRead(key)
    if( !is.null(hit) ){
      polygons <- hit$polygons
    }
  }
  if( is.null(polygons) ){
    polygons <- .ftwQuery(bbox_wgs84 = qb, source = source, crs = source_crs,
                          geometry_col = geometry_col,
                          datetime_col = datetime_col)
    # Do not cache a query with no usable requested source year, so a later
    # call can try again.
    no_year <- !is.null(datetime_col) && all(is.na(polygons$source_year))
    if( !is.null(polygons) && !no_year ){
      .ftwCacheWrite(key, polygons,
                     meta = list(
                       source_id = .safeFtwSourceLabel(source),
                       source_hash = substr(.stableHash(source), 1L, 16L),
                       bbox = qb, query = query))
      .ftwCacheEvictLRU(max_entries = max_cache_entries)
    }
  }

  if( is.null(polygons) || nrow(polygons) == 0L ){
    message("add_ftw_prior(): no field boundaries intersect the analysis ",
            "area. The ocular object is unchanged.")
    return(rs)
  }

  polygons <- .selectFtwField(polygons, rs)
  if( is.null(polygons) || nrow(polygons) == 0L ){
    message("add_ftw_prior(): no field boundary contains the supplied point. ",
            "The ocular object is unchanged.")
    return(rs)
  }

  yr <- if( !is.null(polygons$source_year) ){
    suppressWarnings(stats::median(polygons$source_year, na.rm = TRUE))
  }else{
    NA_real_
  }
  yr <- if( is.finite(yr) ) { as.integer(round(yr)) } else { NULL }

  rs$geom$ftw_prior <- list(polygon = polygons, source_year = yr)
  # Replace seed support with the newly attached field's support.
  if( !is.null(rs$internals$feature_stack) ){
    updated <- tryCatch(.calibrateScene(rs), error = identity)
    if( !inherits(updated, "error") && is.null(updated$error) ){
      rs$internals$calibration <- updated
    }else{
      msg <- if( inherits(updated, "error") ) { conditionMessage(updated) } else { updated$error }
      warning("add_ftw_prior(): FTW calibration was unavailable. The prior ",
              "calibration was preserved. ", msg, call. = FALSE)
    }
  }
  message(sprintf("add_ftw_prior(): seed field attached%s.",
                  if( is.null(yr) ) { "" } else { sprintf(" (source_year %d)", yr) }))
  return(rs)
}
