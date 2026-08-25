# =========================================================================
# FTW query and cache
# =========================================================================
# Auto-populates rs$geom$ftw_prior$polygon from a field boundary GeoParquet.
# The default columns follow the fiboa schema; Fields of The World is one
# compatible dataset. A DuckDB-spatial bbox query is cached to disk and feeds
# the calibration-source dispatch in .calibrateScene
# (rs$geom$ftw_prior$polygon present -> .ftwSupport).
#
# Default fiboa columns used by the query:
#   id                     unique field identifier
#   determination_datetime last timestamp the boundary was observed  -> source_year
#   area                   field area (hectares)
#   geometry               field polygon (WKB). CRS varies by dataset: the FTW
#                          benchmark is EPSG:4326, but source.coop sets keep
#                          their source CRS (the NZ irrigated set is EPSG:2193).
#                          Pass crs= to add_ftw_prior/.ftwQuery when not 4326.
#
# Queries require duckdb, DBI, sf, and locally available DuckDB extensions.
# Extension auto-installation is disabled, so ocular does not download an
# extension when a query is run. Remote reads require httpfs; S3 sources also
# require the user's connection settings. Timestamp extraction requires icu;
# without it, source_year is unavailable. The remote query path is not covered
# by the offline test suite and requires live interoperability testing.
#
# The optional prior currently uses reference age but no external confidence
# layer because no confidence product is configured by the package.
# =========================================================================

.FTW_CACHE_VERSION <- "v2"

#' Redacted identifier for an FTW data source (internal)
#'
#' Keeps error messages and cache metadata useful without retaining URL query
#' credentials, URL userinfo, or a full local path.
#' @noRd
.safeFtwSourceLabel <- function(source){
  no_query <- sub("[?#].*$", "", source)
  if( grepl("^[A-Za-z][A-Za-z0-9+.-]*://", no_query) ){
    scheme <- tolower(sub(":.*$", "", no_query))
    path_part <- sub("^[A-Za-z][A-Za-z0-9+.-]*://[^/]*", "", no_query)
    leaf <- basename(sub("/+$", "", path_part))
    if( !nzchar(leaf) || leaf %in% c(".", "/") ) leaf <- "source"
    return(paste0(scheme, "://.../", leaf))
  }
  leaf <- basename(no_query)
  if( !nzchar(leaf) || leaf %in% c(".", "/") ) "local source" else leaf
}

#' Redact a source string from a caught condition (internal)
#' @noRd
.redactFtwCondition <- function(condition, source){
  msg <- conditionMessage(condition)
  if( nzchar(source) )
    msg <- gsub(source, .safeFtwSourceLabel(source), msg, fixed = TRUE)
  msg
}


# ---- FTW cache primitives ---------------------------------------------------

#' On-disk FTW cache directory (per-user, versioned)
#' @noRd
.ftwCacheDir <- function(create = TRUE) {
  base <- tools::R_user_dir("ocular", which = "cache")
  d <- file.path(base, "ftw", .FTW_CACHE_VERSION)
  if (isTRUE(create))
    d <- .initialiseCacheDir(d, "ftw", .FTW_CACHE_VERSION, ".ftwCacheDir()")
  d
}

#' Quantise a bbox to a fixed 0.01-degree grid for cache-key purposes.
#' Two AOIs whose bboxes round to the same cell share a cache entry; the fetch
#' pulls the quantised (slightly enlarged) bbox and consumers clip at use time.
#' @noRd
.ftwQuantiseBbox <- function(bbox_wgs84, grid_deg = 0.01) {
  required <- c("xmin", "ymin", "xmax", "ymax")
  if (!is.numeric(bbox_wgs84) || !all(required %in% names(bbox_wgs84)) ||
      anyNA(bbox_wgs84[required]) ||
      any(!is.finite(bbox_wgs84[required])) ||
      bbox_wgs84[["xmin"]] >= bbox_wgs84[["xmax"]] ||
      bbox_wgs84[["ymin"]] >= bbox_wgs84[["ymax"]] ||
      bbox_wgs84[["xmin"]] < -180 || bbox_wgs84[["xmax"]] > 180 ||
      bbox_wgs84[["ymin"]] < -90 || bbox_wgs84[["ymax"]] > 90)
    stop("bbox_wgs84 must be a finite named c(xmin, ymin, xmax, ymax) ",
         "with increasing bounds.", call. = FALSE)
  if (!.isScalarNumber(grid_deg) || grid_deg <= 0)
    stop("grid_deg must be a finite positive numeric scalar.", call. = FALSE)
  xmin <- floor(bbox_wgs84[["xmin"]] / grid_deg) * grid_deg
  ymin <- floor(bbox_wgs84[["ymin"]] / grid_deg) * grid_deg
  xmax <- ceiling(bbox_wgs84[["xmax"]] / grid_deg) * grid_deg
  ymax <- ceiling(bbox_wgs84[["ymax"]] / grid_deg) * grid_deg
  c(xmin = unname(xmin), ymin = unname(ymin),
    xmax = unname(xmax), ymax = unname(ymax))
}

#' Stable hexadecimal digest of an R payload (base R).
#' @noRd
.ftwStableHash <- function(payload) {
  .stableHash(payload)
}

#' Content-fingerprint cache key from source, bbox, and query options.
#' @noRd
.ftwCacheKey <- function(source_id, bbox_quantised, query = list()) {
  parts <- list(source_id      = source_id,
                bbox_quantised = round(bbox_quantised, 6L),
                query          = query)
  substr(.ftwStableHash(parts), 1L, 16L)
}


# ---- sf-adapted cache read/write (polygons, not rasters) --------------------
# Cache entries store sf polygons in RDS payloads and record last access for
# least-recently-used eviction.

#' Read a cached polygon set. Returns list(polygons, meta) on hit, NULL on miss.
#' Touches last_accessed_at for LRU bookkeeping.
#' @noRd
.ftwCacheRead <- function(key, cache_dir = .ftwCacheDir()) {
  f <- file.path(cache_dir, paste0(key, ".rds"))
  if (!file.exists(f)) return(NULL)
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (!is.list(obj) || !inherits(obj$polygons, "sf")) return(NULL)
  if (!is.list(obj$meta)) obj$meta <- list()
  obj$meta$last_accessed_at <- Sys.time()
  tryCatch(saveRDS(obj, f), error = function(e) NULL)
  obj
}

#' Write a polygon set + sidecar meta to the cache.
#' @noRd
.ftwCacheWrite <- function(key, polygons, meta = list(),
                           cache_dir = .ftwCacheDir()) {
  f <- file.path(cache_dir, paste0(key, ".rds"))
  meta$created_at       <- Sys.time()
  meta$last_accessed_at <- meta$created_at
  tryCatch(saveRDS(list(polygons = polygons, meta = meta), f),
           error = function(e) warning("FTW cache write failed: ",
                                       conditionMessage(e), call. = FALSE))
  invisible(f)
}

#' @noRd
.ftwCacheList <- function(cache_dir = .ftwCacheDir()) {
  .cacheEntryFiles(cache_dir, "^[0-9a-f]{16}[.]rds$")
}

#' Evict least-recently-accessed entries down to max_entries.
#' @noRd
.ftwCacheEvictLRU <- function(cache_dir = .ftwCacheDir(), max_entries = 200L) {
  if (!.isWholeNumber(max_entries) || max_entries < 1L)
    stop("max_entries must be a positive integer.", call. = FALSE)
  max_entries <- as.integer(max_entries)
  files <- .ftwCacheList(cache_dir)
  if (length(files) <= max_entries) return(invisible(0L))
  atime <- vapply(files, function(f) {
    obj <- tryCatch(readRDS(f), error = function(e) NULL)
    stamp <- if( is.list(obj) && is.list(obj$meta) )
      obj$meta$last_accessed_at else NULL
    value <- suppressWarnings(as.numeric(stamp %||% file.info(f)$mtime))
    if( length(value) != 1L || !is.finite(value) )
      value <- suppressWarnings(as.numeric(file.info(f)$mtime))
    value
  }, numeric(1L))
  drop <- files[order(atime)][seq_len(length(files) - max_entries)]
  unlink(drop)
  invisible(length(drop))
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
clear_ftw_cache <- function() {
  default_dir <- .ftwCacheDir(create = FALSE)
  cache_dir <- .validateCacheClearDir(
    default_dir, "ftw", .FTW_CACHE_VERSION, "clear_ftw_cache()", default_dir)
  entries <- .cacheEntryFiles(cache_dir, "^[0-9a-f]{16}[.]rds$")
  if( length(entries) > 0L ){
    status <- unlink(entries, recursive = FALSE, force = FALSE)
    if( !identical(status, 0L) )
      warning("clear_ftw_cache(): one or more cache entries could not be removed.",
              call. = FALSE)
  }
  invisible(TRUE)
}


# ---- the DuckDB-spatial query ------------------------------------------------

#' Query a field boundary GeoParquet within a bounding box (internal)
#'
#' Uses DuckDB spatial to read polygons intersecting a bounding box from a
#' compatible GeoParquet. The default geometry and date columns follow the
#' fiboa schema. Returns an \code{sf} polygon layer in the source dataset's CRS,
#' with a \code{source_year} column derived from the configured date column, or
#' \code{NULL} when nothing intersects.
#'
#' @param bbox_wgs84 Named numeric vector \code{c(xmin, ymin, xmax, ymax)} in
#'   EPSG:4326.
#' @param source Path or URL to a compatible field boundary GeoParquet: a local
#'   file, an HTTPS object, or an S3 key.
#' @param crs The source dataset's geometry CRS, supplied as an EPSG code or
#'   \code{sf} CRS. When it is not EPSG:4326, the bounding box is reprojected
#'   before filtering and the returned polygons retain the source CRS.
#' @param geometry_col,datetime_col Geometry and date column names. Their
#'   defaults follow the fiboa core schema; override them for other compatible
#'   exports.
#' @noRd
.ftwQuery <- function(bbox_wgs84, source, crs = 4326L,
                      geometry_col = "geometry",
                      datetime_col = "determination_datetime") {
  for (pkg in c("duckdb", "DBI", "sf"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(".ftwQuery() needs the '", pkg, "' package.", call. = FALSE)

  .ftwQuantiseBbox(bbox_wgs84, grid_deg = 1)
  if (!is.character(source) || length(source) != 1L || is.na(source) ||
      !nzchar(trimws(source)))
    stop(".ftwQuery(): source must be one non-empty path or URL.",
         call. = FALSE)
  if (!is.character(geometry_col) || length(geometry_col) != 1L ||
      is.na(geometry_col) || !nzchar(geometry_col))
    stop(".ftwQuery(): geometry_col must be one non-empty column name.",
         call. = FALSE)
  if (!is.null(datetime_col) &&
      (!is.character(datetime_col) || length(datetime_col) != 1L ||
       is.na(datetime_col) || !nzchar(datetime_col)))
    stop(".ftwQuery(): datetime_col must be NULL or one non-empty column name.",
         call. = FALSE)

  # Field boundary datasets keep their source CRS while the query bbox is
  # WGS84. Reproject the bbox when necessary; .ftwPolygonMask later reprojects
  # the selected geometry onto the feature grid.
  out_crs <- tryCatch(sf::st_crs(crs), error = function(e) sf::st_crs(NA))
  if (is.na(out_crs))
    stop(".ftwQuery(): `crs` is missing or invalid; supply the source ",
         "dataset CRS explicitly.", call. = FALSE)
  flt <- bbox_wgs84
  if (out_crs != sf::st_crs(4326)) {
    bb  <- sf::st_bbox(c(xmin = bbox_wgs84[["xmin"]], ymin = bbox_wgs84[["ymin"]],
                         xmax = bbox_wgs84[["xmax"]], ymax = bbox_wgs84[["ymax"]]),
                       crs = 4326)
    bbT <- sf::st_bbox(sf::st_transform(sf::st_as_sfc(bb), out_crs))
    flt <- c(xmin = bbT[["xmin"]], ymin = bbT[["ymin"]],
             xmax = bbT[["xmax"]], ymax = bbT[["ymax"]])
  }

  con <- DBI::dbConnect(duckdb::duckdb())
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # The spatial extension reads GeoParquet geometry, icu supports timestamp
  # extraction from determination_datetime, and httpfs supports remote reads.
  # Disable extension auto-installation before loading any extension.
  for (setting in c("autoinstall_known_extensions",
                    "autoload_known_extensions")) {
    tryCatch(
      DBI::dbExecute(con, paste("SET", setting, "= false")),
      error = function(e)
        stop(".ftwQuery(): could not disable DuckDB extension auto-install ",
             "(`", setting, "`): ", conditionMessage(e), call. = FALSE))
  }
  tryCatch(
    DBI::dbExecute(con, "LOAD spatial"),
    error = function(e)
      stop(".ftwQuery(): DuckDB's spatial extension is not available locally. ",
           "ocular will not download it automatically. ", conditionMessage(e),
           call. = FALSE))
  have_icu <- !is.null(datetime_col) && tryCatch({
    DBI::dbExecute(con, "LOAD icu")
    TRUE
  }, error = function(e) FALSE)
  if (grepl("^(https?|s3)://", source, ignore.case = TRUE))
    tryCatch(
      DBI::dbExecute(con, "LOAD httpfs"),
      error = function(e)
        stop(".ftwQuery(): a remote source needs DuckDB's httpfs extension, ",
             "which is not available locally. ocular will not download it ",
             "automatically. ", conditionMessage(e), call. = FALSE))

  # read_parquet() is a table function -- its path argument cannot be a bind
  # parameter, so inject it as a quoted string literal; the numeric bbox values
  # (scalar-function args) stay parameterised.
  src_lit <- DBI::dbQuoteString(con, source)
  # Apply the bbox prefilter to the geometry extent in the dataset CRS. The
  # consumer performs the precise spatial clipping.
  geom_id <- DBI::dbQuoteIdentifier(con, geometry_col)
  dt_sel  <- if (have_icu && !is.null(datetime_col))
    paste0(", EXTRACT(epoch FROM ",
           DBI::dbQuoteIdentifier(con, datetime_col), ") AS dt_epoch") else ""
  sql <- sprintf(
    "SELECT ST_AsWKB(%1$s) AS wkb%2$s
       FROM read_parquet(%3$s)
      WHERE ST_XMax(%1$s) >= ? AND ST_XMin(%1$s) <= ?
        AND ST_YMax(%1$s) >= ? AND ST_YMin(%1$s) <= ?",
    geom_id, dt_sel, src_lit)

  res <- tryCatch(
    DBI::dbGetQuery(con, sql, params = list(
      flt[["xmin"]], flt[["xmax"]],   # ST_XMax >= xmin ; ST_XMin <= xmax
      flt[["ymin"]], flt[["ymax"]])), # ST_YMax >= ymin ; ST_YMin <= ymax
    error = function(e) {
      stop(".ftwQuery() failed for ", .safeFtwSourceLabel(source), ": ",
           .redactFtwCondition(e, source),
           "\nCheck `geometry_col`, `datetime_col`, and `crs` against the ",
           "source schema.", call. = FALSE)
    })

  if (is.null(res) || nrow(res) == 0L) return(NULL)

  # Convert the WKB BLOB column to an sf geometry in the source CRS.
  wkb <- res$wkb
  if (!is.list(wkb)) wkb <- as.list(wkb)
  geom <- sf::st_as_sfc(structure(lapply(wkb, as.raw), class = "WKB"),
                        EWKB = FALSE)
  sf::st_crs(geom) <- out_crs
  # Derive source_year when icu supplied the epoch column; otherwise use NA.
  yr <- if (have_icu && !is.null(datetime_col) && !is.null(res$dt_epoch))
    suppressWarnings(as.integer(format(
      as.POSIXct(res$dt_epoch, origin = "1970-01-01", tz = "UTC"), "%Y")))
  else rep(NA_integer_, nrow(res))
  if (!is.null(datetime_col) && !isTRUE(have_icu))
    warning("FTW: DuckDB's icu extension is unavailable locally; ",
            "source_year was set to NA.", call. = FALSE)
  sf::st_sf(source_year = yr, geometry = geom)
}


# ---- public: attach an FTW prior to an rs -----------------------------------

#' Attach an optional field boundary prior to an ocular object
#'
#' Queries a user-supplied field boundary GeoParquet, then stores the polygon
#' containing the supplied point and its representative \code{source_year} on
#' the \code{ocular} object. The default column names follow the fiboa schema.
#' Fields of The World (FTW) is one field boundary dataset that can be supplied
#' in a compatible GeoParquet; fiboa and FTW are not interchangeable names.
#'
#' The selected polygon can provide calibration support for
#' \code{segment_area()} and \code{boundary_delineation()}. When a source year
#' is available, it can also act as an age-weighted soft prior during the
#' initial area search. Use \code{diagnose_against_ftw()} to compare the
#' resulting mask with a field boundary reference.
#'
#' An attached polygon can influence the delineation. Agreement with that same
#' polygon is therefore an in-sample diagnostic, not an independent accuracy
#' assessment. Use a separate held-out reference polygon when reporting
#' accuracy.
#'
#' Whether the prior improves field boundary delineation requires independent
#' validation. This is therefore an experimental feature.
#'
#' Queries require the optional \code{duckdb} and \code{DBI} packages and a
#' locally installed DuckDB spatial extension. Remote sources also require the
#' DuckDB \code{httpfs} extension. \code{ocular} does not download or install
#' DuckDB extensions automatically.
#'
#' @param rs An \code{ocular} object returned by \code{get_rs()}.
#' @param source Path or URL to a compatible field boundary GeoParquet covering
#'   the area of interest. Local files, HTTPS objects, and S3 keys are accepted.
#' @param crs The source dataset's geometry CRS, supplied as an EPSG code or
#'   \code{sf} CRS. Set it when the source is not EPSG:4326.
#' @param refresh If \code{TRUE}, bypass the cache and repeat the query.
#' @param max_cache_entries Maximum number of entries retained in the on-disk
#'   cache. The default is 200; least-recently accessed entries are removed
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
                          datetime_col = "determination_datetime") {
  if (!is_rs(rs))
    stop("add_ftw_prior(): rs must be an ocular object.", call. = FALSE)
  if (missing(source) || !is.character(source) || length(source) != 1L ||
      is.na(source) || !nzchar(trimws(source)))
    stop("add_ftw_prior(): `source` must be a single path or URL to a field ",
         "boundary GeoParquet.", call. = FALSE)
  if (!is.logical(refresh) || length(refresh) != 1L || is.na(refresh))
    stop("add_ftw_prior(): `refresh` must be TRUE or FALSE.", call. = FALSE)
  if (!.isWholeNumber(max_cache_entries) || max_cache_entries < 1L)
    stop("add_ftw_prior(): `max_cache_entries` must be a positive integer.",
         call. = FALSE)
  max_cache_entries <- as.integer(max_cache_entries)

  source_crs <- tryCatch(sf::st_crs(crs), error = function(e) sf::st_crs(NA))
  if (is.na(source_crs))
    stop("add_ftw_prior(): `crs` is missing or invalid.", call. = FALSE)
  if (!is.character(geometry_col) || length(geometry_col) != 1L ||
      is.na(geometry_col) || !nzchar(geometry_col))
    stop("add_ftw_prior(): `geometry_col` must be one non-empty name.",
         call. = FALSE)
  if (!is.null(datetime_col) &&
      (!is.character(datetime_col) || length(datetime_col) != 1L ||
       is.na(datetime_col) || !nzchar(datetime_col)))
    stop("add_ftw_prior(): `datetime_col` must be NULL or one non-empty name.",
         call. = FALSE)

  bbox <- rs$geom$bbox
  if (is.null(bbox))
    stop("add_ftw_prior(): rs$geom$bbox is missing.", call. = FALSE)
  bbox <- c(xmin = unname(bbox[["xmin"]]), ymin = unname(bbox[["ymin"]]),
            xmax = unname(bbox[["xmax"]]), ymax = unname(bbox[["ymax"]]))

  qb <- .ftwQuantiseBbox(bbox)
  query <- list(crs_wkt = source_crs$wkt,
                geometry_col = geometry_col,
                datetime_col = datetime_col)
  key <- .ftwCacheKey(source, qb, query = query)

  polygons <- NULL
  if (!refresh) {
    hit <- .ftwCacheRead(key)
    if (!is.null(hit)) polygons <- hit$polygons
  }
  if (is.null(polygons)) {
    polygons <- .ftwQuery(bbox_wgs84 = qb, source = source, crs = source_crs,
                          geometry_col = geometry_col,
                          datetime_col = datetime_col)
    if (!is.null(polygons)) {
      .ftwCacheWrite(key, polygons,
                     meta = list(
                       source_id = .safeFtwSourceLabel(source),
                       source_hash = substr(.stableHash(source), 1L, 16L),
                       bbox = qb, query = query))
      .ftwCacheEvictLRU(max_entries = max_cache_entries)
    }
  }

  if (is.null(polygons) || nrow(polygons) == 0L) {
    message("add_ftw_prior(): no field boundaries intersect the analysis ",
            "area; ocular object unchanged.")
    return(rs)
  }

  polygons <- .selectFtwField(polygons, rs)
  if (is.null(polygons) || nrow(polygons) == 0L) {
    message("add_ftw_prior(): no field boundary contains the supplied point; ",
            "ocular object unchanged.")
    return(rs)
  }

  # The selected seed field supplies the trust-decay year.
  yr <- if (!is.null(polygons$source_year))
    suppressWarnings(stats::median(polygons$source_year, na.rm = TRUE))
  else NA_real_
  yr <- if (is.finite(yr)) as.integer(round(yr)) else NULL

  rs$geom$ftw_prior <- list(polygon = polygons, source_year = yr)
  ## get_rs() normally calibrates eagerly before the optional prior is added.
  ## Re-run now so the documented FTW support actually takes effect.
  if (!is.null(rs$internals$feature_stack)) {
    updated <- tryCatch(.calibrateScene(rs), error = identity)
    if (!inherits(updated, "error") && is.null(updated$error)) {
      rs$internals$calibration <- updated
    } else {
      msg <- if (inherits(updated, "error")) conditionMessage(updated)
             else updated$error
      warning("add_ftw_prior(): FTW calibration was unavailable; the prior ",
              "calibration was preserved. ", msg, call. = FALSE)
    }
  }
  message(sprintf("add_ftw_prior(): seed field attached%s.",
                  if (is.null(yr)) "" else sprintf(" (source_year %d)", yr)))
  rs
}
