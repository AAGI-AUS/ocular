# =========================================================================
# Core: cache ownership and deletion safety
# =========================================================================

#' Validate that a cache path is not a broad filesystem location (internal)
#' @noRd
.normaliseCachePath <- function(path, caller){
  if( !is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(trimws(path)) )
    stop(caller, ": cache directory must be one non-empty path.",
         call. = FALSE)
  out <- normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
  home <- normalizePath(path.expand("~"), winslash = "/", mustWork = FALSE)
  cwd <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  if( identical(dirname(out), out) || out %in% c(home, cwd) )
    stop(caller, ": refusing to use a filesystem root, home directory, or ",
         "working directory as a cache directory.", call. = FALSE)
  out
}

#' Expected cache ownership marker (internal)
#' @noRd
.cacheMarker <- function(cache_dir, kind, version){
  file.path(cache_dir, paste0(".ocular-", kind, "-cache-", version))
}

#' Test a cache ownership marker (internal)
#' @noRd
.cacheMarkerValid <- function(cache_dir, kind, version){
  marker <- .cacheMarker(cache_dir, kind, version)
  if( !file.exists(marker) || dir.exists(marker) ) return(FALSE)
  content <- tryCatch(readLines(marker, warn = FALSE),
                      error = function(e) character(0L))
  identical(content, c("ocular-cache", kind, version))
}

#' Create an owned cache directory and marker (internal)
#' @noRd
.initialiseCacheDir <- function(cache_dir, kind, version, caller){
  cache_dir <- .normaliseCachePath(cache_dir, caller)
  if( !dir.exists(cache_dir) ){
    ok <- dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    if( !isTRUE(ok) && !dir.exists(cache_dir) )
      stop(caller, ": could not create cache directory `", cache_dir, "`.",
           call. = FALSE)
  }
  marker <- .cacheMarker(cache_dir, kind, version)
  if( file.exists(marker) && !.cacheMarkerValid(cache_dir, kind, version) )
    stop(caller, ": cache ownership marker is invalid: `", marker, "`.",
         call. = FALSE)
  if( !file.exists(marker) ){
    ok <- tryCatch({
      writeLines(c("ocular-cache", kind, version), marker, useBytes = TRUE)
      TRUE
    }, error = function(e) FALSE)
    if( !ok || !.cacheMarkerValid(cache_dir, kind, version) )
      stop(caller, ": could not create a cache ownership marker in `",
           cache_dir, "`.", call. = FALSE)
  }
  cache_dir
}

#' Resolve a directory before a cache-clear operation (internal)
#' @noRd
.validateCacheClearDir <- function(cache_dir, kind, version, caller,
                                   default_dir){
  cache_dir <- .normaliseCachePath(cache_dir, caller)
  default_dir <- .normaliseCachePath(default_dir, caller)
  if( !dir.exists(cache_dir) ) return(cache_dir)
  marker <- .cacheMarker(cache_dir, kind, version)
  if( file.exists(marker) && !.cacheMarkerValid(cache_dir, kind, version) )
    stop(caller, ": cache ownership marker is invalid: `", marker, "`.",
         call. = FALSE)
  if( !identical(cache_dir, default_dir) &&
      !.cacheMarkerValid(cache_dir, kind, version) )
    stop(caller, ": refusing to clear an unmarked custom directory. Run the ",
         "corresponding cache-enabled workflow with this directory first.",
         call. = FALSE)
  cache_dir
}

#' List only regular cache entries with the expected key shape (internal)
#' @noRd
.cacheEntryFiles <- function(cache_dir, key_pattern){
  if( !dir.exists(cache_dir) ) return(character(0L))
  candidates <- list.files(cache_dir, pattern = key_pattern, full.names = TRUE,
                           all.files = FALSE, no.. = TRUE)
  if( length(candidates) == 0L ) return(character(0L))
  info <- file.info(candidates)
  candidates[!is.na(info$isdir) & !info$isdir]
}
