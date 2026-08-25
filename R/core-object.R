# =========================================================================
# Core: rs object
# =========================================================================

#' Test whether an object is an ocular object
#' @param x Any R object.
#' @returns Logical.
#' @examples
#' is_rs(list())
#' is_rs(structure(list(), class = "ocular"))
#' @export
is_rs <- function(x) inherits(x, "ocular")

#' Construct an ocular object
#'
#' Constructs an empty rs with the standard grouped structure. Most users
#' build one via \code{get_rs()} which fills in spec, geom, and scenes from
#' a STAC fetch. Direct \code{.newOcular()} construction is useful for tests.
#'
#' Structure:
#' \itemize{
#'   \item \code{spec}: user-supplied query (point, dates, dimensions, index, filters).
#'   \item \code{geom}: derived spatial properties (bbox, source, pixel size, dims, centre cell).
#'   \item \code{scenes}: per-date band rasters + computed index_name index + STAC asset hrefs.
#'   \item \code{internals}: internal detection data (feature_stack, feat_array, mu, detection_scenes).
#'   \item \code{state}: matrix carriers (alive_mat, decided_mat).
#'   \item \code{params}: full populated rs_params list; auto-applied by every pipeline stage.
#' }
#'
#' Pre-cleanup pipeline state is carried as \code{attr(rs, "pending")} --
#' kept off the visible structure because it is execution state rather than
#' user-facing.
#'
#' @returns An empty ocular object whose primary class is \code{"ocular"}.
#' @noRd
.newOcular <- function(){

  obj <- list(
    spec = list(
      point           = list(longitude = NA_real_, latitude = NA_real_),
      start_date      = NA_character_,
      end_date        = NA_character_,
      x_metres        = NA_real_,
      y_metres        = NA_real_,
      index_name      = NA_character_,
      max_cloud_cover = NA_real_,
      scl_classes     = NULL
    ),
    geom = list(
      bbox          = NULL,
      source        = NA_character_,
      pixel_size_m  = NA_real_,
      nr            = NA_integer_,
      nc            = NA_integer_,
      centre_rc     = NULL,
      ## Optional FTW reference geometry (polygon, confidence, source-year).
      ftw_prior     = NULL
    ),
    scenes = NULL,
    internals = list(
      feature_stack    = NULL,
      feat_array       = NULL,
      mu               = NULL,
      detection_scenes = NULL,
      ## Versioned detection-window calibration record: support source,
      ## mu_window/disp_window signature, gated baselines, support provenance,
      ## and warnings.
      calibration      = NULL
    ),
    state = list(
      alive_mat   = NULL,
      decided_mat = NULL
    ),
    params = list()
  )
  class(obj) <- "ocular"
  attr(obj, "pending") <- list(pre_cleanup = NULL)
  return(obj)
}

#' Validate the structural invariants of an ocular object
#'
#' Stops on inconsistency. Checks: presence of every top-level group;
#' shape of point sub-list; alive_mat / decided_mat dimensions match
#' geom$nr x geom$nc when populated; scenes and detection_scenes are
#' lists of \code{list(date, bands, index, item_assets)} when populated;
#' params is a named list containing recognised parameter names.
#'
#' @param rs An ocular object.
#' @returns The rs invisibly, on success.
#' @noRd
.validateOcular <- function(rs){

  required <- c("spec", "geom", "scenes", "internals", "state", "params")
  missing <- setdiff(required, names(rs))
  if( length(missing) )
    stop("rs missing top-level groups: ", paste(missing, collapse = ", "),
         call. = FALSE)

  pt <- rs$spec$point
  if( !is.list(pt) || !all(c("longitude", "latitude") %in% names(pt)) )
    stop("rs$spec$point must be list(longitude, latitude).", call. = FALSE)

  nr <- rs$geom$nr; nc <- rs$geom$nc
  has_grid_state <- !is.null(rs$state$alive_mat) ||
    !is.null(rs$state$decided_mat)
  if( has_grid_state &&
      (!.isWholeNumber(nr) || nr < 1L || !.isWholeNumber(nc) || nc < 1L) )
    stop("rs grid dimensions must be positive integers when state matrices exist.",
         call. = FALSE)
  if( !is.null(rs$state$alive_mat) ){
    if( !is.matrix(rs$state$alive_mat) || !is.logical(rs$state$alive_mat) ||
        anyNA(rs$state$alive_mat) ||
        !identical(dim(rs$state$alive_mat), c(as.integer(nr), as.integer(nc))) )
      stop("rs$state$alive_mat dim does not match geom$nr x geom$nc.", call. = FALSE)
  }
  if( !is.null(rs$state$decided_mat) ){
    if( !is.matrix(rs$state$decided_mat) ||
        !is.logical(rs$state$decided_mat) || anyNA(rs$state$decided_mat) ||
        !identical(dim(rs$state$decided_mat), c(as.integer(nr), as.integer(nc))) )
      stop("rs$state$decided_mat dim does not match geom$nr x geom$nc.", call. = FALSE)
  }

  validate_scenes <- function(x, label){
    if( is.null(x) ) return(invisible())
    if( !is.list(x) )
      stop(label, " must be NULL or a list of scenes.", call. = FALSE)
    required_scene <- c("date", "bands", "index", "item_assets")
    for( i in seq_along(x) ){
      sc <- x[[i]]
      if( !is.list(sc) || !all(required_scene %in% names(sc)) )
        stop(label, "[[", i, "]] must contain date, bands, index, and item_assets.",
             call. = FALSE)
      date_type_ok <- inherits(sc$date, "Date") || is.character(sc$date)
      date_ok <- date_type_ok && tryCatch({
        d <- suppressWarnings(as.Date(sc$date))
        length(d) == 1L && !is.na(d)
      }, error = function(e) FALSE)
      if( !date_ok ||
          !is.list(sc$bands) || !inherits(sc$index, "SpatRaster") ||
          !is.list(sc$item_assets) )
        stop(label, "[[", i, "]] has an invalid scene field.", call. = FALSE)
    }
    invisible()
  }
  validate_scenes(rs$scenes, "rs$scenes")
  validate_scenes(rs$internals$detection_scenes,
                  "rs$internals$detection_scenes")

  if( !is.list(rs$params) )
    stop("rs$params must be a list.", call. = FALSE)
  if( length(rs$params) > 0L ){
    param_names <- names(rs$params)
    if( is.null(param_names) || anyNA(param_names) || any(!nzchar(param_names)) ||
        anyDuplicated(param_names) )
      stop("rs$params must have unique, non-empty names.", call. = FALSE)
    if( is.null(attr(rs$params, "user_set", exact = TRUE)) )
      stop("rs$params is missing explicit parameter ownership metadata. ",
           "Recreate the ocular object with get_rs().", call. = FALSE)
    .paramUserSet(rs$params, ".validateOcular()")
    unknown <- setdiff(param_names, names(rs_params()))
    if( length(unknown) > 0L )
      stop("rs$params contains unknown parameter(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
  }

  if( !is.null(rs$internals$calibration) && !is.list(rs$internals$calibration) )
    stop("rs$internals$calibration must be NULL or a list.", call. = FALSE)
  if( !is.null(rs$geom$ftw_prior) && !is.list(rs$geom$ftw_prior) )
    stop("rs$geom$ftw_prior must be NULL or a list.", call. = FALSE)

  return(invisible(rs))
}

#' Derive a UTM EPSG string from an ocular object's analysis centre
#'
#' Six-degree zone derivation, including the standard Norway and Svalbard
#' exceptions: hemisphere (326x north / 327x south) plus zone number.
#' Computed on demand rather than cached on \code{rs}.
#'
#' @param rs An ocular object.
#' @returns A character EPSG string (e.g. \code{"EPSG:32750"}).
#' @export
rs_utm <- function(rs){

  if( !is_rs(rs) )
    stop("rs_utm(): rs must be an ocular object.", call. = FALSE)
  pt <- rs$spec$point
  if( !is.list(pt) || !.isScalarNumber(pt$longitude) ||
      !.isScalarNumber(pt$latitude) )
    stop("rs_utm(): rs$spec$point must contain finite longitude/latitude scalars.",
         call. = FALSE)
  if( pt$longitude < -180 || pt$longitude > 180 )
    stop("rs_utm(): longitude must be in [-180, 180].", call. = FALSE)
  if( pt$latitude < -80 || pt$latitude > 84 )
    stop("rs_utm(): latitude ", pt$latitude,
         " is outside the UTM-defined band [-80, 84]. Use a polar projection ",
         "for this location.", call. = FALSE)
  zone <- if( pt$longitude == 180 ) 60L
  else as.integer(floor((pt$longitude + 180) / 6) + 1L)
  ## EPSG/UTM special zones for southwest Norway and Svalbard.
  if( pt$latitude >= 56 && pt$latitude < 64 &&
      pt$longitude >= 3 && pt$longitude < 12 )
    zone <- 32L
  if( pt$latitude >= 72 && pt$latitude < 84 ){
    if( pt$longitude >= 0  && pt$longitude < 9  ) zone <- 31L
    if( pt$longitude >= 9  && pt$longitude < 21 ) zone <- 33L
    if( pt$longitude >= 21 && pt$longitude < 33 ) zone <- 35L
    if( pt$longitude >= 33 && pt$longitude < 42 ) zone <- 37L
  }
  return(paste0("EPSG:",
                if( pt$latitude >= 0 ) "326" else "327",
                sprintf("%02d", zone)))
}

#' Derive analysis centre xy in target UTM (internal)
#'
#' On-demand replacement for the dropped \code{rs$centre_v} cache.
#' Projects spec$point to rs's UTM and returns the xy coordinates.
#'
#' @noRd
.centreXY <- function(rs){

  pt <- rs$spec$point
  centre_v <- terra::project(
    terra::vect(cbind(pt$longitude, pt$latitude), crs = "EPSG:4326"),
    rs_utm(rs))
  return(terra::geom(centre_v)[1L, c("x", "y")])
}

#' Print method for an ocular object
#'
#' Compact summary: query line, source line, state line, params line.
#' Avoids dumping full rasters / matrices / arrays, which the default
#' list printer otherwise floods the terminal with.
#'
#' @param x An ocular object.
#' @param ... Ignored.
#' @returns \code{x} invisibly.
#' @method print ocular
#' @export
print.ocular <- function(x, ...){

  cat("<ocular \u2014 sample-free field boundary delineation>\n")

  pt <- x$spec$point
  if( .isScalarNumber(pt$longitude) && .isScalarNumber(pt$latitude) ){
    cat(sprintf("  query:    %.4f, %.4f | %sm x %sm | %s to %s\n",
                pt$longitude, pt$latitude,
                format(x$spec$x_metres), format(x$spec$y_metres),
                x$spec$start_date, x$spec$end_date))
  }

  if( is.character(x$geom$source) && length(x$geom$source) == 1L &&
      !is.na(x$geom$source) && nzchar(x$geom$source) ){
    n_scenes <- length(x$scenes %||% list())
    cat(sprintf("  source:   %s | %s | %d scenes\n",
                x$geom$source,
                if( is.na(x$spec$index_name) ) "(no index)" else x$spec$index_name,
                n_scenes))
  }

  state_parts <- character(0L)
  if( !is.null(x$state$alive_mat) ){
    state_parts <- c(state_parts, sprintf("%d alive px", sum(x$state$alive_mat)))
  }
  if( !is.null(attr(x, "pending")$pre_cleanup) ){
    state_parts <- c(state_parts, "pre-cleanup pending")
  }
  if( length(state_parts) > 0L )
    cat(sprintf("  state:    %s\n", paste(state_parts, collapse = " | ")))

  ## Identify non-default params by diffing against rs_params() defaults.
  if( length(x$params) > 0L ){
    defaults  <- rs_params()
    differs   <- vapply(names(x$params), function(nm){
      !identical(x$params[[nm]], defaults[[nm]])
    }, logical(1L))
    user_mods <- names(x$params)[differs]
    if( length(user_mods) > 0L ){
      preview <- if( length(user_mods) <= 4L ) paste(user_mods, collapse = ", ")
      else paste(c(user_mods[1:4], sprintf("...(+%d)",
                                           length(user_mods) - 4L)),
                 collapse = ", ")
      cat(sprintf("  params:   %d non-default (%s)\n",
                  length(user_mods), preview))
    }
  }

  return(invisible(x))
}

#' Detailed summary for an ocular object
#'
#' Verbose breakdown of every populated group. Use when \code{print.ocular}'s
#' compact form is not enough.
#'
#' @param object An ocular object.
#' @param ... Ignored.
#' @returns \code{object} invisibly.
#' @method summary ocular
#' @export
summary.ocular <- function(object, ...){

  print(object)
  cat("\n--- spec ---\n")
  utils::str(object$spec, max.level = 2L, give.attr = FALSE)
  cat("\n--- geom ---\n")
  utils::str(object$geom, max.level = 1L, give.attr = FALSE)
  if( length(object$scenes %||% list()) > 0L )
    cat(sprintf("\n--- scenes --- %d entries; first scene bands: %s\n",
                length(object$scenes),
                paste(names(object$scenes[[1L]]$bands), collapse = ", ")))
  cat("\n--- internals --- ")
  cat(sprintf("feature_stack: %s | feat_array: %s | mu: %s | detection_scenes: %s\n",
              if( is.null(object$internals$feature_stack) )    "NULL" else "set",
              if( is.null(object$internals$feat_array) )       "NULL" else "set",
              if( is.null(object$internals$mu) )               "NULL" else "set",
              if( is.null(object$internals$detection_scenes) ) "NULL" else
                sprintf("%d entries", length(object$internals$detection_scenes))))
  cat("\n--- state --- ")
  cat(sprintf("alive_mat: %s | decided_mat: %s\n",
              if( is.null(object$state$alive_mat) )   "NULL" else
                sprintf("%dx%d (%d alive)",
                        nrow(object$state$alive_mat), ncol(object$state$alive_mat),
                        sum(object$state$alive_mat)),
              if( is.null(object$state$decided_mat) ) "NULL" else "set"))
  if( length(object$params) > 0L ){
    defaults <- rs_params()
    differs  <- vapply(names(object$params), function(nm){
      !identical(object$params[[nm]], defaults[[nm]])
    }, logical(1L))
    if( any(differs) ){
      cat("\n--- params (non-default) ---\n")
      utils::str(object$params[differs], max.level = 1L, give.attr = FALSE)
    }else{
      cat("\n--- params --- all at rs_params() defaults\n")
    }
  }

  return(invisible(object))
}
