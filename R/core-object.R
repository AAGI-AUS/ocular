# =========================================================================
# Core: rs object
# =========================================================================

#' Test whether an object is an ocular object
#' @param x Any R object.
#' @returns \code{TRUE} for an \code{ocular} object, otherwise \code{FALSE}.
#' @examples
#' is_rs(list())
#' is_rs(structure(list(), class = "ocular"))
#' @export
is_rs <- function(x){

  return(inherits(x, "ocular"))
}

#' Create an empty ocular object
#'
#' \code{get_rs()} fills this structure from a query. Tests can fill it directly.
#' The groups hold:
#' \itemize{
#'   \item \code{spec}: the query point, dates, dimensions, index and filters.
#'   \item \code{geom}: spatial properties and an optional field prior.
#'   \item \code{scenes}: dated bands, output indices and STAC asset URLs.
#'   \item \code{internals}: features and calibration used for delineation.
#'   \item \code{state}: the current mask and cell decisions.
#'   \item \code{params}: stored settings for later stages.
#' }
#' Cleanup requested before a mask exists is stored in \code{attr(rs, "pending")}.
#'
#' @returns An empty \code{ocular} object.
#' @noRd
.newOcular <- function(){

  obj <- list(
    spec = list(
      point = list(longitude = NA_real_, latitude = NA_real_),
      start_date = NA_character_,
      end_date = NA_character_,
      x_metres = NA_real_,
      y_metres = NA_real_,
      index_name = NA_character_,
      max_cloud_cover = NA_real_,
      scl_classes = NULL
    ),
    geom = list(
      bbox = NULL,
      source = NA_character_,
      pixel_size_m = NA_real_,
      nr = NA_integer_,
      nc = NA_integer_,
      centre_rc = NULL,
      ftw_prior = NULL
    ),
    scenes = NULL,
    internals = list(
      feature_stack = NULL,
      feat_array = NULL,
      mu = NULL,
      detection_scenes = NULL,
      # Calibration includes its schema version, support, signatures and warnings.
      calibration = NULL
    ),
    state = list(
      alive_mat = NULL,
      decided_mat = NULL
    ),
    params = list()
  )
  class(obj) <- "ocular"
  attr(obj, "pending") <- list(pre_cleanup = NULL)
  return(obj)
}

#' Check an ocular object's structure
#'
#' Require its standard groups, the query-point structure, masks matching the grid,
#' valid scene lists and recognised parameter names. Stop on inconsistency.
#'
#' @param rs An ocular object.
#' @returns The object invisibly.
#' @noRd
.validateOcular <- function(rs){

  required <- c("spec", "geom", "scenes", "internals", "state", "params")
  missing <- setdiff(required, names(rs))
  if( length(missing) ){
    stop("rs missing top-level groups: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  pt <- rs$spec$point
  if( !is.list(pt) || !all(c("longitude", "latitude") %in% names(pt)) ){
    stop("rs$spec$point must be list(longitude, latitude).", call. = FALSE)
  }

  nr <- rs$geom$nr
  nc <- rs$geom$nc
  has_grid_state <- !is.null(rs$state$alive_mat) ||
    !is.null(rs$state$decided_mat)
  if( has_grid_state &&
      (!.isWholeNumber(nr) || nr < 1L || !.isWholeNumber(nc) || nc < 1L) ){
    stop("rs grid dimensions must be positive integers when state matrices exist.",
         call. = FALSE)
  }
  if( !is.null(rs$state$alive_mat) ){
    if( !is.matrix(rs$state$alive_mat) || !is.logical(rs$state$alive_mat) ||
        anyNA(rs$state$alive_mat) ||
        !identical(dim(rs$state$alive_mat), c(as.integer(nr), as.integer(nc))) ){
      stop("rs$state$alive_mat dim does not match geom$nr x geom$nc.", call. = FALSE)
    }
  }
  if( !is.null(rs$state$decided_mat) ){
    if( !is.matrix(rs$state$decided_mat) ||
        !is.logical(rs$state$decided_mat) || anyNA(rs$state$decided_mat) ||
        !identical(dim(rs$state$decided_mat), c(as.integer(nr), as.integer(nc))) ){
      stop("rs$state$decided_mat dim does not match geom$nr x geom$nc.", call. = FALSE)
    }
  }

  validate_scenes <- function(x, label){

    if( is.null(x) ){
      return(invisible())
    }
    if( !is.list(x) ){
      stop(label, " must be NULL or a list of scenes.", call. = FALSE)
    }
    required_scene <- c("date", "bands", "index", "item_assets")
    for( i in seq_along(x) ){
      sc <- x[[i]]
      if( !is.list(sc) || !all(required_scene %in% names(sc)) ){
        stop(label, "[[", i, "]] must contain date, bands, index, and item_assets.",
             call. = FALSE)
      }
      date_type_ok <- inherits(sc$date, "Date") || is.character(sc$date)
      date_ok <- date_type_ok && tryCatch({
        d <- suppressWarnings(as.Date(sc$date))
        length(d) == 1L && !is.na(d)
      }, error = function(e){
        return(FALSE)
      })
      if( !date_ok ||
          !is.list(sc$bands) || !inherits(sc$index, "SpatRaster") ||
          !is.list(sc$item_assets) ){
        stop(label, "[[", i, "]] has an invalid scene field.", call. = FALSE)
      }
    }
    return(invisible())
  }
  validate_scenes(rs$scenes, "rs$scenes")
  validate_scenes(rs$internals$detection_scenes,
                  "rs$internals$detection_scenes")

  if( !is.list(rs$params) ){
    stop("rs$params must be a list.", call. = FALSE)
  }
  if( length(rs$params) > 0L ){
    param_names <- names(rs$params)
    if( is.null(param_names) || anyNA(param_names) || any(!nzchar(param_names)) ||
        anyDuplicated(param_names) ){
      stop("rs$params must have unique, non-empty names.", call. = FALSE)
    }
    if( is.null(attr(rs$params, "user_set", exact = TRUE)) ){
      stop("rs$params is missing explicit parameter ownership metadata. ",
           "Recreate the ocular object with get_rs().", call. = FALSE)
    }
    .paramUserSet(rs$params, ".validateOcular()")
    unknown <- setdiff(param_names, names(rs_params()))
    if( length(unknown) > 0L ){
      stop("rs$params contains unknown parameter(s): ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
  }

  if( !is.null(rs$internals$calibration) && !is.list(rs$internals$calibration) ){
    stop("rs$internals$calibration must be NULL or a list.", call. = FALSE)
  }
  if( !is.null(rs$geom$ftw_prior) && !is.list(rs$geom$ftw_prior) ){
    stop("rs$geom$ftw_prior must be NULL or a list.", call. = FALSE)
  }

  return(invisible(rs))
}

#' Derive a UTM EPSG string from an ocular object's analysis centre
#'
#' Finds the UTM zone for the analysis centre, including the Norway and
#' Svalbard exceptions, and returns its EPSG identifier.
#'
#' @param rs An ocular object.
#' @returns A character EPSG string (e.g. \code{"EPSG:32750"}).
#' @export
rs_utm <- function(rs){

  if( !is_rs(rs) ){
    stop("rs_utm(): rs must be an ocular object.", call. = FALSE)
  }
  pt <- rs$spec$point
  if( !is.list(pt) || !.isScalarNumber(pt$longitude) ||
      !.isScalarNumber(pt$latitude) ){
    stop("rs_utm(): rs$spec$point must contain finite longitude/latitude scalars.",
         call. = FALSE)
  }
  if( pt$longitude < -180 || pt$longitude > 180 ){
    stop("rs_utm(): longitude must be in [-180, 180].", call. = FALSE)
  }
  if( pt$latitude < -80 || pt$latitude > 84 ){
    stop("rs_utm(): latitude ", pt$latitude,
         " is outside the UTM-defined band [-80, 84]. Use a polar projection ",
         "for this location.", call. = FALSE)
  }
  zone <- if( pt$longitude == 180 ){
    60L
  }else{
    as.integer(floor((pt$longitude + 180) / 6) + 1L)
  }
  # EPSG/UTM special zones for southwest Norway and Svalbard.
  if( pt$latitude >= 56 && pt$latitude < 64 &&
      pt$longitude >= 3 && pt$longitude < 12 ){
    zone <- 32L
  }
  if( pt$latitude >= 72 && pt$latitude < 84 ){
    if( pt$longitude >= 0 && pt$longitude < 9 ){
      zone <- 31L
    }
    if( pt$longitude >= 9 && pt$longitude < 21 ){
      zone <- 33L
    }
    if( pt$longitude >= 21 && pt$longitude < 33 ){
      zone <- 35L
    }
    if( pt$longitude >= 33 && pt$longitude < 42 ){
      zone <- 37L
    }
  }
  return(paste0("EPSG:",
                if( pt$latitude >= 0 ) { "326" } else { "327" },
                sprintf("%02d", zone)))
}

#' Project the query point into the analysis UTM zone
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
#' Shows the query, imagery source and delineation mask state without printing
#' the underlying raster data. Stored parameters supplied by the user or changed
#' from defaults are shown. User-supplied settings are marked as user input.
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
                if( is.na(x$spec$index_name) ) { "(no index)" } else { x$spec$index_name },
                n_scenes))
  }

  state_parts <- if( is.null(x$state$alive_mat) ){
    "no delineation mask"
  }else if( !any(x$state$alive_mat) ){
    "empty delineation mask"
  }else{
    sprintf("%d cells in delineation mask", sum(x$state$alive_mat))
  }
  if( isTRUE(attr(x, "pending")$pre_cleanup$active) ){
    state_parts <- c(state_parts, "pre-cleanup pending")
  }
  cat(sprintf("  state:    %s\n", paste(state_parts, collapse = " | ")))

  if( length(x$params) > 0L ){
    defaults <- rs_params()
    differs <- vapply(names(x$params), function(nm){
      return(!identical(x$params[[nm]], defaults[[nm]]))
    }, logical(1L))
    user_set <- attr(x$params, "user_set", exact = TRUE) %||% character(0L)
    user_mods <- names(x$params)[differs | names(x$params) %in% user_set]
    if( length(user_mods) > 0L ){
      user_mods <- ifelse(user_mods %in% user_set,
                          paste0(user_mods, " (user input)"), user_mods)
      preview <- if( length(user_mods) <= 4L ){
        paste(user_mods, collapse = ", ")
      }else{
        paste(c(user_mods[1:4], sprintf("...(+%d)",
                                        length(user_mods) - 4L)),
              collapse = ", ")
      }
      cat(sprintf("  stored parameters: %s\n", preview))
    }
  }

  return(invisible(x))
}

#' Detailed summary for an ocular object
#'
#' Adds query details, geometry, scene counts, processing state and stored
#' parameter values to the compact \code{print()} summary. Stored parameters
#' do not record every override supplied to later calls.
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
  if( length(object$scenes %||% list()) > 0L ){
    cat(sprintf("\n--- scenes --- %d entries. First scene bands: %s\n",
                length(object$scenes),
                paste(names(object$scenes[[1L]]$bands), collapse = ", ")))
  }
  cat("\n--- internals --- ")
  cat(sprintf("feature_stack: %s | feat_array: %s | mu: %s | detection_scenes: %s\n",
              if( is.null(object$internals$feature_stack) ) { "NULL" } else { "set" },
              if( is.null(object$internals$feat_array) ) { "NULL" } else { "set" },
              if( is.null(object$internals$mu) ) { "NULL" } else { "set" },
              if( is.null(object$internals$detection_scenes) ){
                "NULL"
              }else{
                sprintf("%d entries", length(object$internals$detection_scenes))
              }))
  cat("\n--- state --- ")
  cat(sprintf("alive_mat: %s | decided_mat: %s\n",
              if( is.null(object$state$alive_mat) ){
                "NULL"
              }else{
                sprintf("%dx%d (%d alive)",
                        nrow(object$state$alive_mat), ncol(object$state$alive_mat),
                        sum(object$state$alive_mat))
              },
              if( is.null(object$state$decided_mat) ) { "NULL" } else { "set" }))
  if( length(object$params) > 0L ){
    defaults <- rs_params()
    differs <- vapply(names(object$params), function(nm){
      return(!identical(object$params[[nm]], defaults[[nm]]))
    }, logical(1L))
    user_set <- attr(object$params, "user_set", exact = TRUE)
    shown <- names(object$params)[differs | names(object$params) %in% user_set]
    cat("\n--- stored parameters ---\n")
    if( length(shown) > 0L ){
      for( nm in shown ){
        value <- paste(deparse(object$params[[nm]], width.cutoff = 100L), collapse = " ")
        flag <- if( nm %in% user_set ) { " (user input)" } else { "" }
        cat(sprintf("%s: %s%s\n", nm, value, flag))
      }
    }else if( is.null(user_set) ){
      cat("All stored values match defaults. User input was not recorded.\n")
    }else{
      cat("All defaults, no user input.\n")
    }
    cat("Stored parameters do not record every override supplied to later calls.\n")
  }

  return(invisible(object))
}
