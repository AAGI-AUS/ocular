# =========================================================================
# Reusable processing pipelines
# =========================================================================

#' Build a reusable ocular pipeline
#'
#' Turns a chain of processing and output calls into a reusable function.
#' Retrieve scenes with \code{get_rs()}, then pass the resulting \code{ocular}
#' object to this function. Stage arguments are evaluated when the pipeline is
#' built where possible, then reused on each run.
#'
#' Supported stages are \code{trace_perimeter()}, \code{segment_area()},
#' \code{segment_interior()}, \code{split_area()},
#' \code{boundary_delineation()}, \code{add_modis()}, \code{as_raster()}, and
#' \code{as_time_series()}.
#'
#' @param pipeline A native-pipe expression using the supported stages.
#'   It may start with a stage, \code{rs}, or a name holding an \code{ocular}
#'   object. Each form uses the input passed when the pipeline is called.
#'   Other unbound names at the head are errors. A positional stage setting
#'   follows the supplied input. Retrieval calls and raster inputs are not
#'   supported.
#' @param params Optional default \code{rs_params()} settings to store with the
#'   pipeline. \code{NULL} uses the settings on each input object. A
#'   \code{params} argument supplied when the resulting function is called
#'   replaces these defaults for stages that do not supply their own.
#' @returns A callable function with class \code{ocular_pipeline}.
#' @seealso \code{\link{get_rs}}, \code{\link{rs_params}}
#' @export
build_pipeline <- function(pipeline, params = NULL){

  pipeline_expr <- substitute(pipeline)
  capture_envir <- parent.frame()
  if( missing(pipeline) ||
      identical(pipeline_expr, quote(expr = )) ){
    stop("build_pipeline: pipeline argument is required.", call. = FALSE)
  }

  stage_fns <- c("trace_perimeter", "segment_area", "segment_interior",
                 "split_area", "boundary_delineation", "add_modis",
                 "as_raster", "as_time_series")

  if( !is.call(pipeline_expr) ||
      !is.symbol(pipeline_expr[[1L]]) ||
      !as.character(pipeline_expr[[1L]]) %in% stage_fns ){
    stop("build_pipeline: pipeline must be a chain of stage functions (",
         paste(stage_fns, collapse = ", "), "). Got: ",
         deparse1(pipeline_expr), call. = FALSE)
  }

  # Keep NULL so each stage uses the input object's stored settings.
  captured_params <- params
  spec_expr <- .buildPipelineSpec(pipeline_expr, stage_fns, capture_envir)

  fn <- function(rs, params = NULL){

    if( !is_rs(rs) ){
      stop("ocular pipeline input must be an ocular object.", call. = FALSE)
    }
    invocation <- .beginCalibrationInvocation(rs)
    rs <- invocation$rs
    if( invocation$owner ){
      on.exit({
        invocation$context$active <- FALSE
        invocation$context$failures <- NULL
      }, add = TRUE)
    }

    eval_params <- if( is.null(params) ) { captured_params } else { params }

    # Use package functions when the caller has no binding for a stage name.
    # This also works when ocular is loaded but not attached.
    eval_list <- list(rs = rs, params = eval_params)
    for( stage_fn in stage_fns ){
      if( !exists(stage_fn, envir = parent.frame(), mode = "function") ){
        eval_list[[stage_fn]] <- get(stage_fn, envir = asNamespace("ocular"),
                                     mode = "function")
      }
    }

    result <- eval(spec_expr,
                   envir = eval_list,
                   enclos = parent.frame())

    return(.endCalibrationInvocation(result, invocation))
  }
  attr(fn, "input_expr") <- pipeline_expr
  attr(fn, "spec_expr") <- spec_expr
  attr(fn, "captured_params") <- captured_params
  class(fn) <- c("ocular_pipeline", "function")

  return(fn)
}

#' Prepare a stage expression for reuse
#'
#' Capture argument values when the pipeline is built. Insert \code{rs} at
#' the first stage and \code{params} wherever a stage does not supply its own.
#'
#' @returns The rewritten expression.
#' @noRd
.buildPipelineSpec <- function(expr, stage_fns, capture_envir){

  injected_rs <- FALSE

  walk <- function(e){

    if( !is.call(e) ){
      return(e)
    }
    fn_name <- if( is.symbol(e[[1L]]) ) { as.character(e[[1L]]) } else { NULL }
    is_stage <- !is.null(fn_name) && fn_name %in% stage_fns

    if( !is_stage ){
      # Capture expressions inside stage arguments now, such as mean(x).
      return(tryCatch(eval(e, envir = capture_envir),
                      error = function(err){
                        return(e)
                      }))
    }

    args <- as.list(e)[-1L]
    arg_names <- if( is.null(names(args)) ) { rep("", length(args)) } else { names(args) }
    positional <- which(arg_names == "")
    head_idx <- if( length(positional) ) { positional[1L] } else { 0L }
    head_expr <- if( head_idx > 0L ) { args[[head_idx]] } else { NULL }
    has_inner_stage <- any(vapply(args, function(a){
      return(is.call(a) && is.symbol(a[[1L]]) &&
               as.character(a[[1L]]) %in% stage_fns)
    }, logical(1L)))

    args <- lapply(seq_along(args), function(i){
      a <- args[[i]]
      is_pos <- arg_names[i] == ""
      if( i == head_idx && !has_inner_stage && is.symbol(a) ){
        if( identical(a, quote(rs)) ){
          return(a)
        }
        if( !exists(as.character(a), envir = capture_envir, inherits = TRUE) ){
          stop("build_pipeline: '", as.character(a),
               "' has no value when the pipeline is built. ",
               "Use rs for the input, or define the stage setting first.", call. = FALSE)
        }
      }
      if( is.call(a) ){
        inner_fn <- if( is.symbol(a[[1L]]) ) { as.character(a[[1L]]) } else { NULL }
        is_inner_stage <- !is.null(inner_fn) && inner_fn %in% stage_fns
        if( is_inner_stage ){
          return(walk(a))
        }
        if( is_pos ){
          stop("build_pipeline: non-stage call '", deparse1(a[[1L]]),
               "()' in piping position. Allowed stages: ",
               paste(stage_fns, collapse = ", "), call. = FALSE)
        }
      }
      return(tryCatch(eval(a, envir = capture_envir),
                      error = function(err){
                        return(a)
                      }))
    })
    names(args) <- arg_names

    if( !has_inner_stage && !injected_rs ){
      if( any(arg_names %in% c("longitude", "latitude", "start_date", "end_date")) ){
        stop("build_pipeline: retrieval arguments are not supported. ",
             "Retrieve an ocular object before building the pipeline.", call. = FALSE)
      }
      head_value <- if( head_idx > 0L ) { args[[head_idx]] } else { NULL }
      if( inherits(head_value, "SpatRaster") ){
        stop("build_pipeline: a raster cannot be the pipeline input. ",
             "Apply trace_perimeter() to the raster directly.", call. = FALSE)
      }
      if( is.symbol(head_expr) &&
          (identical(head_expr, quote(rs)) || is_rs(head_value)) ){
        args[head_idx] <- list(quote(rs))
      }else{
        args <- c(list(quote(rs)), args)
      }
      injected_rs <<- TRUE
    }
    if( !"params" %in% names(args) ){
      args[["params"]] <- quote(params)
    }

    return(as.call(c(list(e[[1L]]), args)))
  }
  spec <- walk(expr)
  return(spec)
}

#' Pretty-print a captured pipeline expression with pipe syntax
#' @noRd
.deparsePipeline <- function(expr){

  if( !is.call(expr) ){
    return(deparse1(expr))
  }
  args <- as.list(expr)[-1L]
  arg_names <- if( is.null(names(args)) ) { rep("", length(args)) } else { names(args) }
  pos_idx <- which(arg_names == "")

  if( length(pos_idx) >= 1L && is.call(args[[pos_idx[1L]]]) ){
    lhs <- .deparsePipeline(args[[pos_idx[1L]]])
    rhs_args <- args[-pos_idx[1L]]
    rhs <- deparse1(as.call(c(list(expr[[1L]]), rhs_args)))
    return(paste0(lhs, " |> ", rhs))
  }
  return(deparse1(expr))
}

#' Format a parameter list with one entry per line
#' @noRd
.formatParamsList <- function(p){

  if( !is.list(p) || length(p) == 0L ){
    return(paste0("    ", deparse1(p)))
  }
  nms <- names(p)
  entries <- vapply(seq_along(p), function(i){
    return(paste0("      ", nms[i], " = ", deparse1(p[[i]])))
  }, character(1L))
  return(paste(c("    list(",
                 paste(entries, collapse = ",\n"),
                 "    )"),
               collapse = "\n"))
}

#' Print an ocular pipeline
#'
#' Displays the captured stage chain, the expression used when the pipeline is
#' called, and any stored parameter settings. The pipeline is not evaluated.
#'
#' @param x An \code{ocular_pipeline} object.
#' @param ... Ignored.
#' @returns \code{x} invisibly.
#' @method print ocular_pipeline
#' @export
print.ocular_pipeline <- function(x, ...){

  cat("<ocular_pipeline>\n")
  cat("  input:\n    ", .deparsePipeline(attr(x, "input_expr")), "\n", sep = "")
  cat("  spec:\n    ",  .deparsePipeline(attr(x, "spec_expr")),  "\n", sep = "")
  cat("  params:\n",   .formatParamsList(attr(x, "captured_params")), "\n", sep = "")
  return(invisible(x))
}
