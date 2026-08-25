# =========================================================================
# Pipeline construction -- user-defined reusable rs pipelines
# =========================================================================

#' Build a reusable ocular pipeline
#'
#' Captures a chain of supported processing and output calls as a callable
#' \code{ocular_pipeline} object. Retrieval is not part of the captured chain:
#' first create an \code{ocular} object with \code{get_rs()}, then pass that
#' object to the resulting function. Values that can be evaluated when the
#' pipeline is built are stored and reused on each run.
#'
#' Supported stages are \code{trace_perimeter()}, \code{segment_area()},
#' \code{segment_interior()}, \code{split_area()},
#' \code{boundary_delineation()}, \code{add_modis()}, \code{as_raster()}, and
#' \code{as_time_series()}.
#'
#' @param pipeline A native-pipe expression composed of supported stage calls.
#'   Leave the initial \code{ocular} object out of the expression; it is
#'   supplied when the resulting function is called.
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
      identical(pipeline_expr, quote(expr =)) )
    stop("build_pipeline: pipeline argument is required.", call. = FALSE)

  stage_fns <- c("trace_perimeter", "segment_area", "segment_interior",
                 "split_area", "boundary_delineation", "add_modis",
                 "as_raster", "as_time_series")

  if( !is.call(pipeline_expr) ||
      !is.symbol(pipeline_expr[[1L]]) ||
      !as.character(pipeline_expr[[1L]]) %in% stage_fns )
    stop("build_pipeline: pipeline must be a chain of stage functions (",
         paste(stage_fns, collapse = ", "), "). Got: ",
         deparse1(pipeline_expr), call. = FALSE)

  ## When the user does not supply params, do not materialise a fresh
  ## rs_params() snapshot -- let each stage fall back to rs$params at call
  ## time so stored settings and their explicit-ownership metadata are used.
  captured_params <- params  # NULL stays NULL
  spec_expr <- .buildPipelineSpec(pipeline_expr, stage_fns, capture_envir)

  fn <- function(rs, params = NULL){
    if( !is_rs(rs) )
      stop("ocular pipeline input must be an ocular object.", call. = FALSE)
    eval_params <- if( is.null(params) ) captured_params else params

    result <- eval(spec_expr,
                   envir  = list(rs = rs, params = eval_params),
                   enclos = parent.frame())

    result
  }
  attr(fn, "input_expr")      <- pipeline_expr
  attr(fn, "spec_expr")       <- spec_expr
  attr(fn, "captured_params") <- captured_params
  class(fn) <- c("ocular_pipeline", "function")

  return(fn)
}

#' Walk a pipeline expression: capture-time eval of non-stage args,
#' inject rs at the deepest stage call, inject params = params at every
#' stage call that doesn't already specify one.
#'
#' Returns the rewritten expression.
#'
#' @noRd
.buildPipelineSpec <- function(expr, stage_fns, capture_envir){

  injected_rs <- FALSE

  walk <- function(e){
    if( !is.call(e) ) return(e)
    fn_name  <- if( is.symbol(e[[1L]]) ) as.character(e[[1L]]) else NULL
    is_stage <- !is.null(fn_name) && fn_name %in% stage_fns

    if( !is_stage ){
      ## Non-stage call inside a stage's named arg (e.g. mean(x)) --
      ## evaluate at capture time to freeze the value.
      return( tryCatch(eval(e, envir = capture_envir),
                       error = function(err) e) )
    }

    args <- as.list(e)[-1L]
    arg_names <- if( is.null(names(args)) )
      rep("", length(args)) else names(args)
    has_positional <- any(arg_names == "")

    args <- lapply(seq_along(args), function(i){
      a <- args[[i]]
      is_pos <- arg_names[i] == ""
      if( is.call(a) ){
        inner_fn <- if( is.symbol(a[[1L]]) )
          as.character(a[[1L]]) else NULL
        is_inner_stage <- !is.null(inner_fn) && inner_fn %in% stage_fns
        if( is_inner_stage ) return(walk(a))
        if( is_pos )
          stop("build_pipeline: non-stage call '", deparse1(a[[1L]]),
               "()' in piping position. Allowed stages: ",
               paste(stage_fns, collapse = ", "), call. = FALSE)
      }
      tryCatch(eval(a, envir = capture_envir),
               error = function(err) a)
    })
    names(args) <- arg_names

    if( !has_positional && !injected_rs ){
      args <- c(list(quote(rs)), args)
      injected_rs <<- TRUE
    }
    if( !"params" %in% names(args) )
      args[["params"]] <- quote(params)

    as.call(c(list(e[[1L]]), args))
  }
  spec <- walk(expr)
  return(spec)
}

#' Pretty-print a captured pipeline expression with pipe syntax
#' @noRd
.deparsePipeline <- function(expr){

  if( !is.call(expr) ) return(deparse1(expr))
  args <- as.list(expr)[-1L]
  arg_names <- if( is.null(names(args)) )
    rep("", length(args)) else names(args)
  pos_idx <- which(arg_names == "")

  ## Recursive case: first positional arg is itself a call -> emit pipe.
  if( length(pos_idx) >= 1L && is.call(args[[pos_idx[1L]]]) ){
    lhs <- .deparsePipeline(args[[pos_idx[1L]]])
    rhs_args <- args[-pos_idx[1L]]
    rhs <- deparse1(as.call(c(list(expr[[1L]]), rhs_args)))
    return( paste0(lhs, " |> ", rhs) )
  }
  return(deparse1(expr))
}

#' Format a params list with one entry per line, R-standard indented
#' @noRd
.formatParamsList <- function(p){

  if( !is.list(p) || length(p) == 0L )
    return(paste0("    ", deparse1(p)))
  nms <- names(p)
  entries <- vapply(seq_along(p), function(i){
    paste0("      ", nms[i], " = ", deparse1(p[[i]]))
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
