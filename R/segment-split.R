# =========================================================================
# Segmentation: area splitting
# =========================================================================

#' Linear (LOS) cut along envelope axis (internal)
#'
#' Returns the LOS row through commit_pos, terminated at first-shot
#' envelope on each side.
#'
#' @noRd
.linearCut <- function(commit_pos, envelope, nr, nc){

  empty_rc <- matrix(integer(0), ncol = 2L,
                     dimnames = list(NULL, c("row", "col")))
  near <- envelope$near_rc; far <- envelope$far_rc
  return(if( envelope$axis == "y" ){
    idx <- which(near[, "row"] == commit_pos[1L])
    if( length(idx) == 0L ) return(empty_rc)
    idx <- idx[1L]
    lo <- min(near[idx, "col"], far[idx, "col"])
    hi <- max(near[idx, "col"], far[idx, "col"])
    cbind(row = rep(commit_pos[1L], hi - lo + 1L), col = lo:hi)
  }else{
    idx <- which(near[, "col"] == commit_pos[2L])
    if( length(idx) == 0L ) return(empty_rc)
    idx <- idx[1L]
    lo <- min(near[idx, "row"], far[idx, "row"])
    hi <- max(near[idx, "row"], far[idx, "row"])
    cbind(row = lo:hi, col = rep(commit_pos[2L], hi - lo + 1L))
  })
}

#' Nonlinear cut envelope (internal)
#'
#' Spine perpendicular to commit_dir; per spine step, scan
#' perpendicular-left for first alive cell within \code{sens} of mu_seed;
#' trail cell = one-removed back. Terminates at first-shot envelope.
#'
#' @noRd
.nonlinearCut <- function(commit_pos, commit_dir, envelope,
                          alive_mat, feat_array, mu_seed, sens, nr, nc){

  empty_rc <- matrix(integer(0), ncol = 2L,
                     dimnames = list(NULL, c("row", "col")))
  dir_delta <- function(d) switch(EXPR = d,
                                  N = c(-1L, 0L), S = c(1L, 0L), E = c(0L, 1L), W = c(0L, -1L))
  left_of <- function(d) switch(EXPR = d, N = "W", S = "E", E = "N", W = "S")

  ## fs_mat: TRUE on first-shot envelope cells (near or far, both axes).
  fs_mat <- matrix(FALSE, nr, nc)
  fs_mat[envelope$near_rc] <- TRUE
  fs_mat[envelope$far_rc]  <- TRUE

  vi_at <- function(r, c){
    if( length(dim(feat_array)) == 3L ){
      strip <- feat_array[r, c, ]
      if( length(strip) == 1L ) strip else median(strip, na.rm = TRUE)
    }else{
      feat_array[r, c]
    }
  }

  sd_d <- dir_delta(left_of(commit_dir))
  sc_d <- dir_delta(left_of(left_of(commit_dir)))

  cut_rc <- empty_rc
  spine_pos <- as.integer(commit_pos)
  max_steps <- max(nr, nc) + 4L
  for( s in seq_len(max_steps) ){
    spine_pos <- spine_pos + sd_d
    if( spine_pos[1L] < 1L || spine_pos[1L] > nr ||
        spine_pos[2L] < 1L || spine_pos[2L] > nc ) break
    if( fs_mat[spine_pos[1L], spine_pos[2L]] ){
      cut_rc <- rbind(cut_rc, matrix(spine_pos, nrow = 1L,
                                     dimnames = list(NULL, c("row", "col"))))
      break
    }
    scan_pos <- spine_pos
    terminated <- FALSE
    for( k in seq_len(max_steps) ){
      scan_pos <- scan_pos + sc_d
      if( scan_pos[1L] < 1L || scan_pos[1L] > nr ||
          scan_pos[2L] < 1L || scan_pos[2L] > nc ) break
      if( fs_mat[scan_pos[1L], scan_pos[2L]] ){
        cut_rc <- rbind(cut_rc, matrix(scan_pos, nrow = 1L,
                                       dimnames = list(NULL, c("row", "col"))))
        terminated <- TRUE; break
      }
      if( alive_mat[scan_pos[1L], scan_pos[2L]] ){
        v <- vi_at(scan_pos[1L], scan_pos[2L])
        if( is.finite(v) && abs(v - mu_seed) <= sens ){
          tcell <- scan_pos - sc_d
          cut_rc <- rbind(cut_rc, matrix(tcell, nrow = 1L,
                                         dimnames = list(NULL, c("row", "col"))))
          terminated <- TRUE; break
        }
      }
    }
    if( !terminated ) break
    if( fs_mat[cut_rc[nrow(cut_rc), "row"], cut_rc[nrow(cut_rc), "col"]] ) break
  }
  return(unique(cut_rc))
}

#' Reduce a split seed signature to one finite scalar (internal)
#'
#' @noRd
.splitSeedSignature <- function(mu){
  mu_values <- as.numeric(mu)
  if( length(mu_values) == 0L || !any(is.finite(mu_values)) )
    return(NULL)
  mu_seed <- if( length(mu_values) > 1L )
    median(mu_values, na.rm = TRUE) else mu_values
  if( length(mu_seed) != 1L || !is.finite(mu_seed) ) return(NULL)
  mu_seed
}

#' Protected point-connected split component (internal)
#'
#' Marks finite alive cells whose scalar VI value is within \code{sens} of the
#' seed signature, then returns the complete 8-connected marked component
#' containing \code{centre_rc}. Matching is inclusive. Returns \code{NULL} when
#' the supplied-point cell cannot anchor that component without fabricating
#' spectral evidence.
#'
#' @noRd
.splitProtectedComponent <- function(alive_mat, feat_array, mu_seed, sens,
                                     centre_rc){

  if( !is.matrix(alive_mat) || !is.logical(alive_mat) || anyNA(alive_mat) )
    return(NULL)
  feat_dim <- dim(feat_array)
  if( !(length(feat_dim) %in% c(2L, 3L)) ||
      !identical(as.integer(feat_dim[1:2]), as.integer(dim(alive_mat))) )
    return(NULL)
  if( length(mu_seed) != 1L || !is.finite(mu_seed) ||
      length(sens) != 1L || !is.finite(sens) || sens < 0 )
    return(NULL)
  if( length(centre_rc) != 2L || anyNA(centre_rc) ||
      !all(vapply(centre_rc, .isWholeNumber, logical(1L))) )
    return(NULL)

  centre_rc <- as.integer(centre_rc)
  nr <- nrow(alive_mat); nc <- ncol(alive_mat)
  if( centre_rc[1L] < 1L || centre_rc[1L] > nr ||
      centre_rc[2L] < 1L || centre_rc[2L] > nc ||
      !isTRUE(alive_mat[centre_rc[1L], centre_rc[2L]]) )
    return(NULL)

  cells <- which(alive_mat, arr.ind = TRUE)
  values <- vapply(seq_len(nrow(cells)), function(i){
    .viAt(feat_array, cells[i, "row"], cells[i, "col"])
  }, numeric(1L))
  finite <- is.finite(values)
  matching <- matrix(FALSE, nr, nc)
  if( any(finite) ){
    finite_cells <- cells[finite, , drop = FALSE]
    matching[finite_cells] <-
      abs(values[finite] - mu_seed) <= sens
  }
  if( !isTRUE(matching[centre_rc[1L], centre_rc[2L]]) ) return(NULL)

  info <- .enumerateAreas(matching, centre_rc, strict = FALSE)
  centre_label <- info$label_mat[centre_rc[1L], centre_rc[2L]]
  if( centre_label < 1L ) return(NULL)
  return(info$label_mat == centre_label)
}

#' Split-area rule module (internal)
#'
#' Rule-module form of split_area's perimeter-walking logic.
#' Declares \code{axes = c("y", "x")} with \code{axis_mode = "independent"}:
#' the driver snapshots alive_mat/decided_mat at axis-loop entry and
#' runs each axis against the snapshot, merging per-axis cuts back at
#' the end. y-axis cuts do not influence the x-axis pass's input.
#' Per-axis state (consecutive divergence counter, flag list) lives in
#' \code{rule_state$by_axis[[axis]]}.
#'
#' Single-area commit semantics: returns \code{action$next_si =
#' first$step_idx} so the driver rewinds the step iterator; the next
#' iteration resumes one past the first flag. Multi-area commits
#' fall through without rewind.
#'
#' Cuts mutate both alive_mat (FALSE) and decided_mat (TRUE) along
#' the cut row, threaded back through \code{action = "both"}. Both
#' operations are monotonic so the per-axis cell-wise merge is
#' conflict-free (a cell cut by either accepted axis is dead and decided).
#' Each axis is provisional: \code{accept_axis()} rejects an axis that removes
#' any member of the protected supplied-point component.
#'
#' @noRd
.splitAreaRule <- function(){
  list(
    name      = "split_area",
    axes      = c("y", "x"),
    axis_mode = "independent",

    build_walk = function(area, alive_mat, decided_mat, rs, params, walker_state){
      if( area$size < 3L ) return(NULL)
      .areaWalk(area, rs$geom$nr, rs$geom$nc)
    },

    setup = function(area, walk, alive_mat, decided_mat, rs, params, walker_state, log_msg){
      if( is.null(walk) || nrow(walk$trail) < 2L )
        return(list(skip = TRUE))
      mu_seed <- .splitSeedSignature(rs$internals$mu)
      if( is.null(mu_seed) )
        return(list(skip = TRUE))
      protected <- .splitProtectedComponent(
        alive_mat, rs$internals$feat_array, mu_seed,
        params$split_sensitivity, rs$geom$centre_rc
      )
      if( is.null(protected) )
        return(list(skip = TRUE))
      list(
        skip        = FALSE,
        envelopes   = .firstShotEnvelope(area$cells, rs$geom$nr, rs$geom$nc),
        mu_seed     = mu_seed,
        sens        = params$split_sensitivity,
        thresh      = params$split_gate,
        is_multi    = isTRUE(params$multiple_areas),
        linear_mode = isTRUE(params$split_linear),
        nr          = rs$geom$nr,
        nc          = rs$geom$nc,
        feat_array  = rs$internals$feat_array,
        protected_component = protected,
        ## Per-axis state slots. Each axis traverses independently;
        ## the y-axis loop only writes to by_axis$y, etc. Empty at
        ## start of each axis loop because nothing has populated them.
        by_axis = list(y = list(consecutive = 0L, flags = list()),
                       x = list(consecutive = 0L, flags = list())),
        log_msg = log_msg
      )
    },

    apply_cell = function(si, axis, walk, state, alive_mat, decided_mat, rs, params, walker_state){
      ## Direction requires a preceding position, so traversal begins at the
      ## second perimeter step.
      if( si < 2L ) return(list(action = "noop"))
      env <- state$envelopes[[axis]]
      if( is.null(env) ) return(list(action = "noop"))

      clear_streak <- function(){
        ax <- state$by_axis[[axis]]
        if( ax$consecutive == 0L && length(ax$flags) == 0L )
          return(list(action = "noop"))
        state$by_axis[[axis]]$consecutive <- 0L
        state$by_axis[[axis]]$flags       <- list()
        list(action = "noop", state = state)
      }

      trail <- walk$trail
      pos_now  <- trail[si, ]
      pos_prev <- trail[si - 1L, ]
      dir_now  <- .deltaToDir(pos_now[1L] - pos_prev[1L],
                              pos_now[2L] - pos_prev[2L])
      if( is.na(dir_now) ) return(clear_streak())

      los <- .scanLOS(env, pos_now, dir_now, alive_mat,
                      state$feat_array, state$nr, state$nc)
      if( los$skip ) return(clear_streak())

      diverged <- abs(los$median_alive - state$mu_seed) >= state$sens
      ax <- state$by_axis[[axis]]

      if( !diverged ){
        ## A measurable non-divergent profile breaks the qualifying streak.
        return(clear_streak())
      }

      ## Accumulate flag.
      ax$consecutive <- ax$consecutive + 1L
      ax$flags <- c(ax$flags, list(list(rc = as.integer(pos_now),
                                        dir = dir_now,
                                        step_idx = si)))
      if( ax$consecutive < state$thresh ){
        state$by_axis[[axis]] <- ax
        return(list(action = "noop", state = state))
      }

      ## Threshold reached -- commit.
      n_flags <- length(ax$flags)
      first <- ax$flags[[1L]]
      last  <- ax$flags[[n_flags]]

      commit_one <- function(pos, dir){
        if( state$linear_mode ) .linearCut(pos, env, state$nr, state$nc)
        else                    .nonlinearCut(pos, dir, env,
                                              alive_mat, state$feat_array,
                                              state$mu_seed, state$sens,
                                              state$nr, state$nc)
      }

      changed <- FALSE

      if( !state$is_multi ){
        ## Single-area: one cut at the first flag, then rewind.
        c1 <- commit_one(first$rc, first$dir)
        if( nrow(c1) > 0L ){
          alive_mat[c1]   <- FALSE
          decided_mat[c1] <- TRUE
          changed <- TRUE
        }
        state$by_axis[[axis]]$consecutive <- 0L
        state$by_axis[[axis]]$flags       <- list()
        if( changed )
          return(list(action  = "both",
                      alive_mat = alive_mat,
                      decided_mat = decided_mat,
                      state    = state,
                      next_si  = first$step_idx))
        return(list(action  = "noop",
                    state   = state,
                    next_si = first$step_idx))
      }

      ## Multi-area: 1st + nth + intermediates (always linear). No rewind.
      c1 <- commit_one(first$rc, first$dir)
      if( nrow(c1) > 0L ){
        alive_mat[c1]   <- FALSE
        decided_mat[c1] <- TRUE
        changed <- TRUE
      }
      if( n_flags > 1L ){
        cN <- commit_one(last$rc, last$dir)
        if( nrow(cN) > 0L ){
          alive_mat[cN]   <- FALSE
          decided_mat[cN] <- TRUE
          changed <- TRUE
        }
        if( n_flags > 2L ){
          for( i in 2:(n_flags - 1L) ){
            mid <- ax$flags[[i]]
            cI  <- .linearCut(mid$rc, env, state$nr, state$nc)
            if( nrow(cI) > 0L ){
              alive_mat[cI]   <- FALSE
              decided_mat[cI] <- TRUE
              changed <- TRUE
            }
          }
        }
      }
      state$log_msg(sprintf("axis=%s commit @ step %d (%d flags)",
                            axis, si, n_flags))
      state$by_axis[[axis]]$consecutive <- 0L
      state$by_axis[[axis]]$flags       <- list()
      if( changed )
        return(list(action      = "both",
                    alive_mat   = alive_mat,
                    decided_mat = decided_mat,
                    state       = state))
      return(list(action = "noop", state = state))
    },

    accept_axis = function(axis, alive_before, decided_before,
                           alive_after, decided_after, state, area, walk,
                           rs, params, walker_state){
      protected <- state$protected_component
      if( is.null(protected) ||
          !identical(dim(protected), dim(alive_before)) )
        return(FALSE)
      !any(protected & alive_before & !alive_after)
    },

    finalize = function(state, alive_mat, decided_mat, area, walk, rs, params, walker_state, log_msg){
      return(NULL)
    },

    postprocess = function(rs, params, log_msg){
      ## Rank-1 trim: cuts can split the area into multiple components.
      ## Honours area_separation_strict.
      if( isFALSE(params$multiple_areas) && any(rs$state$alive_mat) ){
        info <- .applyStrictSeparation(rs$state$alive_mat, rs$geom$centre_rc,
                                       params$area_separation_strict)
        if( info$n_areas > 1L ){
          keep_label <- info$area_order[1L]
          rs$state$alive_mat <- info$label_mat == keep_label
          log_msg(sprintf("Single-area: keeping rank-1 area (%d px)",
                          sum(rs$state$alive_mat)))
        }
      }
      log_msg(sprintf("split_area: %d alive", sum(rs$state$alive_mat)))
      return(rs)
    }
  )
}

#' Separate adjoining candidate fields
#'
#' Looks for consecutive measurable line-of-sight profiles whose alive-cell
#' medians differ from the seed signature. A profile qualifies when its
#' difference meets or exceeds \code{split_sensitivity}; \code{split_gate}
#' controls how many qualifying profiles are required before a cut. Each axis
#' is evaluated provisionally. An axis is rejected if it removes any cell in
#' the finite, 8-connected seed-signature-matching component containing the
#' supplied point. If a present detection state cannot establish that
#' component, the candidate mask is returned unchanged. Use
#' \code{multiple_areas = TRUE} to retain disconnected components in the
#' output.
#'
#' Parameters stored in \code{rs$params} provide the base settings. Named
#' arguments supplied here override the corresponding setting for this call
#' only. The protected-component rule applies to provisional axis cuts. An
#' explicit stricter connectivity setting can repartition diagonally connected
#' cells during later single-area postprocessing, while the component
#' containing the supplied-point cell remains rank 1.
#'
#' @param x An \code{ocular} object in piped form, or longitude in standalone
#'   form.
#' @param params Optional, possibly partial \code{rs_params()} list. Overlays
#'   \code{rs$params} for this call. \code{NULL} keeps the stored settings.
#' @param split_sensitivity Minimum difference between a line-of-sight
#'   alive-pixel median and the seed signature that counts as divergence.
#'   Finite alive cells within this distance, using an inclusive comparison,
#'   can form the protected point-connected component. Exact equality therefore
#'   qualifies for both comparisons. \code{NULL} disables splitting.
#' @param split_gate Number of consecutive measurable divergent perimeter
#'   steps required before a cut is made. A measurable non-divergent or
#'   unmeasurable profile clears the streak.
#' @param split_linear Whether to use straight line-of-sight cuts rather than
#'   the default curved cuts.
#' @param strictness Calibration control. It affects initial area growth when
#'   \code{segment_area()} must run and updates the calibration record for
#'   subsequent stages; it does not directly change the split test. See
#'   \code{?rs_params}.
#' @param ... Standalone mode forwards remaining arguments to \code{get_rs()}.
#'   Piped mode accepts further named \code{rs_params()} overrides for this
#'   call, including \code{multiple_areas}; unknown names error.
#' @returns An \code{ocular} object with its revised candidate field mask stored
#'   in \code{state$alive_mat}.
#' @export
split_area <- function(x,
                       split_sensitivity = NULL,
                       split_gate        = NULL,
                       strictness        = NULL,
                       split_linear      = NULL,
                       params            = NULL,
                       ...){
  disable_split <- !missing(split_sensitivity) &&
    is.null(split_sensitivity)
  formal_overrides <- list(
    split_sensitivity = split_sensitivity,
    split_gate        = split_gate,
    split_linear      = split_linear,
    strictness        = strictness
  )
  if( is_rs(x) ){
    rs <- x
    params <- .resolveParams(rs, params, c(formal_overrides, list(...)))
  }else{
    rs <- get_rs(longitude = x, params = params, ...)
    params <- .resolveParams(rs, NULL, formal_overrides)
  }
  if( disable_split ){
    params["split_sensitivity"] <- list(NULL)
    attr(params, "user_set") <- union(
      attr(params, "user_set", exact = TRUE), "split_sensitivity")
  }
  rs <- .applyScenePriors(rs, params)
  if( is.null(rs$state$alive_mat) ){
    rs <- segment_area(rs, params = params)
  }
  if( is.null(params$split_sensitivity) ) return(rs)
  if( !any(rs$state$alive_mat) ) return(rs)
  if( is.null(rs$internals$feat_array) || is.null(rs$internals$mu) )
    stop("split_area(): detection state is unavailable; run segment_area() ",
         "first.", call. = FALSE)

  mu_seed <- .splitSeedSignature(rs$internals$mu)
  if( is.null(mu_seed) ) return(rs)
  protected <- .splitProtectedComponent(
    rs$state$alive_mat, rs$internals$feat_array, mu_seed,
    params$split_sensitivity, rs$geom$centre_rc
  )
  if( is.null(protected) ) return(rs)

  log_msg <- function(...) invisible()

  rule <- .splitAreaRule()
  rs   <- .runStagedWalk(rs, rule, params, log_msg)
  rs   <- rule$postprocess(rs, params, log_msg)
  return(rs)
}
