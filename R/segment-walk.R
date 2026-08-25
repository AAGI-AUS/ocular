# =========================================================================
# Shared perimeter-walking driver for segmentation stages
# Each cleanup, interior, or split rule supplies setup(), apply_cell(), and
# finalize() callbacks. .runStagedWalk applies one rule to each candidate area.
#
# Action shape returned by apply_cell:
#   list(action = "noop")
#   list(action = "alive_mat", new = updated_alive_mat)
#   list(action = "decided_mat", new = updated_decided_mat)
#   list(action = "both", alive_mat = ..., decided_mat = ...)
# A rule can return a "state" entry to update its rule_state.
#
# Delta shape returned by finalize:
#   list(alive_mat = updated_alive_mat or NULL,
#        decided_mat = updated_decided_mat or NULL)
# =========================================================================

#' Apply a per-cell rule action to alive_mat / decided_mat (internal)
#'
#' Returns updated mats plus a \code{changed} flag derived from the action
#' type: any action that mutates alive_mat or decided_mat counts as a
#' change for convergence purposes; pure-state \code{noop} actions do not.
#'
#' @noRd
.applyAction <- function(action, alive_mat, decided_mat){
  if( is.null(action) ) return(list(alive_mat = alive_mat,
                                    decided_mat = decided_mat,
                                    state   = NULL,
                                    changed = FALSE))
  a <- action$action
  if( is.null(a) || identical(a, "noop") )
    return(list(alive_mat = alive_mat, decided_mat = decided_mat,
                state = action$state, changed = FALSE))
  if( identical(a, "alive_mat") )
    return(list(alive_mat = action$new, decided_mat = decided_mat,
                state = action$state, changed = TRUE))
  if( identical(a, "decided_mat") )
    return(list(alive_mat = alive_mat, decided_mat = action$new,
                state = action$state, changed = TRUE))
  if( identical(a, "both") )
    return(list(alive_mat = action$alive_mat,
                decided_mat = action$decided_mat,
                state = action$state, changed = TRUE))
  stop(".applyAction: unknown action '", a, "'.", call. = FALSE)
}

#' Apply a finalize delta to alive_mat / decided_mat (internal)
#'
#' Returns updated mats plus a \code{changed} flag: TRUE iff the delta
#' supplied a non-NULL alive_mat or decided_mat replacement.
#'
#' @noRd
.applyDelta <- function(delta, alive_mat, decided_mat){
  if( is.null(delta) ) return(list(alive_mat = alive_mat,
                                   decided_mat = decided_mat,
                                   changed = FALSE))
  changed <- FALSE
  if( !is.null(delta$alive_mat) ){
    alive_mat <- delta$alive_mat
    changed <- TRUE
  }
  if( !is.null(delta$decided_mat) ){
    decided_mat <- delta$decided_mat
    changed <- TRUE
  }
  return(list(alive_mat = alive_mat, decided_mat = decided_mat,
              changed = changed))
}

#' Drive a single rule through its per-area perimeter walks (internal)
#'
#' Applies a segmentation rule across candidate-area perimeter walks.
#'
#' \strong{Iter loop and convergence.} Each rule declares
#' \code{max_iter(params)} (defaults to 1L). The driver re-enumerates
#' connected components at the top of each iter and breaks when an iter
#' completes with no changes (any action returning "alive_mat" /
#' "decided_mat" / "both", or any finalize delta with non-NULL mats,
#' counts as a change). For \code{max_iter = 1L} the loop runs once
#' regardless of changes -- the convergence check is benign.
#'
#' \strong{Walker state.} Rules can declare \code{walker_setup(rs, params)}
#' and \code{walker_teardown(state, alive_mat, decided_mat, rs, params)}
#' callbacks to manage state that lives across all iters and all areas
#' (for example, cleanup's shared area_mat). The driver passes the resulting
#' \code{walker_state} to setup, apply_cell, and finalize.
#'
#' \strong{Axes and rewinds.} A rule may declare \code{axes} (character
#' vector) -- the per-step apply_cell is then called once per
#' (axis, step) in axis-outer order. Per-step actions may include
#' \code{next_si} to control the step iterator (used by split_area's
#' single-area rewind). In independent mode, an optional
#' \code{accept_axis()} callback evaluates each completed provisional pass.
#' A result other than scalar \code{TRUE} restores the pre-axis matrices and
#' rule state and does not count as a convergence change.
#'
#' \strong{Finalize-only rules.} A rule may set \code{apply_cell = NULL};
#' the driver skips the per-cell loop and calls setup followed by finalize.
#' Cleanup uses this form because \code{.classifyTrail} evaluates a complete
#' perimeter at once.
#'
#' @noRd
.runStagedWalk <- function(rs, rule, params, log_msg = function(...) invisible()){

  if( !any(rs$state$alive_mat) ) return(rs)

  alive_mat   <- rs$state$alive_mat
  decided_mat <- if( is.null(rs$state$decided_mat) )
    matrix(FALSE, rs$geom$nr, rs$geom$nc) else rs$state$decided_mat

  axes       <- if( is.null(rule$axes) )      NA_character_ else rule$axes
  axis_mode  <- if( is.null(rule$axis_mode) ) "sequential"   else rule$axis_mode
  max_iter   <- if( is.null(rule$max_iter) )   1L
  else max(1L, as.integer(rule$max_iter(params)))

  walker_state <- if( !is.null(rule$walker_setup) )
    rule$walker_setup(rs, params)
  else NULL

  ## Inner loop body -- runs one axis pass to completion. Closed over so
  ## the sequential and independent branches can both invoke it without
  ## duplication. Returns the post-pass (alive, decided, state, changed).
  run_axis_pass <- function(axis, walk, n_steps, alive_mat, decided_mat,
                            rule_state){
    changed <- FALSE
    si <- 0L
    while( si < n_steps ){
      si <- si + 1L
      action <- rule$apply_cell(si, axis, walk, rule_state,
                                alive_mat, decided_mat, rs, params,
                                walker_state)
      out <- .applyAction(action, alive_mat, decided_mat)
      alive_mat   <- out$alive_mat
      decided_mat <- out$decided_mat
      if( out$changed ) changed <- TRUE
      if( !is.null(out$state) ) rule_state <- out$state
      ## Optional rule-controlled rewind (split_area single-area branch).
      if( !is.null(action$next_si) ){
        ns <- as.integer(action$next_si)
        if( !is.na(ns) && ns >= 0L && ns < n_steps )
          si <- ns
      }
    }
    list(alive_mat = alive_mat, decided_mat = decided_mat,
         rule_state = rule_state, changed = changed)
  }

  for( iter in seq_len(max_iter) ){
    info <- .enumerateAreas(alive_mat, rs$geom$centre_rc, strict = FALSE)
    if( info$n_areas == 0L ) break
    iter_changed <- FALSE

    for( ai in seq_len(info$n_areas) ){
      area <- info$areas[[ai]]
      walk <- if( !is.null(rule$build_walk) )
        rule$build_walk(area, alive_mat, decided_mat, rs, params,
                        walker_state)
      else NULL
      rule_state <- rule$setup(area, walk, alive_mat, decided_mat,
                               rs, params, walker_state, log_msg)
      if( isTRUE(rule_state$skip) ) next
      ## walker_state may be mutated by setup (e.g. cleanup's prev_cells
      ## bookkeeping). Pull updated copy if rule_state surfaces it.
      if( !is.null(rule_state$.walker_state) )
        walker_state <- rule_state$.walker_state

      ## Per-cell loop only when apply_cell is declared and there's a walk.
      if( !is.null(rule$apply_cell) && !is.null(walk) &&
          !is.null(walk$trail) && nrow(walk$trail) > 0L ){
        n_steps <- nrow(walk$trail)
        if( axis_mode == "independent" ){
          ## Snapshot at axis-loop entry; each axis evaluates against the
          ## snapshot. Per-axis deltas are merged back via cell-wise
          ## difference from the snapshot. Per-cell agreement is
          ## guaranteed for monotonic rules (split_area only sets cells
          ## FALSE / decided TRUE; segment_interior is handled outside the
          ## driver via .runDetectObjects).
          alive_snap   <- alive_mat
          decided_snap <- decided_mat
          per_axis_alive   <- list()
          per_axis_decided <- list()
          for( axis in axes ){
            state_snap <- rule_state
            pass <- run_axis_pass(axis, walk, n_steps,
                                  alive_snap, decided_snap, rule_state)
            if( pass$changed && !is.null(rule$accept_axis) ){
              accepted <- isTRUE(rule$accept_axis(
                axis = axis,
                alive_before = alive_snap,
                decided_before = decided_snap,
                alive_after = pass$alive_mat,
                decided_after = pass$decided_mat,
                state = pass$rule_state,
                area = area,
                walk = walk,
                rs = rs,
                params = params,
                walker_state = walker_state
              ))
              if( !accepted ){
                pass$alive_mat   <- alive_snap
                pass$decided_mat <- decided_snap
                pass$rule_state  <- state_snap
                pass$changed     <- FALSE
              }
            }
            per_axis_alive[[axis]]   <- pass$alive_mat
            per_axis_decided[[axis]] <- pass$decided_mat
            rule_state               <- pass$rule_state
            if( pass$changed ) iter_changed <- TRUE
          }
          for( axis in axes ){
            am <- per_axis_alive[[axis]]
            dm <- per_axis_decided[[axis]]
            a_chg <- am != alive_snap
            d_chg <- dm != decided_snap
            alive_mat[a_chg]   <- am[a_chg]
            decided_mat[d_chg] <- dm[d_chg]
          }
        }else{
          ## Sequential: each axis runs against the current alive_mat;
          ## mutations from one axis are visible to the next.
          for( axis in axes ){
            pass <- run_axis_pass(axis, walk, n_steps,
                                  alive_mat, decided_mat, rule_state)
            alive_mat   <- pass$alive_mat
            decided_mat <- pass$decided_mat
            rule_state  <- pass$rule_state
            if( pass$changed ) iter_changed <- TRUE
          }
        }
      }

      delta <- rule$finalize(rule_state, alive_mat, decided_mat,
                             area, walk, rs, params, walker_state, log_msg)
      out <- .applyDelta(delta, alive_mat, decided_mat)
      alive_mat   <- out$alive_mat
      decided_mat <- out$decided_mat
      if( out$changed ) iter_changed <- TRUE
      ## finalize may also surface walker_state updates.
      if( !is.null(delta$.walker_state) )
        walker_state <- delta$.walker_state
    }

    if( !iter_changed && max_iter > 1L ){
      log_msg(sprintf("converged at iter %d", iter))
      break
    }
  }

  if( !is.null(rule$walker_teardown) ){
    td <- rule$walker_teardown(walker_state, alive_mat, decided_mat,
                               rs, params)
    if( !is.null(td$alive_mat) )   alive_mat   <- td$alive_mat
    if( !is.null(td$decided_mat) ) decided_mat <- td$decided_mat
  }

  rs$state$alive_mat   <- alive_mat
  rs$state$decided_mat <- decided_mat
  return(rs)
}
