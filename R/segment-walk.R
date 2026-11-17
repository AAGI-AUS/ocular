# =========================================================================
# Shared perimeter walker for cleanup and splitting
# =========================================================================
# The walker handles movement and iteration. Stage rules decide how each
# visited cell changes the mask.

#' Apply a rule's cell action
#'
#' Return the masks and a change flag. A matrix replacement counts as a change
#' for iteration control. A \code{noop} that only updates rule state does not.
#' @noRd
.applyAction <- function(action, alive_mat, decided_mat){

  if( is.null(action) ){
    return(list(alive_mat = alive_mat,
                decided_mat = decided_mat,
                state = NULL,
                changed = FALSE))
  }
  a <- action$action
  if( is.null(a) || identical(a, "noop") ){
    return(list(alive_mat = alive_mat, decided_mat = decided_mat,
                state = action$state, changed = FALSE))
  }
  if( identical(a, "alive_mat") ){
    return(list(alive_mat = action$new, decided_mat = decided_mat,
                state = action$state, changed = TRUE))
  }
  if( identical(a, "decided_mat") ){
    return(list(alive_mat = alive_mat, decided_mat = action$new,
                state = action$state, changed = TRUE))
  }
  if( identical(a, "both") ){
    return(list(alive_mat = action$alive_mat,
                decided_mat = action$decided_mat,
                state = action$state, changed = TRUE))
  }
  stop(".applyAction: unknown action '", a, "'.", call. = FALSE)
}

#' Apply a rule's final matrix changes
#'
#' Return the masks and a change flag. Any supplied matrix counts as a change.
#' @noRd
.applyDelta <- function(delta, alive_mat, decided_mat){

  if( is.null(delta) ){
    return(list(alive_mat = alive_mat,
                decided_mat = decided_mat,
                changed = FALSE))
  }
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

#' Run a rule around each candidate area's perimeter
#'
#' Rebuild connected areas on each iteration. Stop when no rule action replaces
#' a matrix, or when \code{max_iter(params)} is reached (one pass by default).
#'
#' Rules can set up shared \code{walker_state} for all areas and iterations.
#' Each area's \code{setup()}, \code{apply_cell()} and \code{finalize()}
#' callbacks receive it. Cleanup omits \code{apply_cell()} because it classifies
#' a complete trail at once.
#'
#' For rules with axes, finish one axis before the next. An action can set
#' \code{next_si} to resume from another step. Independent axes start from the
#' same matrices. If \code{accept_axis()} does not return scalar \code{TRUE},
#' restore that axis's matrices and rule state without recording a change.
#' @noRd
.runStagedWalk <- function(rs, rule, params, log_msg = function(...) invisible()){

  if( !any(rs$state$alive_mat) ){
    return(rs)
  }

  alive_mat <- rs$state$alive_mat
  decided_mat <- if( is.null(rs$state$decided_mat) ){
    matrix(FALSE, rs$geom$nr, rs$geom$nc)
  }else{
    rs$state$decided_mat
  }

  axes <- if( is.null(rule$axes) ) { NA_character_ } else { rule$axes }
  axis_mode <- if( is.null(rule$axis_mode) ) { "sequential" } else { rule$axis_mode }
  max_iter <- if( is.null(rule$max_iter) ){
    1L
  }else{
    max(1L, as.integer(rule$max_iter(params)))
  }

  walker_state <- if( !is.null(rule$walker_setup) ) { rule$walker_setup(rs, params) } else { NULL }

  # Both axis modes use the same pass, returning masks, rule state and change flag.
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
      alive_mat <- out$alive_mat
      decided_mat <- out$decided_mat
      if( out$changed ){
        changed <- TRUE
      }
      if( !is.null(out$state) ){
        rule_state <- out$state
      }
      if( !is.null(action$next_si) ){
        ns <- as.integer(action$next_si)
        if( !is.na(ns) && ns >= 0L && ns < n_steps ){
          si <- ns
        }
      }
    }
    return(list(alive_mat = alive_mat, decided_mat = decided_mat,
                rule_state = rule_state, changed = changed))
  }

  for( iter in seq_len(max_iter) ){
    info <- .enumerateAreas(alive_mat, rs$geom$centre_rc, strict = FALSE)
    if( info$n_areas == 0L ){
      break
    }
    iter_changed <- FALSE

    for( ai in seq_len(info$n_areas) ){
      area <- info$areas[[ai]]
      walk <- if( !is.null(rule$build_walk) ){
        rule$build_walk(area, alive_mat, decided_mat, rs, params,
                        walker_state)
      }else{
        NULL
      }
      rule_state <- rule$setup(area, walk, alive_mat, decided_mat,
                               rs, params, walker_state, log_msg)
      if( isTRUE(rule_state$skip) ){
        next
      }
      # setup() can return changes to the shared walker state.
      if( !is.null(rule_state$.walker_state) ){
        walker_state <- rule_state$.walker_state
      }

      if( !is.null(rule$apply_cell) && !is.null(walk) &&
          !is.null(walk$trail) && nrow(walk$trail) > 0L ){
        n_steps <- nrow(walk$trail)
        if( axis_mode == "independent" ){
          # Each axis starts from the same matrices. Split only removes cells and
          # marks decisions, so accepted changes cannot conflict when combined.
          alive_snap <- alive_mat
          decided_snap <- decided_mat
          per_axis_alive <- list()
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
                pass$alive_mat <- alive_snap
                pass$decided_mat <- decided_snap
                pass$rule_state <- state_snap
                pass$changed <- FALSE
              }
            }
            per_axis_alive[[axis]] <- pass$alive_mat
            per_axis_decided[[axis]] <- pass$decided_mat
            rule_state <- pass$rule_state
            if( pass$changed ){
              iter_changed <- TRUE
            }
          }
          for( axis in axes ){
            am <- per_axis_alive[[axis]]
            dm <- per_axis_decided[[axis]]
            a_chg <- am != alive_snap
            d_chg <- dm != decided_snap
            alive_mat[a_chg] <- am[a_chg]
            decided_mat[d_chg] <- dm[d_chg]
          }
        }else{
          # Each axis sees the changes made by the preceding axis.
          for( axis in axes ){
            pass <- run_axis_pass(axis, walk, n_steps,
                                  alive_mat, decided_mat, rule_state)
            alive_mat <- pass$alive_mat
            decided_mat <- pass$decided_mat
            rule_state <- pass$rule_state
            if( pass$changed ){
              iter_changed <- TRUE
            }
          }
        }
      }

      delta <- rule$finalize(rule_state, alive_mat, decided_mat,
                             area, walk, rs, params, walker_state, log_msg)
      out <- .applyDelta(delta, alive_mat, decided_mat)
      alive_mat <- out$alive_mat
      decided_mat <- out$decided_mat
      if( out$changed ){
        iter_changed <- TRUE
      }
      # finalize() can also return changes to shared state.
      if( !is.null(delta$.walker_state) ){
        walker_state <- delta$.walker_state
      }
    }

    if( !iter_changed && max_iter > 1L ){
      log_msg(sprintf("converged at iter %d", iter))
      break
    }
  }

  if( !is.null(rule$walker_teardown) ){
    td <- rule$walker_teardown(walker_state, alive_mat, decided_mat,
                               rs, params)
    if( !is.null(td$alive_mat) ){
      alive_mat <- td$alive_mat
    }
    if( !is.null(td$decided_mat) ){
      decided_mat <- td$decided_mat
    }
  }

  rs$state$alive_mat <- alive_mat
  rs$state$decided_mat <- decided_mat
  return(rs)
}
