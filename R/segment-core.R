# =========================================================================
# Segmentation: signature learning and flood fill (matrix-based)
# =========================================================================

#' Learn per-window reference signatures from the centre neighbourhood (internal)
#'
#' Centre-pixel value as fallback. Refined to median of self-consistent
#' neighbours when at least two pass both the distance test and the
#' area_threshold range \code{[area_threshold[1L], area_threshold[2L]]}.
#'
#' @noRd
.learnSignatures <- function(feat_array, centre_rc, area_sensitivity, area_threshold){

  nr <- dim(feat_array)[1L]
  nc <- dim(feat_array)[2L]
  nw <- dim(feat_array)[3L]
  r0 <- centre_rc[1L]
  c0 <- centre_rc[2L]
  thr_lo <- area_threshold[1L]; thr_hi <- area_threshold[2L]
  mu <- numeric(nw)
  for( w in seq_len(nw) ){
    sigma0 <- feat_array[r0, c0, w]
    if( !is.finite(sigma0) ){ mu[w] <- NA_real_; next }
    passing <- numeric(0)
    for( dr in -1:1 ) for( dc in -1:1 ){
      rn <- r0 + dr; cn <- c0 + dc
      if( rn < 1L || rn > nr || cn < 1L || cn > nc ) next
      val <- feat_array[rn, cn, w]
      if( !is.finite(val) ) next
      if( abs(val - sigma0) < area_sensitivity &&
          val >= thr_lo && val <= thr_hi )
        passing <- c(passing, val)
    }
    mu[w] <- if( length(passing) >= 2L ) median(passing) else sigma0
  }
  return(mu)
}

#' Constrained BFS flood fill (internal, matrix-based)
#'
#' Returns the alive_mat. Predicate precomputed once into mask_alive.
#' Reserve geometry: outer 1/8 ring is a reserve; in-reserve admissions
#' stash until the first non-reserve admission borders the reserve.
#'
#' @noRd
.floodFill <- function(feat_array, centre_rc, mu,
                       area_sensitivity, area_threshold,
                       min_windows_alive,
                       reserve_rows, reserve_cols,
                       nr, nc,
                       local_alive_density = 0.50,
                       prior_field = NULL, prior_strength = 0,
                       prior_cross = 0.75){

  nw <- dim(feat_array)[3L]
  r0 <- centre_rc[1L]; c0 <- centre_rc[2L]

  ## A scalar area_sensitivity applies to every window; a vector supplies one
  ## matching tolerance per window.
  tol_w <- if( length(area_sensitivity) == 1L ) rep(area_sensitivity, nw)
  else area_sensitivity
  if( length(tol_w) != nw )
    stop(".floodFill(): area_sensitivity must be length 1 or dim(feat_array)[3].",
         call. = FALSE)

  if( !is.null(prior_field) &&
      ( !is.matrix(prior_field) ||
        nrow(prior_field) != nr || ncol(prior_field) != nc ) )
    stop(".floodFill(): prior_field must be an (nr x nc) matrix.",
         call. = FALSE)

  ## Predicate precompute: mask_alive[r,c] TRUE where pixel passes alive test.
  ## Vectorised over (nr, nc, nw).
  thr_lo <- area_threshold[1L]; thr_hi <- area_threshold[2L]
  any_veg <- matrix(FALSE, nr, nc)
  n_match <- matrix(0L,    nr, nc)
  for( w in seq_len(nw) ){
    plane <- feat_array[, , w]
    if( is.finite(mu[w]) ){
      n_match <- n_match + (is.finite(plane) & abs(plane - mu[w]) < tol_w[w])
    }
    any_veg <- any_veg | (is.finite(plane) & plane >= thr_lo & plane <= thr_hi)
  }
  mask_alive <- any_veg & (n_match >= min_windows_alive)

  alive_mat           <- matrix(FALSE, nr, nc)
  visited_mat         <- matrix(FALSE, nr, nc)
  predicate_alive_mat <- matrix(FALSE, nr, nc)

  if( !mask_alive[r0, c0] )
    return(list(alive_mat = alive_mat, n_flagged = 0L))

  ## Reserve helpers
  in_reserve <- function(r, c){
    r < reserve_rows[1L] || r > reserve_rows[2L] ||
      c < reserve_cols[1L] || c > reserve_cols[2L]
  }
  borders_reserve <- function(r, c){
    (r - 1L) < reserve_rows[1L] || (r + 1L) > reserve_rows[2L] ||
      (c - 1L) < reserve_cols[1L] || (c + 1L) > reserve_cols[2L]
  }
  unlocked <- in_reserve(r0, c0)  # auto-unlock if seed in reserve

  ## Pre-allocated queue + stash
  max_px     <- nr * nc
  q_r        <- integer(max_px); q_c <- integer(max_px)
  q_head     <- 1L; q_tail <- 1L
  q_r[1L]    <- r0; q_c[1L] <- c0
  visited_mat[r0, c0]         <- TRUE
  alive_mat[r0, c0]           <- TRUE
  predicate_alive_mat[r0, c0] <- TRUE

  s_r <- integer(max_px); s_c <- integer(max_px); n_stash <- 0L

  hit_n <- FALSE; hit_s <- FALSE; hit_e <- FALSE; hit_w <- FALSE

  while( q_head <= q_tail ){
    r <- q_r[q_head]; c <- q_c[q_head]; q_head <- q_head + 1L

    if( !unlocked && !in_reserve(r, c) && borders_reserve(r, c) ){
      unlocked <- TRUE
      if( n_stash > 0L ){
        for( i in seq_len(n_stash) ){
          q_tail <- q_tail + 1L
          q_r[q_tail] <- s_r[i]; q_c[q_tail] <- s_c[i]
          alive_mat[s_r[i], s_c[i]] <- TRUE
        }
        n_stash <- 0L
      }
    }

    for( dr in -1:1 ) for( dc in -1:1 ){
      if( dr == 0L && dc == 0L ) next
      rn <- r + dr; cn <- c + dc
      if( rn < 1L || rn > nr || cn < 1L || cn > nc ){
        if( rn < 1L ) hit_n <- TRUE
        if( rn > nr ) hit_s <- TRUE
        if( cn < 1L ) hit_w <- TRUE
        if( cn > nc ) hit_e <- TRUE
        next
      }
      if( visited_mat[rn, cn] ) next
      visited_mat[rn, cn] <- TRUE
      if( !mask_alive[rn, cn] ) next
      predicate_alive_mat[rn, cn] <- TRUE  # phantom-alive: counts toward density

      ## Local-density admission. Denominator: visited 8-neighbours.
      ## Numerator: predicate-passing 8-neighbours (incl. density-rejected).
      r_lo <- max(1L, rn - 1L); r_hi <- min(nr, rn + 1L)
      c_lo <- max(1L, cn - 1L); c_hi <- min(nc, cn + 1L)
      vis_block   <- visited_mat        [r_lo:r_hi, c_lo:c_hi]
      alive_block <- predicate_alive_mat[r_lo:r_hi, c_lo:c_hi]
      vis_block[rn - r_lo + 1L, cn - c_lo + 1L]   <- FALSE  # exclude self
      alive_block[rn - r_lo + 1L, cn - c_lo + 1L] <- FALSE
      vis_total <- sum(vis_block)
      if( vis_total > 0L ){
        nd_ratio <- sum(alive_block) / vis_total
        admit_bar <- local_alive_density
        if( prior_strength > 0 && !is.null(prior_field) ){
          pf <- prior_field[rn, cn]
          if( !is.na(pf) && !pf && ( n_match[rn, cn] / nw ) < prior_cross )
            admit_bar <- local_alive_density + prior_strength
        }
        if( nd_ratio < admit_bar ) next
      }

      if( in_reserve(rn, cn) && !unlocked ){
        n_stash <- n_stash + 1L
        s_r[n_stash] <- rn; s_c[n_stash] <- cn
        next
      }
      q_tail <- q_tail + 1L
      q_r[q_tail] <- rn; q_c[q_tail] <- cn
      alive_mat[rn, cn] <- TRUE
    }
  }

  return(list(alive_mat = alive_mat,
              n_flagged = sum(hit_n, hit_s, hit_e, hit_w)))
}

#' Search the alive set across a sqrt(2) sub-window ladder (internal)
#'
#' Runs .floodFill on progressively larger centred sub-windows of the
#' feature stack, ascending in *sqrt(2)* steps from a 10 px floor up to
#' the supplied stack dims. Each step ~doubles the search area, matching
#' the search area by approximately a factor of two at each step. Terminates
#' early as soon as the alive set no longer touches the sub-window
#' boundary; otherwise returns the largest rung's result.
#'
#' Reserve geometry is recomputed per rung from the actual (clamped)
#' sub-window dims. Sub-window slicing is a pure-R array slice on
#' \code{feat_array}; no terra calls.
#'
#' @noRd
.searchAreaWithWidens <- function(feat_array, centre_rc, mu,
                                  area_sensitivity, area_threshold,
                                  min_windows_alive, local_alive_density,
                                  nr_full, nc_full,
                                  log_msg = function(...) invisible(),
                                  prior_field = NULL, prior_strength = 0,
                                  prior_cross = 0.75){

  sqrt2    <- sqrt(2)
  floor_px <- 10L

  ## Build ascending ladder by paired descent from supplied dims.
  ## Stop when either axis would fall below floor_px.
  rungs_nr <- as.integer(nr_full)
  rungs_nc <- as.integer(nc_full)
  nr_curr  <- rungs_nr[1L]
  nc_curr  <- rungs_nc[1L]
  while( TRUE ){
    nxt_nr <- as.integer(trunc(nr_curr / sqrt2))
    nxt_nc <- as.integer(trunc(nc_curr / sqrt2))
    if( nxt_nr < floor_px || nxt_nc < floor_px ) break
    rungs_nr <- c(nxt_nr, rungs_nr)
    rungs_nc <- c(nxt_nc, rungs_nc)
    nr_curr  <- nxt_nr
    nc_curr  <- nxt_nc
  }
  n_rungs <- length(rungs_nr)

  alive_full <- matrix(FALSE, nr_full, nc_full)
  r0 <- centre_rc[1L]; c0 <- centre_rc[2L]

  for( ri in seq_len(n_rungs) ){
    target_nr <- rungs_nr[ri]
    target_nc <- rungs_nc[ri]

    ## Centred sub-window, clamped to [1, nr_full] x [1, nc_full].
    half_r_lo <- as.integer(trunc((target_nr - 1L) / 2L))
    half_r_hi <- target_nr - 1L - half_r_lo
    half_c_lo <- as.integer(trunc((target_nc - 1L) / 2L))
    half_c_hi <- target_nc - 1L - half_c_lo
    r_lo <- max(1L, r0 - half_r_lo); r_hi <- min(nr_full, r0 + half_r_hi)
    c_lo <- max(1L, c0 - half_c_lo); c_hi <- min(nc_full, c0 + half_c_hi)
    nr_sub <- r_hi - r_lo + 1L
    nc_sub <- c_hi - c_lo + 1L

    sub_feat   <- feat_array[r_lo:r_hi, c_lo:c_hi, , drop = FALSE]
    sub_prior <- NULL
    if( !is.null(prior_field) )
      sub_prior <- prior_field[r_lo:r_hi, c_lo:c_hi, drop = FALSE]
    sub_centre <- c(r0 - r_lo + 1L, c0 - c_lo + 1L)

    reserve_rows <- c(trunc(nr_sub * 0.125) + 1L, nr_sub - trunc(nr_sub * 0.125))
    reserve_cols <- c(trunc(nc_sub * 0.125) + 1L, nc_sub - trunc(nc_sub * 0.125))

    fill <- .floodFill(feat_array          = sub_feat,
                       centre_rc           = sub_centre,
                       mu                  = mu,
                       area_sensitivity      = area_sensitivity,
                       area_threshold        = area_threshold,
                       min_windows_alive   = min_windows_alive,
                       reserve_rows        = reserve_rows,
                       reserve_cols        = reserve_cols,
                       nr                  = nr_sub,
                       nc                  = nc_sub,
                       local_alive_density = local_alive_density,
                       prior_field        = sub_prior,
                       prior_strength      = prior_strength,
                       prior_cross         = prior_cross)

    alive_full[] <- FALSE
    alive_full[r_lo:r_hi, c_lo:c_hi] <- fill$alive_mat
    log_msg(sprintf("Rung %d/%d (%dx%d): %d alive, %d edge-hits",
                    ri, n_rungs, nr_sub, nc_sub,
                    sum(fill$alive_mat), fill$n_flagged))

    if( fill$n_flagged == 0L ) return(alive_full)
  }
  return(alive_full)
}

# =========================================================================
# Connected components (matrix-based, 8 or 4 connected)
# =========================================================================

#' Enumerate connected components (internal)
#'
#' BFS over alive_mat. strict = FALSE -> 8-connected; TRUE -> 4-connected.
#'
#' @returns list(n_areas, label_mat, areas, area_order)
#'   \code{areas[[i]]}: list(rank, cells (2-col rc), size, min_dist)
#' @noRd
.enumerateAreas <- function(alive_mat, centre_rc, strict = FALSE){

  nr <- nrow(alive_mat); nc <- ncol(alive_mat)
  label_mat <- matrix(0L, nr, nc)
  if( !any(alive_mat) ){
    return(list(n_areas = 0L, label_mat = label_mat,
                areas = list(), area_order = integer(0)))
  }

  deltas <- if( isTRUE(strict) ){
    list(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L))
  }else{
    out <- list()
    for( dr in -1:1 ) for( dc in -1:1 ){
      if( dr == 0L && dc == 0L ) next
      out[[length(out) + 1L]] <- c(dr, dc)
    }
    out
  }

  alive_idx <- which(alive_mat, arr.ind = TRUE)
  n <- nrow(alive_idx)
  q_r <- integer(n); q_c <- integer(n)
  next_label <- 0L

  for( i in seq_len(n) ){
    r0 <- alive_idx[i, 1L]; c0 <- alive_idx[i, 2L]
    if( label_mat[r0, c0] != 0L ) next
    next_label <- next_label + 1L
    label_mat[r0, c0] <- next_label
    q_r[1L] <- r0; q_c[1L] <- c0
    head <- 1L; tail <- 1L
    while( head <= tail ){
      r <- q_r[head]; c <- q_c[head]; head <- head + 1L
      for( d in deltas ){
        rn <- r + d[1L]; cn <- c + d[2L]
        if( rn < 1L || rn > nr || cn < 1L || cn > nc ) next
        if( !alive_mat[rn, cn] ) next
        if( label_mat[rn, cn] != 0L ) next
        label_mat[rn, cn] <- next_label
        tail <- tail + 1L
        q_r[tail] <- rn; q_c[tail] <- cn
      }
    }
  }

  n_areas <- next_label
  areas <- vector("list", n_areas)
  for( lab in seq_len(n_areas) ){
    cells <- which(label_mat == lab, arr.ind = TRUE)
    colnames(cells) <- c("row", "col")
    dists <- pmax(abs(cells[, "row"] - centre_rc[1L]),
                  abs(cells[, "col"] - centre_rc[2L]))
    areas[[lab]] <- list(rank     = NA_integer_,
                         cells    = cells,
                         size     = nrow(cells),
                         min_dist = min(dists))
  }
  dist_vec <- vapply(areas, `[[`, numeric(1), "min_dist")
  size_vec <- vapply(areas, `[[`, integer(1), "size")
  area_order <- order(dist_vec, -size_vec)
  for( pos in seq_along(area_order) ) areas[[area_order[pos]]]$rank <- pos

  return(list(n_areas = n_areas, label_mat = label_mat,
              areas = areas, area_order = area_order))
}

#' Apply area_separation_strict (internal)
#'
#' FALSE -> identity 8-conn enumeration.
#' TRUE  -> 4-conn enumeration over all alive cells.
#' c(0, n) -> 8-conn baseline; top-n centre-most areas re-enumerated under
#'            4-conn, sub-components each become separate ranked areas.
#'
#' @noRd
.applyStrictSeparation <- function(alive_mat, centre_rc, strict_param){

  if( is.null(alive_mat) || !any(alive_mat) )
    return(.enumerateAreas(alive_mat, centre_rc, strict = FALSE))

  if( isFALSE(strict_param) )
    return(.enumerateAreas(alive_mat, centre_rc, strict = FALSE))
  if( isTRUE(strict_param) )
    return(.enumerateAreas(alive_mat, centre_rc, strict = TRUE))

  ## c(0, n)
  base8 <- .enumerateAreas(alive_mat, centre_rc, strict = FALSE)
  if( base8$n_areas <= 1L )
    return(.enumerateAreas(alive_mat, centre_rc, strict = TRUE))

  n <- as.integer(strict_param[2L])
  target_n <- min(n, base8$n_areas)
  target_labels <- base8$area_order[1:target_n]

  ## Build per-target sub-area set; non-targets pass through unchanged.
  new_areas <- list()
  for( pos in seq_along(base8$area_order) ){
    lab <- base8$area_order[pos]
    if( lab %in% target_labels ){
      sub_mat <- matrix(FALSE, nrow(alive_mat), ncol(alive_mat))
      cells <- base8$areas[[lab]]$cells
      sub_mat[cells] <- TRUE
      sub_info <- .enumerateAreas(sub_mat, centre_rc, strict = TRUE)
      for( sai in seq_len(sub_info$n_areas) ){
        new_areas[[length(new_areas) + 1L]] <- sub_info$areas[[sai]]
      }
    }else{
      new_areas[[length(new_areas) + 1L]] <- base8$areas[[lab]]
    }
  }
  dist_vec <- vapply(new_areas, `[[`, numeric(1), "min_dist")
  size_vec <- vapply(new_areas, `[[`, integer(1), "size")
  ord <- order(dist_vec, -size_vec)
  reordered <- new_areas[ord]
  for( pos in seq_along(reordered) ) reordered[[pos]]$rank <- pos

  ## Rebuild label_mat from the new partition.
  label_mat <- matrix(0L, nrow(alive_mat), ncol(alive_mat))
  for( pos in seq_along(reordered) ){
    cells <- reordered[[pos]]$cells
    label_mat[cells] <- pos
  }
  return(list(n_areas    = length(reordered),
              label_mat  = label_mat,
              areas      = reordered,
              area_order = seq_along(reordered)))
}

# =========================================================================
# Perimeter walk (matrix-based)
# =========================================================================

#' Walk the hypothetical perimeter (internal)
#'
#' Traces a closed perimeter around the alive set using a 1x1 walker
#' positioned in dead-space adjacent to the guide perimeter edge.
#'
#' Direction codes: 1 = N, 2 = E, 3 = S, 4 = W. Guide is on walker's left.
#'
#' @returns list(trail, pivots, edges, pre_scouted_empty)
#' @noRd
.walkPerimeter <- function(alive_mat, log_msg = function(...) invisible()){

  empty_trail <- matrix(integer(0), ncol = 2L,
                        dimnames = list(NULL, c("row", "col")))
  empty_out <- list(trail = empty_trail, pivots = integer(0),
                    edges = list(), pre_scouted_empty = FALSE)
  if( !any(alive_mat) ) return(empty_out)
  cells <- which(alive_mat, arr.ind = TRUE)
  if( nrow(cells) < 3L ) return(empty_out)

  ## Pad one dead ring so the NW-most alive cell always has dead space at
  ## (-1, -1): removes the row-1/col-1 start clamp that could place the
  ## walker ON the alive set for edge-touching areas. All emitted trail
  ## coordinates are de-offset by 1 on return, so callers see original-grid
  ## indices; ring positions appear as 0 / n+1 and denote dead space
  ## outside the original grid.
  alive_mat <- rbind(FALSE, cbind(FALSE, alive_mat, FALSE), FALSE)
  cells <- cells + 1L  ## shift alive coordinates into the padded grid

  nr <- nrow(alive_mat); nc <- ncol(alive_mat)
  dr_vec <- c(-1L, 0L, 1L,  0L)
  dc_vec <- c( 0L, 1L, 0L, -1L)
  right_of <- function(d) ((d %% 4L) + 1L)
  left_of  <- function(d) (((d - 2L) %% 4L) + 1L)

  is_alive <- function(r, c){
    if( r < 1L || r > nr || c < 1L || c > nc ) return(FALSE)
    alive_mat[r, c]
  }

  ## Initialise at NW-most alive + (-1, -1); face S; guide on left (E).
  ord <- order(cells[, "row"], cells[, "col"])
  p0  <- c(cells[ord[1L], "row"], cells[ord[1L], "col"])
  W0  <- c(max(1L, p0[1L] - 1L), max(1L, p0[2L] - 1L))
  init_dir <- 3L

  global_los <- function(pos, d){
    gd <- left_of(d)
    out_alive <- logical(0)
    out_r <- integer(0); out_c <- integer(0)
    limit <- max(nr, nc)
    for( s in seq_len(limit) ){
      rc <- pos[1L] + s * dr_vec[gd]
      cc <- pos[2L] + s * dc_vec[gd]
      if( rc < 1L || rc > nr || cc < 1L || cc > nc ) break
      out_alive <- c(out_alive, alive_mat[rc, cc])
      out_r <- c(out_r, rc); out_c <- c(out_c, cc)
    }
    list(alive = out_alive, r = out_r, c = out_c)
  }
  los_furthest_alive <- function(los){
    if( length(los$alive) == 0L ) return(0L)
    idx <- which(los$alive)
    if( length(idx) == 0L ) return(0L)
    max(idx)
  }
  los_any_alive_upto <- function(los, k){
    if( length(los$alive) == 0L || k < 1L ) return(FALSE)
    any(los$alive[1:min(k, length(los$alive))])
  }
  corner_dead <- function(pos, d){
    gd <- left_of(d)
    cr <- pos[1L] + dr_vec[d] + dr_vec[gd]
    cc <- pos[2L] + dc_vec[d] + dc_vec[gd]
    !is_alive(cr, cc)
  }

  trail <- matrix(W0, nrow = 1L, byrow = TRUE,
                  dimnames = list(NULL, c("row", "col")))
  pivots <- 1L
  pos <- W0
  d   <- init_dir
  G_curr  <- global_los(pos, d)
  G_prev  <- NULL
  G_pprev <- NULL
  pre_scouted_empty <- (los_furthest_alive(G_curr) == 0L)

  state <- "walking"
  scout_pos_before  <- NULL
  scout_flagged_dir <- NA_integer_

  edges <- list()
  curr_edge_has_alive <- FALSE

  pivot_seen <- matrix(FALSE, nr, nc)
  pivot_seen[W0[1L], W0[2L]] <- TRUE
  pivot_count <- 1L

  register_pivot <- function(new_pos){
    new_pivot_idx <- nrow(trail)
    edges[[length(edges) + 1L]] <<- list(from_idx  = pivots[length(pivots)],
                                         to_idx    = new_pivot_idx,
                                         has_alive = curr_edge_has_alive)
    pivots <<- c(pivots, new_pivot_idx)
    if( pivot_seen[new_pos[1L], new_pos[2L]] && pivot_count >= 3L ) return(TRUE)
    pivot_seen[new_pos[1L], new_pos[2L]] <<- TRUE
    pivot_count <<- pivot_count + 1L
    curr_edge_has_alive <<- FALSE
    FALSE
  }

  max_steps <- nrow(cells) * 12L + 200L
  for( step in seq_len(max_steps) ){
    if( state == "walking" ){
      if( corner_dead(pos, d) ){
        ar <- pos[1L] + dr_vec[d]; ac <- pos[2L] + dc_vec[d]
        if( ar < 1L || ar > nr || ac < 1L || ac > nc ){
          if( register_pivot(pos) ) break
          d <- left_of(d)
          G_pprev <- G_prev; G_prev <- NULL
          G_curr  <- global_los(pos, d)
          next
        }
        scout_pos_before  <- pos
        scout_flagged_dir <- left_of(d)
        pos <- c(ar, ac)
        trail <- rbind(trail, pos)
        G_pprev <- G_prev; G_prev <- G_curr
        G_curr  <- global_los(pos, d)
        state <- "scouting"
      }else{
        ar <- pos[1L] + dr_vec[d]; ac <- pos[2L] + dc_vec[d]
        if( ar < 1L || ar > nr || ac < 1L || ac > nc ) break
        pos <- c(ar, ac)
        trail <- rbind(trail, pos)
        if( is_alive(pos[1L], pos[2L]) ) curr_edge_has_alive <- TRUE
        G_pprev <- G_prev; G_prev <- G_curr
        G_curr  <- global_los(pos, d)
      }
    }else if( state == "scouting" ){
      cr <- pos[1L] + dr_vec[d]; cc <- pos[2L] + dc_vec[d]
      pd1 <- left_of(d); pd2 <- right_of(d)
      a_alive <- is_alive(cr, cc) ||
        is_alive(cr + dr_vec[pd1], cc + dc_vec[pd1]) ||
        is_alive(cr + dr_vec[pd2], cc + dc_vec[pd2])
      k_out  <- if( !is.null(G_pprev) ) los_furthest_alive(G_pprev) else 0L
      b_alive <- los_any_alive_upto(G_curr, k_out)
      if( !a_alive && !b_alive ){
        trail <- trail[-nrow(trail), , drop = FALSE]
        pos <- scout_pos_before
        if( register_pivot(pos) ) break
        d <- scout_flagged_dir
        G_pprev <- NULL; G_prev <- NULL
        G_curr  <- global_los(pos, d)
        state <- "walking"
      }else{
        state <- "pending"
      }
    }else if( state == "pending" ){
      any12 <- los_any_alive_upto(G_curr, 2L)
      if( any12 ){
        if( is_alive(pos[1L], pos[2L]) ) curr_edge_has_alive <- TRUE
        state <- "walking"
      }else{
        trail <- trail[-nrow(trail), , drop = FALSE]
        pos <- scout_pos_before
        if( register_pivot(pos) ) break
        d <- scout_flagged_dir
        G_pprev <- NULL; G_prev <- NULL
        G_curr  <- global_los(pos, d)
        state <- "walking"
      }
    }
  }
  if( length(pivots) > 0L && pivots[length(pivots)] < nrow(trail) ){
    edges[[length(edges) + 1L]] <- list(from_idx  = pivots[length(pivots)],
                                        to_idx    = nrow(trail),
                                        has_alive = curr_edge_has_alive)
    pivots <- c(pivots, nrow(trail))
  }

  ## De-offset trail back to original-grid coordinates (the pad ring maps
  ## to 0 / n+1: genuine dead-space walker positions outside the grid).
  ## pivots and edges index INTO the trail and are offset-free.
  trail_out <- trail
  if( nrow(trail_out) > 0L ) trail_out <- trail_out - 1L
  return(list(trail             = trail_out,
              pivots            = pivots,
              edges             = edges,
              pre_scouted_empty = pre_scouted_empty))
}

#' Perimeter cells of an area (internal)
#'
#' Per-row leftmost/rightmost + per-column topmost/bottommost. Returns
#' a 2-col rc matrix.
#'
#' @noRd
.areaPerimeter <- function(area_cells){

  if( nrow(area_cells) == 0L ) return(area_cells)
  rows <- area_cells[, "row"]; cols <- area_cells[, "col"]
  perim_idx <- logical(nrow(area_cells))
  ## per row: leftmost / rightmost col
  for( rv in unique(rows) ){
    in_row <- which(rows == rv)
    cols_r <- cols[in_row]
    perim_idx[in_row[which.min(cols_r)]] <- TRUE
    perim_idx[in_row[which.max(cols_r)]] <- TRUE
  }
  for( cv in unique(cols) ){
    in_col <- which(cols == cv)
    rows_c <- rows[in_col]
    perim_idx[in_col[which.min(rows_c)]] <- TRUE
    perim_idx[in_col[which.max(rows_c)]] <- TRUE
  }
  return(area_cells[perim_idx, , drop = FALSE])
}

#' First-shot envelope per axis (internal)
#'
#' For one area: y-axis envelope (per-row near/far cols) and x-axis
#' envelope (per-col near/far rows). Either may be NULL.
#'
#' @noRd
.firstShotEnvelope <- function(area_cells, nr, nc){

  if( nrow(area_cells) == 0L ) return(list(y = NULL, x = NULL))

  rows <- area_cells[, "row"]; cols <- area_cells[, "col"]
  unique_rows <- sort(unique(rows))
  unique_cols <- sort(unique(cols))

  y_env <- NULL
  if( length(unique_rows) > 0L ){
    y_near <- matrix(integer(0), ncol = 2L, dimnames = list(NULL, c("row", "col")))
    y_far  <- y_near
    for( rv in unique_rows ){
      cv <- cols[rows == rv]
      y_near <- rbind(y_near, c(rv, max(1L,  min(cv) - 1L)))
      y_far  <- rbind(y_far,  c(rv, min(nc, max(cv) + 1L)))
    }
    y_env <- list(axis = "y", near_rc = y_near, far_rc = y_far)
  }

  x_env <- NULL
  if( length(unique_cols) > 0L ){
    x_near <- matrix(integer(0), ncol = 2L, dimnames = list(NULL, c("row", "col")))
    x_far  <- x_near
    for( cv in unique_cols ){
      rv <- rows[cols == cv]
      x_near <- rbind(x_near, c(max(1L,  min(rv) - 1L), cv))
      x_far  <- rbind(x_far,  c(min(nr, max(rv) + 1L), cv))
    }
    x_env <- list(axis = "x", near_rc = x_near, far_rc = x_far)
  }
  return(list(y = y_env, x = x_env))
}

# =========================================================================
# Boundary trail classification (matrix-based)
# =========================================================================

#' Threshold predicate constructor (internal)
#' @noRd
.thresholdPredicates <- function(thresholds){

  low  <- thresholds[1L]; high <- thresholds[2L]
  bias <- if( length(thresholds) >= 3L ) thresholds[3L] else NA_real_
  is_below <- if( !is.na(bias) && bias == 0 ) function(r) r <= low
  else                            function(r) r < low
  is_above <- if( !is.na(bias) && bias == 1 ) function(r) r >= high
  else                            function(r) r > high
  return(list(is_below = is_below, is_above = is_above))
}

#' Classify the perimeter trail into revive/remove decisions (internal)
#'
#' Two-pass: edge classification + spur pass on alive cells not on the
#' trail. Returns 2-col matrices to_revive, to_remove. Cells with
#' decided_mat == TRUE are skipped.
#'
#' @noRd
.classifyTrail <- function(walk, alive_mat, decided_mat, perim_cells, area_cells,
                           is_below, is_above){

  empty_rc <- matrix(integer(0), ncol = 2L,
                     dimnames = list(NULL, c("row", "col")))
  trail <- walk$trail; edges <- walk$edges
  if( length(edges) == 0L || nrow(trail) == 0L )
    return(list(to_revive = empty_rc, to_remove = empty_rc))

  nr <- nrow(alive_mat); nc <- ncol(alive_mat)

  ## Trail coordinates may include padded-ring positions (row/col 0 or
  ## n+1) from .walkPerimeter for edge-touching areas. Those are dead
  ## space by construction: treat them as alive = FALSE and never as
  ## revive/remove candidates.
  in_grid <- function(r, c) r >= 1L & r <= nr & c >= 1L & c <= nc
  alive_at <- function(r, c) in_grid(r, c) && alive_mat[r, c]

  ## trail_mat: TRUE where trail visited (in-grid positions only).
  trail_mat <- matrix(FALSE, nr, nc)
  tr_keep <- in_grid(trail[, "row"], trail[, "col"])
  if( any(tr_keep) ) trail_mat[trail[tr_keep, , drop = FALSE]] <- TRUE

  ## Pre-allocated accumulators sized to a loose upper bound: every trail
  ## cell could become a revive or remove candidate in the edge pass, plus
  ## every perimeter cell could become a remove candidate in the spur pass.
  cap <- nrow(trail) + nrow(perim_cells)
  to_revive_r <- integer(cap); to_revive_c <- integer(cap); n_rev <- 0L
  to_remove_r <- integer(cap); to_remove_c <- integer(cap); n_rem <- 0L

  pending_alive_mat <- matrix(FALSE, nr, nc)

  for( ei in seq_along(edges) ){
    e <- edges[[ei]]
    fi <- e$from_idx; ti <- e$to_idx
    if( ti < fi || ti > nrow(trail) || fi < 1L ) next

    interior <- if( ti > fi + 1L ) (fi + 1L):(ti - 1L) else integer(0)
    interior_alive <- 0L
    if( length(interior) > 0L ){
      tr_int <- trail[interior, , drop = FALSE]
      int_keep <- in_grid(tr_int[, "row"], tr_int[, "col"])
      ## ring cells are dead: filtering them out leaves the count exact
      interior_alive <- sum(alive_mat[tr_int[int_keep, , drop = FALSE]])
    }

    preceding_has_alive <- if( ei > 1L ) edges[[ei - 1L]]$has_alive
    else (!isTRUE(walk$pre_scouted_empty))
    if( is.na(preceding_has_alive) ) preceding_has_alive <- FALSE
    current_has_alive <- e$has_alive
    if( is.na(current_has_alive) ) current_has_alive <- FALSE

    total_len   <- length(interior)
    total_alive <- interior_alive
    if( preceding_has_alive ){
      total_len <- total_len + 1L
      if( alive_at(trail[fi, "row"], trail[fi, "col"]) )
        total_alive <- total_alive + 1L
    }
    if( current_has_alive ){
      total_len <- total_len + 1L
      if( alive_at(trail[ti, "row"], trail[ti, "col"]) )
        total_alive <- total_alive + 1L
    }
    if( total_len == 0L ) next
    ratio <- total_alive / total_len

    cand_idx <- interior
    if( preceding_has_alive ) cand_idx <- c(fi, cand_idx)
    if( current_has_alive )   cand_idx <- c(cand_idx, ti)
    if( length(cand_idx) == 0L ) next
    cand_rc <- trail[cand_idx, , drop = FALSE]

    if( is_below(ratio) ){
      for( j in seq_len(nrow(cand_rc)) ){
        r <- cand_rc[j, "row"]; c <- cand_rc[j, "col"]
        if( !in_grid(r, c) ) next  ## padded-ring trail cell: dead space
        if( decided_mat[r, c] ) next
        if( alive_mat[r, c] ){
          n_rem <- n_rem + 1L
          to_remove_r[n_rem] <- r; to_remove_c[n_rem] <- c
        }
      }
    }else if( is_above(ratio) ){
      for( j in seq_len(nrow(cand_rc)) ){
        r <- cand_rc[j, "row"]; c <- cand_rc[j, "col"]
        if( !in_grid(r, c) ) next  ## padded-ring trail cell: dead space
        if( decided_mat[r, c] ) next
        if( alive_mat[r, c] ) next
        n_nb <- 0L
        for( delta in list(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L)) ){
          rn <- r + delta[1L]; cn <- c + delta[2L]
          if( rn < 1L || rn > nr || cn < 1L || cn > nc ) next
          if( alive_mat[rn, cn] || pending_alive_mat[rn, cn] ){
            n_nb <- n_nb + 1L
            if( n_nb >= 2L ) break
          }
        }
        if( n_nb >= 2L ){
          n_rev <- n_rev + 1L
          to_revive_r[n_rev] <- r; to_revive_c[n_rev] <- c
          pending_alive_mat[r, c] <- TRUE
        }
      }
    }
  }

  ## Spur pass: perimeter cells only. Interior cells are never candidates
  ## for removal here -- they belong to the area's body and pass through
  ## cleanup unchanged. Range / column-count / row-count statistics are
  ## still derived from the full area for correct denominators.
  if( nrow(area_cells) == 0L ){
    area_row_range <- 0L; area_col_range <- 0L
    col_counts <- integer(nc); row_counts <- integer(nr)
  }else{
    area_row_range <- max(area_cells[, "row"]) - min(area_cells[, "row"]) + 1L
    area_col_range <- max(area_cells[, "col"]) - min(area_cells[, "col"]) + 1L
    col_counts <- tabulate(area_cells[, "col"], nbins = nc)
    row_counts <- tabulate(area_cells[, "row"], nbins = nr)
  }

  for( mi in seq_len(nrow(perim_cells)) ){
    r <- perim_cells[mi, "row"]; c <- perim_cells[mi, "col"]
    if( trail_mat[r, c] ) next
    if( decided_mat[r, c] ) next

    ## Local 3x3 alive ratio, including the candidate cell.
    ## With cleanup thresholds = c(0.2, 0.9), a 1-alive-neighbour cell
    ## scores 2/9 ~= 0.222 (kept) under self-inclusion vs 1/8 = 0.125
    ## (pruned) under self-exclusion. Preserves 1-neighbour spurs and
    ## 1-cell islands -- robust on noisy / marginal alive sets.
    n_total <- 0L; n_alive <- 0L
    for( dr in -1:1 ) for( dc in -1:1 ){
      rn <- r + dr; cn <- c + dc
      if( rn < 1L || rn > nr || cn < 1L || cn > nc ) next
      n_total <- n_total + 1L
      if( alive_mat[rn, cn] ) n_alive <- n_alive + 1L
    }
    local_ratio <- if( n_total > 0L ) n_alive / n_total else 1.0
    spur_ratio <- min(local_ratio,
                      col_counts[c] / area_row_range,
                      row_counts[r] / area_col_range)
    if( is_below(spur_ratio) ){
      n_rem <- n_rem + 1L
      to_remove_r[n_rem] <- r; to_remove_c[n_rem] <- c
    }
  }

  to_revive <- if( n_rev > 0L )
    cbind(row = to_revive_r[seq_len(n_rev)],
          col = to_revive_c[seq_len(n_rev)]) else empty_rc
  to_remove <- if( n_rem > 0L )
    cbind(row = to_remove_r[seq_len(n_rem)],
          col = to_remove_c[seq_len(n_rem)]) else empty_rc

  to_revive <- unique(to_revive)
  to_remove <- unique(to_remove)
  if( nrow(to_revive) > 0L && nrow(to_remove) > 0L ){
    rev_keys <- paste(to_revive[, "row"], to_revive[, "col"], sep = "_")
    rem_keys <- paste(to_remove[, "row"], to_remove[, "col"], sep = "_")
    to_revive <- to_revive[!rev_keys %in% rem_keys, , drop = FALSE]
  }
  return(list(to_revive = to_revive, to_remove = to_remove))
}

# =========================================================================
## Merge separable areas (merge_separate_areas = TRUE)
# =========================================================================

#' Bridge separable areas via 1-pixel dead-space connectors (internal)
#'
#' Identifies dead cells 8-adjacent to two or more areas; tests each
#' candidate against the alive predicate. Revived candidates are
#' written to alive_mat and locked in decided_mat.
#'
#' @noRd
.mergeSeparableAreas <- function(alive_mat, decided_mat, areas_info,
                                 feat_array, mu,
                                 area_sensitivity, area_threshold){

  nr <- nrow(alive_mat); nc <- ncol(alive_mat)
  if( areas_info$n_areas <= 1L ) return(list(alive_mat = alive_mat,
                                             decided_mat = decided_mat,
                                             revived = 0L))
  label_mat <- areas_info$label_mat

  ## Per dead cell: count distinct labels among 8-neighbours.
  dead_idx <- which(!alive_mat & !decided_mat, arr.ind = TRUE)
  if( nrow(dead_idx) == 0L ) return(list(alive_mat = alive_mat,
                                         decided_mat = decided_mat,
                                         revived = 0L))

  nw <- dim(feat_array)[3L]
  thr_lo <- area_threshold[1L]; thr_hi <- area_threshold[2L]
  pixel_alive <- function(r, c){
    n_match <- 0L; any_veg <- FALSE
    for( w in seq_len(nw) ){
      val <- feat_array[r, c, w]
      if( !is.finite(val) || !is.finite(mu[w]) ) next
      if( val >= thr_lo && val <= thr_hi ) any_veg <- TRUE
      if( abs(val - mu[w]) < area_sensitivity ) n_match <- n_match + 1L
    }
    any_veg && n_match >= 1L
  }

  revived <- 0L
  for( i in seq_len(nrow(dead_idx)) ){
    r <- dead_idx[i, 1L]; c <- dead_idx[i, 2L]
    r_lo <- max(1L, r - 1L); r_hi <- min(nr, r + 1L)
    c_lo <- max(1L, c - 1L); c_hi <- min(nc, c + 1L)
    block <- label_mat[r_lo:r_hi, c_lo:c_hi]
    labs <- unique(as.integer(block))
    labs <- labs[labs != 0L]
    if( length(labs) < 2L ) next
    if( pixel_alive(r, c) ){
      alive_mat[r, c]   <- TRUE
      decided_mat[r, c] <- TRUE
      revived <- revived + 1L
    }
  }
  return(list(alive_mat = alive_mat, decided_mat = decided_mat, revived = revived))
}

# =========================================================================
# Per-area perimeter walk helper
# =========================================================================

#' Build a perimeter walk for an area (internal)
#'
#' Returns NULL if the area has too few cells to walk.
#'
#' @noRd
.areaWalk <- function(area, nr, nc){

  cells <- area$cells
  if( nrow(cells) == 0L ) return(NULL)
  perim <- .areaPerimeter(cells)
  if( nrow(perim) < 3L ) return(NULL)
  area_mat <- matrix(FALSE, nr, nc)
  area_mat[cells] <- TRUE
  walk <- .walkPerimeter(area_mat)
  walk$perim_rc   <- perim
  walk$area_cells <- cells
  return(walk)
}
