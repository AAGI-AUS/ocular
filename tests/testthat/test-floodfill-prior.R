# Soft-prior admission penalty: mechanism contract (.floodFill) and the
# per-rung slice + default no-op through the widen ladder (.searchAreaWithWidens).
# Pure base R -- no terra/sf needed (the prior is a logical matrix; the
# polygon -> matrix rasterise that does need terra is covered in test-ftw-core.R).

# Synthetic 21x21 scene shared by the .floodFill block: a solid 0.5 field on
# a 0.0 background with a non-veg hole, seed left of the hole inside the field.
.mk_scene <- function() {
  nr <- 21L; nc <- 21L; nw <- 3L
  feat <- array(0.0, dim = c(nr, nc, nw))
  for (w in 1:nw) feat[5:17, 5:17, w] <- 0.5
  for (w in 1:nw) feat[10:12, 13:15, w] <- 0.0   # non-veg hole inside
  list(feat = feat, nr = nr, nc = nc,
       mu = c(0.5, 0.5, 0.5), thr = c(0.3, 0.9), tol = 0.05,
       centre = c(11L, 7L))
}

test_that(".floodFill soft-prior penalty honours its admission contract", {
  s  <- .mk_scene()
  ff <- function(prior = NULL, strength = 0, cross = 0.75)
    ocular:::.floodFill(s$feat, s$centre, s$mu, s$tol, s$thr, 2L,
                        c(1L, s$nr), c(1L, s$nc), s$nr, s$nc,
                        local_alive_density = 0.50,
                        prior_field = prior, prior_strength = strength,
                        prior_cross = cross)

  base <- ff()
  expect_gt(sum(base$alive_mat), 50L)                       # (1) seed fills a field

  # Prior marking the 0.5 pixels adjacent to the hole as 'outside'.
  outside <- matrix(TRUE, s$nr, s$nc)
  for (r in 9:13) for (cc in 12:16)
    if (s$feat[r, cc, 1] == 0.5) outside[r, cc] <- FALSE

  expect_identical(ff(outside, 0, 0.75)$alive_mat, base$alive_mat)   # (2) strength 0 -> no-op
  expect_identical(ff(matrix(TRUE, s$nr, s$nc), 1, 1.01)$alive_mat,
                   base$alive_mat)                                    # (3) all-inside -> no-op

  pen <- ff(outside, 1, 1.01)                                        # (4) uncrossable penalty bites
  expect_lt(sum(pen$alive_mat), sum(base$alive_mat))
  expect_true(all(which(pen$alive_mat) %in% which(base$alive_mat)))  # (5) penalty only removes (monotone)

  expect_identical(ff(outside, 1, 0.0)$alive_mat, base$alive_mat)    # (6) evidence crossing restores baseline

  na_prior <- matrix(NA, s$nr, s$nc); na_prior[5:17, 5:17] <- TRUE
  expect_identical(ff(na_prior, 1, 1.01)$alive_mat, base$alive_mat)  # (7) NA entries impose no penalty

  expect_error(ff(matrix(TRUE, 3, 3), 1),                            # (8) dimension guard
               "prior_field must be an")
})

test_that(".searchAreaWithWidens per-rung prior slice is dimensionally sound", {
  # nr=nc=30 -> sqrt(2) ladder rungs 14,21,30 (>=2 descents); the largest
  # rung is the full window (== the .floodFill block above), so the NON-
  # trivial coverage is the clamped sub-windows. All-veg feat so every
  # rung's fill reaches its boundary and the loop does NOT early-return,
  # exercising the slice at every rung. proto_prior.R never hits this path.
  nr <- 30L; nc <- 30L; nw <- 2L
  feat <- array(0.5, dim = c(nr, nc, nw)); mu <- rep(0.5, nw)
  saw <- function(centre, prior, strength, cross)
    ocular:::.searchAreaWithWidens(
      feat_array = feat, centre_rc = centre, mu = mu,
      area_sensitivity = 0.05, area_threshold = c(0.3, 0.9),
      min_windows_alive = 1L, local_alive_density = 0.50,
      nr_full = nr, nc_full = nc,
      prior_field = prior, prior_strength = strength, prior_cross = cross)

  full_prior <- matrix(TRUE, nr, nc)
  seeds <- list(centre = c(15L, 15L), corner = c(2L, 2L),
                edge = c(2L, 15L), offcentre = c(6L, 25L))
  for (nm in names(seeds)) {
    res <- saw(seeds[[nm]], full_prior, 1, 1.01)          # uncrossable, slice entered each rung
    expect_true(is.matrix(res) && all(dim(res) == c(nr, nc)),
                info = paste("seed", nm))
  }

  # Interior-FALSE prior, off-centre seed (asymmetric clamps): penalised set
  # is a strict subset of the no-prior set -> penalty runs in sub-window
  # coordinates and is monotone.
  hole <- matrix(TRUE, nr, nc); hole[10:20, 10:20] <- FALSE
  base_np <- saw(c(6L, 25L), NULL, 0, 1.01)
  pen     <- saw(c(6L, 25L), hole, 1, 1.01)
  expect_lt(sum(pen), sum(base_np))
  expect_true(all(which(pen) %in% which(base_np)))
})

test_that(".searchAreaWithWidens soft prior is inert on the default path", {
  # The exact no-op argument states the wiring sends when no ftw_prior is
  # attached: prior_field NULL OR strength 0. Inert through the whole widen
  # ladder, not just .floodFill at full window.
  nr <- 30L; nc <- 30L; nw <- 2L
  feat <- array(0.5, dim = c(nr, nc, nw)); mu <- rep(0.5, nw)
  hole <- matrix(TRUE, nr, nc); hole[10:20, 10:20] <- FALSE
  saw <- function(centre, prior, strength)
    ocular:::.searchAreaWithWidens(
      feat_array = feat, centre_rc = centre, mu = mu,
      area_sensitivity = 0.05, area_threshold = c(0.3, 0.9),
      min_windows_alive = 1L, local_alive_density = 0.50,
      nr_full = nr, nc_full = nc,
      prior_field = prior, prior_strength = strength, prior_cross = 0.75)

  for (centre in list(c(15L, 15L), c(2L, 2L), c(6L, 25L))) {
    a  <- saw(centre, NULL, 0)                 # default
    b  <- saw(centre, matrix(TRUE, nr, nc), 0) # field present, strength 0
    c_ <- saw(centre, hole, 0)                 # FALSE field, strength 0
    expect_identical(a, b)
    expect_identical(a, c_)
  }
})
