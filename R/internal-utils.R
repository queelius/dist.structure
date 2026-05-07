# ==========================================================================
# Internal helpers (not exported)
# ==========================================================================
#
# Shared utilities used by the topology defaults in defaults.R and the
# closed-form specializations in exp_*.R, wei_*.R, gamma_series.R, and
# lognormal_series.R. None of these are part of the public API.
# ==========================================================================


# Remove non-minimal sets from a list of integer vectors. Duplicates are
# collapsed first (the strict-subset minimality check would not remove
# equal-length duplicates on its own); this is important in Berge
# transversal iteration, where independent path extensions routinely
# produce repeated transversals.
minimize_sets <- function(sets) {
  sets <- unique(sets)
  n <- length(sets)
  if (n == 0L) return(sets)
  is_min <- rep(TRUE, n)
  for (i in seq_len(n)) {
    if (!is_min[i]) next
    for (j in seq_len(n)) {
      if (i == j || !is_min[j]) next
      if (all(sets[[i]] %in% sets[[j]]) &&
          length(sets[[i]]) < length(sets[[j]])) {
        is_min[j] <- FALSE
      }
    }
  }
  sets[is_min]
}


# Aggregate survival for a k-of-m system given the per-component survival
# probabilities at a single time point. Computes P(at least k of m
# components are alive), summing over all subsets A with |A| >= k of
# (prod S_j for j in A) * (prod (1 - S_j) for j not in A). Shared by
# exp_kofn and wei_kofn surv methods.
kofn_surv_probability <- function(comp_surv, k) {
  m <- length(comp_surv)
  comp_fail <- 1 - comp_surv
  all_idx <- seq_len(m)
  total <- 0
  for (sz in seq.int(k, m)) {
    for (A in utils::combn(m, sz, simplify = FALSE)) {
      alive <- prod(comp_surv[A])
      failed <- setdiff(all_idx, A)
      dead <- if (length(failed) == 0L) 1 else prod(comp_fail[failed])
      total <- total + alive * dead
    }
  }
  total
}


# Aggregate density for a k-of-m system at a single time point, given
# per-component densities f_j(t), survivals S_j(t), and CDFs F_j(t) =
# 1 - S_j(t). Uses the critical-state formula:
#   f_sys(t) = sum_{j=1}^m f_j(t) * P(component j is critical at t)
# where j is critical iff exactly (k - 1) of the other components are
# alive (so j's failure at t drops the alive count from k to k - 1,
# triggering system failure). Shared by exp_kofn and wei_kofn density
# methods.
kofn_density_value <- function(comp_dens, comp_surv, k) {
  m <- length(comp_surv)
  comp_fail <- 1 - comp_surv
  total <- 0
  for (j in seq_len(m)) {
    others <- setdiff(seq_len(m), j)
    p_crit <- 0
    if (k == 1L) {
      # Need 0 others alive.
      p_crit <- prod(comp_fail[others])
    } else if (k - 1L > length(others)) {
      p_crit <- 0
    } else {
      subsets <- utils::combn(others, k - 1L, simplify = FALSE)
      for (B in subsets) {
        alive_term <- if (length(B) == 0L) 1 else prod(comp_surv[B])
        dead_idx <- setdiff(others, B)
        dead_term <- if (length(dead_idx) == 0L) 1 else prod(comp_fail[dead_idx])
        p_crit <- p_crit + alive_term * dead_term
      }
    }
    total <- total + comp_dens[j] * p_crit
  }
  total
}


# Collect n samples from each of m component samplers into an n-by-m
# matrix. Handles the n == 1 / m == 1 edge cases where vapply would
# otherwise drop dimensions.
sample_component_matrix <- function(samplers, n) {
  m <- length(samplers)
  mat <- vapply(samplers, function(samp) samp(n), numeric(n))
  if (!is.matrix(mat)) dim(mat) <- c(n, m)
  mat
}


# Build a list of per-component sampler closures from a base stats
# rXxx function and parallel parameter vectors. The resulting samplers
# each take a single arg `n` and draw n independent values for that
# component. Used by the closed-form series/parallel/kofn constructors.
make_component_samplers <- function(rfun, ...) {
  args <- list(...)
  m <- length(args[[1L]])
  lapply(seq_len(m), function(j) {
    per_args <- lapply(args, `[[`, j)
    function(n) do.call(rfun, c(list(n = n), per_args))
  })
}


# Product-of-per-component-survivals closure for a series system where
# each component's survival at time ti has the form
# `pXxx(ti, ..., lower.tail = FALSE)`. `params` is a named list of
# parameter vectors (each of length m) to pass to `pfun`. Returns a
# function(t) that vectorises over `t`.
series_surv_product <- function(pfun, params) {
  function(t, ...) {
    vapply(t, function(ti) {
      args <- c(list(q = ti), params, list(lower.tail = FALSE))
      prod(do.call(pfun, args))
    }, numeric(1L))
  }
}


# 2^n binary grid as an integer matrix (no col names, no expand.grid
# attributes). Rows enumerate states in {0, 1}^n.
binary_grid <- function(n) {
  if (n == 0L) return(matrix(integer(0), nrow = 1L, ncol = 0L))
  grid <- expand.grid(rep(list(c(0L, 1L)), n), KEEP.OUT.ATTRS = FALSE)
  grid <- as.matrix(grid)
  colnames(grid) <- NULL
  grid
}


# All permutations of 1:n as a list of integer vectors. Base R has no
# such helper; this recursive implementation is adequate for the small m
# values at which enumeration-based signature computation is feasible.
permutations <- function(n) {
  if (n <= 1L) return(list(seq_len(n)))
  smaller <- permutations(n - 1L)
  result <- vector("list", n * length(smaller))
  idx <- 1L
  for (i in seq_len(n)) {
    rest <- seq_len(n)[-i]
    for (p in smaller) {
      result[[idx]] <- c(i, rest[p])
      idx <- idx + 1L
    }
  }
  result
}
