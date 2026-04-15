# ==========================================================================
# Topology defaults on `dist_structure`
# ==========================================================================
#
# Each default is a direct composition of other generics. Implementors
# override when a closed-form or specialized algorithm is faster.
# ==========================================================================


# --------------------------------------------------------------------------
# Internal helpers (not exported)
# --------------------------------------------------------------------------

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
  for (sz in k:m) {
    for (A in utils::combn(m, sz, simplify = FALSE)) {
      alive <- prod(comp_surv[A])
      failed <- setdiff(all_idx, A)
      dead <- if (length(failed) == 0L) 1 else prod(comp_fail[failed])
      total <- total + alive * dead
    }
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
    rest <- if (i == 1L) 2:n else if (i == n) 1:(n - 1L) else c(1:(i - 1L), (i + 1L):n)
    for (p in smaller) {
      result[[idx]] <- c(i, rest[p])
      idx <- idx + 1L
    }
  }
  result
}


# --------------------------------------------------------------------------
# S3 methods
# --------------------------------------------------------------------------


#' @rdname phi
#' @export
phi.dist_structure <- function(x, state) {
  for (P in min_paths(x)) {
    if (all(state[P] == 1L)) return(1L)
  }
  0L
}


#' @rdname min_cuts
#' @export
min_cuts.dist_structure <- function(x) {
  # Berge transversal: iteratively build minimal hitting sets of min_paths.
  paths <- min_paths(x)
  if (length(paths) == 0L) return(list())
  transversals <- list(integer(0))
  for (P in paths) {
    new_trans <- list()
    for (T in transversals) {
      for (p in P) {
        new_trans[[length(new_trans) + 1L]] <- sort(unique(c(T, p)))
      }
    }
    transversals <- minimize_sets(new_trans)
  }
  transversals
}


#' @rdname critical_states
#' @export
critical_states.dist_structure <- function(x, j) {
  m <- ncomponents(x)
  stopifnot(length(j) == 1L, j >= 1L, j <= m)
  other <- setdiff(seq_len(m), j)
  grid <- binary_grid(m - 1L)
  keep <- logical(nrow(grid))
  x1 <- x0 <- integer(m)
  for (i in seq_len(nrow(grid))) {
    s <- grid[i, ]
    x1[other] <- s; x1[j] <- 1L
    x0[other] <- s; x0[j] <- 0L
    keep[i] <- (phi(x, x1) == 1L) && (phi(x, x0) == 0L)
  }
  grid[keep, , drop = FALSE]
}


#' @rdname system_lifetime
#' @export
system_lifetime.dist_structure <- function(x, times) {
  m <- ncomponents(x)
  stopifnot(length(times) == m, all(times >= 0))
  ord <- order(times)
  sorted <- times[ord]
  state <- rep(1L, m)
  for (k in seq_len(m)) {
    state[ord[k]] <- 0L
    if (phi(x, state) == 0L) return(sorted[k])
  }
  sorted[m]
}


#' @rdname system_censoring
#' @export
system_censoring.dist_structure <- function(x, times) {
  t_sys <- system_lifetime(x, times)
  status <- rep("exact", length(times))
  status[times < t_sys] <- "left"
  status[times > t_sys] <- "right"
  list(system_time = t_sys, component_status = status)
}


#' @rdname is_coherent
#' @export
is_coherent.dist_structure <- function(x) {
  m <- ncomponents(x)
  grid <- binary_grid(m - 1L)
  rows <- seq_len(nrow(grid))
  for (j in seq_len(m)) {
    other <- setdiff(seq_len(m), j)
    any_differ <- FALSE
    x0 <- x1 <- integer(m)
    for (i in rows) {
      s <- grid[i, ]
      x0[other] <- s; x0[j] <- 0L
      x1[other] <- s; x1[j] <- 1L
      p0 <- phi(x, x0); p1 <- phi(x, x1)
      if (p1 < p0) return(FALSE)
      if (p1 != p0) any_differ <- TRUE
    }
    if (!any_differ) return(FALSE)
  }
  TRUE
}


#' @rdname structural_importance
#' @export
structural_importance.dist_structure <- function(x, j) {
  m <- ncomponents(x)
  crit <- critical_states(x, j)
  nrow(crit) / (2^(m - 1L))
}


#' @rdname system_signature
#' @export
system_signature.dist_structure <- function(x) {
  m <- ncomponents(x)
  if (m > 9L) {
    warning("system_signature via permutation enumeration is ",
            "expensive for m > 9 (m! grows quickly); override with a ",
            "specialized method when possible.")
  }
  sig_counts <- integer(m)
  for (perm in permutations(m)) {
    state <- rep(1L, m)
    for (i in seq_len(m)) {
      state[perm[[i]]] <- 0L
      if (phi(x, state) == 0L) {
        sig_counts[i] <- sig_counts[i] + 1L
        break
      }
    }
  }
  sig_counts / factorial(m)
}


#' @rdname reliability
#' @export
reliability.dist_structure <- function(x, p) {
  m <- ncomponents(x)
  if (length(p) == 1L) p <- rep(p, m)
  stopifnot(length(p) == m, all(p >= 0 & p <= 1))
  grid <- binary_grid(m)
  # Probability of each 2^m component state vector: state j contributes
  # p[j] if state[j] == 1 else (1 - p[j]). Vectorized across rows.
  prob_matrix <- grid * rep(p, each = nrow(grid)) +
    (1L - grid) * rep(1 - p, each = nrow(grid))
  state_probs <- apply(prob_matrix, 1L, prod)
  phi_values <- apply(grid, 1L, function(s) phi(x, s))
  sum(phi_values * state_probs)
}
