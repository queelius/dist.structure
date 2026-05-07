# ==========================================================================
# Topology defaults on `dist_structure`
# ==========================================================================
#
# Each default is a direct composition of other generics. Implementors
# override when a closed-form or specialized algorithm is faster.
# Internal helpers used here live in R/internal-utils.R.
# ==========================================================================


#' @rdname phi
#' @export
phi.dist_structure <- function(x, state) {
  for (P in min_paths(x)) {
    if (all(state[P] == 1L)) return(1L)
  }
  0L
}


# Default min_paths derived by enumeration from the structure function.
# Implementors who provide only phi.X (no min_paths.X) get this default;
# implementors who provide min_paths.X override it. The enumeration is
# O(2^m) and therefore most appropriate for small to moderate m or as
# a fallback during prototyping; specialized subclasses with closed-form
# minimal paths should provide their own min_paths method.
#' @rdname min_paths
#' @export
min_paths.dist_structure <- function(x) {
  m <- ncomponents(x)
  grid <- binary_grid(m)
  good <- vapply(seq_len(nrow(grid)),
                 function(i) phi(x, grid[i, ]), integer(1L))
  active <- which(good == 1L)
  if (length(active) == 0L) return(list())
  candidate_sets <- lapply(active, function(i) which(grid[i, ] == 1L))
  minimize_sets(candidate_sets)
}


#' @rdname min_cuts
#' @export
min_cuts.dist_structure <- function(x) {
  # Berge transversal: iteratively build minimal hitting sets of min_paths.
  paths <- min_paths(x)
  if (length(paths) == 0L) return(list())
  transversals <- list(integer(0))
  for (P in paths) {
    # Pre-size: every (transversal, path-element) pair generates one
    # candidate. Avoids the O(n^2) cost of list[[length(.) + 1L]] <- x.
    new_trans <- vector("list", length(transversals) * length(P))
    idx <- 1L
    for (tr in transversals) {
      for (p in P) {
        new_trans[[idx]] <- sort(unique(c(tr, p)))
        idx <- idx + 1L
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
  # Signature is undefined for non-functioning systems (no min_paths
  # means phi is identically 0). Without this guard, the algorithm
  # below silently returns (1, 0, ..., 0) for any always-zero phi.
  if (length(min_paths(x)) == 0L || phi(x, rep(1L, m)) == 0L) {
    stop("system_signature is undefined for systems with no minimal ",
         "path sets (phi identically zero).", call. = FALSE)
  }
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
