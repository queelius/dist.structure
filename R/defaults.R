# ==========================================================================
# Topology defaults on `dist_structure`
# ==========================================================================
#
# Each default is a direct composition of other generics. Implementors
# override when a closed-form or specialized algorithm is faster.
# ==========================================================================


#' @rdname phi
#' @export
phi.dist_structure <- function(x, state) {
  paths <- min_paths(x)
  for (P in paths) {
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


# Internal helper: remove non-minimal sets from a list of integer vectors.
minimize_sets <- function(sets) {
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


#' @rdname critical_states
#' @export
critical_states.dist_structure <- function(x, j) {
  m <- ncomponents(x)
  stopifnot(length(j) == 1L, j >= 1L, j <= m)
  other <- setdiff(seq_len(m), j)
  grid <- expand.grid(rep(list(c(0L, 1L)), m - 1L), KEEP.OUT.ATTRS = FALSE)
  grid <- as.matrix(grid)
  colnames(grid) <- NULL
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
  status <- ifelse(times < t_sys, "left",
            ifelse(times > t_sys, "right", "exact"))
  list(system_time = t_sys, component_status = status)
}


#' @rdname is_coherent
#' @export
is_coherent.dist_structure <- function(x) {
  m <- ncomponents(x)
  for (j in seq_len(m)) {
    other <- setdiff(seq_len(m), j)
    grid <- expand.grid(rep(list(c(0L, 1L)), m - 1L), KEEP.OUT.ATTRS = FALSE)
    grid <- as.matrix(grid); colnames(grid) <- NULL
    any_differ <- FALSE
    x0 <- x1 <- integer(m)
    for (i in seq_len(nrow(grid))) {
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


#' @rdname reliability
#' @export
reliability.dist_structure <- function(x, p) {
  m <- ncomponents(x)
  if (length(p) == 1L) p <- rep(p, m)
  stopifnot(length(p) == m, all(p >= 0 & p <= 1))
  grid <- expand.grid(rep(list(c(0L, 1L)), m), KEEP.OUT.ATTRS = FALSE)
  grid <- as.matrix(grid); colnames(grid) <- NULL
  total <- 0
  for (i in seq_len(nrow(grid))) {
    state <- grid[i, ]
    pv <- prod(ifelse(state == 1L, p, 1 - p))
    total <- total + phi(x, state) * pv
  }
  total
}
