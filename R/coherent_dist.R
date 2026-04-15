# ==========================================================================
# General coherent_dist + topology shortcut constructors
# ==========================================================================


#' Coherent system distribution from minimal path sets
#'
#' General-purpose constructor. Users supply a list of minimal path sets
#' (each a vector of component indices) and a list of component
#' distributions (each an `algebraic.dist::dist` object with parameters
#' baked in). The resulting object is a `dist_structure` and `dist`.
#'
#' @param min_paths List of integer vectors; each is a minimal path set.
#' @param components List of `dist` objects, length `m`. Each is a
#'   fully-parameterized component lifetime distribution.
#' @param m Optional integer. Inferred from `components` and `min_paths`
#'   if omitted.
#' @return An object of class
#'   `c("coherent_dist", "dist_structure", "univariate_dist", "dist")`.
#' @examples
#' # A bridge network with exponential components
#' sys <- coherent_dist(
#'   min_paths = list(c(1, 4), c(2, 5), c(1, 3, 5), c(2, 3, 4)),
#'   components = replicate(5, algebraic.dist::exponential(1), simplify = FALSE)
#' )
#' reliability(sys, 0.9)
#' @export
coherent_dist <- function(min_paths, components, m = NULL) {
  stopifnot(is.list(min_paths), is.list(components))
  if (is.null(m)) {
    all_idx <- unlist(min_paths)
    m <- max(length(components),
             if (length(all_idx)) max(all_idx) else 0L)
  }
  stopifnot(length(components) == m)
  structure(
    list(
      min_paths = lapply(min_paths, as.integer),
      m = as.integer(m),
      components = components
    ),
    class = c("coherent_dist", "dist_structure",
              "univariate_dist", "continuous_dist", "dist")
  )
}


#' @export
ncomponents.coherent_dist <- function(x) x$m


#' @export
component.coherent_dist <- function(x, j, ...) {
  stopifnot(j >= 1L, j <= x$m)
  x$components[[j]]
}


#' @export
min_paths.coherent_dist <- function(x) x$min_paths


#' Dual of a coherent_dist: swap cuts and paths
#'
#' Overrides the lazy-wrapper default with a proper `coherent_dist`:
#' `min_paths(dual(x)) = min_cuts(x)`.
#'
#' @rdname dual
#' @export
dual.coherent_dist <- function(x) {
  coherent_dist(
    min_paths = min_cuts(x),
    components = x$components,
    m = x$m
  )
}


# ==========================================================================
# Topology shortcut constructors
# ==========================================================================


#' Series system distribution
#'
#' A series system fails if any single component fails. Equivalent to
#' `min(components)` as random variables; unlike `min` in the base algebra,
#' `series_dist` preserves the component decomposition so topology queries
#' work.
#'
#' @param components List of `dist` objects.
#' @return A `series_dist` inheriting from `coherent_dist`.
#' @examples
#' sys <- series_dist(replicate(3, algebraic.dist::exponential(1), simplify = FALSE))
#' algebraic.dist::surv(sys)(0.5)
#' @export
series_dist <- function(components) {
  stopifnot(is.list(components), length(components) >= 1L)
  m <- length(components)
  obj <- coherent_dist(
    min_paths = list(seq_len(m)),
    components = components,
    m = m
  )
  class(obj) <- c("series_dist", class(obj))
  obj
}


#' @export
phi.series_dist <- function(x, state) as.integer(all(state == 1L))

#' @export
min_paths.series_dist <- function(x) list(seq_len(x$m))


#' Parallel system distribution
#'
#' A parallel system fails only when all components fail. Equivalent to
#' `max(components)` but preserves topology.
#'
#' @param components List of `dist` objects.
#' @return A `parallel_dist` inheriting from `coherent_dist`.
#' @export
parallel_dist <- function(components) {
  stopifnot(is.list(components), length(components) >= 1L)
  m <- length(components)
  obj <- coherent_dist(
    min_paths = lapply(seq_len(m), function(j) j),
    components = components,
    m = m
  )
  class(obj) <- c("parallel_dist", class(obj))
  obj
}


#' @export
phi.parallel_dist <- function(x, state) as.integer(any(state == 1L))

#' @export
min_paths.parallel_dist <- function(x) lapply(seq_len(x$m), as.integer)


#' k-out-of-n system distribution
#'
#' A k-out-of-n system functions if at least `k` of its `m` components
#' function. Equivalent to the `(m - k + 1)`-th order statistic of
#' component lifetimes.
#'
#' @param k Minimum functioning components for system operation.
#' @param components List of `dist` objects (length `m`).
#' @return A `kofn_dist` inheriting from `coherent_dist`.
#' @export
kofn_dist <- function(k, components) {
  stopifnot(is.list(components), length(components) >= 1L)
  m <- length(components)
  stopifnot(k >= 1L, k <= m)
  paths <- utils::combn(m, k, simplify = FALSE)
  obj <- coherent_dist(
    min_paths = paths,
    components = components,
    m = m
  )
  obj$k <- as.integer(k)
  class(obj) <- c("kofn_dist", class(obj))
  obj
}


#' @export
phi.kofn_dist <- function(x, state) as.integer(sum(state) >= x$k)

#' @export
min_paths.kofn_dist <- function(x) {
  lapply(utils::combn(x$m, x$k, simplify = FALSE), as.integer)
}


#' Bridge system distribution
#'
#' The classical 5-component bridge reliability network with minimal
#' path sets `{1,4}`, `{2,5}`, `{1,3,5}`, `{2,3,4}`.
#'
#' @param components List of 5 `dist` objects.
#' @return A `bridge_dist` inheriting from `coherent_dist`.
#' @export
bridge_dist <- function(components) {
  stopifnot(is.list(components), length(components) == 5L)
  paths <- list(
    c(1L, 4L), c(2L, 5L),
    c(1L, 3L, 5L), c(2L, 3L, 4L)
  )
  obj <- coherent_dist(
    min_paths = paths,
    components = components,
    m = 5L
  )
  class(obj) <- c("bridge_dist", class(obj))
  obj
}


#' Consecutive-k-out-of-n system distribution
#'
#' The system fails if any `k` consecutive components fail. The structure
#' functions when at least one block of `k` consecutive components all
#' function; minimal path sets are the `n - k + 1` consecutive blocks of
#' size `k`.
#'
#' @param k Block size.
#' @param components List of `dist` objects (length `n`).
#' @return A `consecutive_k_dist` inheriting from `coherent_dist`.
#' @export
consecutive_k_dist <- function(k, components) {
  n <- length(components)
  stopifnot(k >= 1L, k <= n)
  paths <- lapply(seq_len(n - k + 1L), function(i) as.integer(i:(i + k - 1L)))
  obj <- coherent_dist(
    min_paths = paths,
    components = components,
    m = n
  )
  class(obj) <- c("consecutive_k_dist", class(obj))
  obj
}
