# ==========================================================================
# IID convenience constructors
# ==========================================================================
#
# Shortcuts for building dist_structure objects when all components share
# the same distribution.
# ==========================================================================


#' Series system of m iid components
#'
#' Equivalent to the minimum of `m` iid random variables from `d`. Unlike
#' `min` in the base distribution algebra, `min_iid` preserves the series
#' topology so structural queries (phi, min_paths, structural_importance)
#' remain available.
#'
#' @param d A `dist` object (the common component distribution).
#' @param m Number of components (positive integer).
#' @return A `series_dist`.
#' @examples
#' sys <- min_iid(algebraic.dist::exponential(1), m = 3)
#' # System survival at t=0.5 equals (exp(-t))^3 = exp(-1.5)
#' @export
min_iid <- function(d, m) {
  stopifnot(m >= 1L)
  series_dist(replicate(m, d, simplify = FALSE))
}


#' Parallel system of m iid components
#'
#' Equivalent to the maximum of `m` iid random variables from `d`.
#' Preserves the parallel topology for structural queries.
#'
#' @param d A `dist` object.
#' @param m Number of components.
#' @return A `parallel_dist`.
#' @export
max_iid <- function(d, m) {
  stopifnot(m >= 1L)
  parallel_dist(replicate(m, d, simplify = FALSE))
}


#' k-th order statistic of m iid components
#'
#' Constructs a `kofn_dist` whose system lifetime equals `T_(k)`, the
#' k-th order statistic of `m` iid draws from `d`. Under the k-of-m
#' parametrization, the system fails at the `(m - k + 1)`-th component
#' failure; setting the threshold to `m - k + 1` makes the system lifetime
#' equal `T_(k)`.
#'
#' @param d A `dist` object.
#' @param k The order statistic index (1 = min, m = max).
#' @param m Number of iid components.
#' @return A `kofn_dist`.
#' @examples
#' # Median of 5 iid exponentials
#' sys <- order_statistic(algebraic.dist::exponential(1), k = 3, m = 5)
#' @export
order_statistic <- function(d, k, m) {
  stopifnot(k >= 1L, k <= m)
  kofn_dist(m - k + 1L, replicate(m, d, simplify = FALSE))
}
