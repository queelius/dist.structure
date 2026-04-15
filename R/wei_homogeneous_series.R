# ==========================================================================
# Closed-form series of Weibulls with common shape
# ==========================================================================
#
# When all components share shape k, the system cumulative hazard is:
#   sum_j (t / s_j)^k  =  t^k * sum_j (1 / s_j^k)
# Let s_agg = ( sum_j (1 / s_j^k) )^{-1/k}
# Then the series lifetime is itself Weibull(shape = k, scale = s_agg).
#
# This identity means wei_homogeneous_series reduces exactly to a single
# Weibull. The constructor preserves the original component information
# (so component(x, j) and topology queries still work) while exposing the
# aggregate Weibull for dist-level queries.
# ==========================================================================


#' Series of Weibull components with common shape (closed form as a Weibull)
#'
#' Constructs a `dist_structure` for a series of Weibull components sharing
#' a common shape parameter. By the standard identity, the system lifetime
#' is itself Weibull with the common shape and an aggregate scale
#' `(sum(1 / scale^shape))^(-1 / shape)`. Methods for `surv`, `cdf`,
#' `sampler`, and `mean` forward to this aggregate Weibull, giving exact
#' closed-form values.
#'
#' @param shape Positive scalar: common Weibull shape.
#' @param scales Positive numeric vector: per-component Weibull scales.
#' @return An object of class
#'   `c("wei_homogeneous_series", "wei_series", "series_dist",
#'   "coherent_dist", "dist_structure", "univariate_dist",
#'   "continuous_dist", "dist")`.
#' @examples
#' sys <- wei_homogeneous_series(shape = 2, scales = c(1, 2, 3))
#' # System lifetime is Weibull(shape = 2, scale = aggregate_scale)
#' algebraic.dist::surv(sys)(1)
#' @export
wei_homogeneous_series <- function(shape, scales) {
  stopifnot(length(shape) == 1L, shape > 0, all(scales > 0))
  shapes <- rep(shape, length(scales))
  obj <- wei_series(shapes = shapes, scales = scales)
  obj$shape <- as.numeric(shape)
  obj$aggregate_scale <- (sum(1 / scales^shape))^(-1 / shape)
  class(obj) <- c("wei_homogeneous_series", class(obj))
  obj
}


#' @rdname wei_homogeneous_series
#' @param x A `wei_homogeneous_series` object.
#' @param ... Ignored.
#' @export
surv.wei_homogeneous_series <- function(x, ...) {
  k <- x$shape
  s_agg <- x$aggregate_scale
  function(t, ...) exp(-(t / s_agg)^k)
}


#' @rdname wei_homogeneous_series
#' @export
cdf.wei_homogeneous_series <- function(x, ...) {
  k <- x$shape
  s_agg <- x$aggregate_scale
  function(t, ...) 1 - exp(-(t / s_agg)^k)
}


#' @rdname wei_homogeneous_series
#' @export
sampler.wei_homogeneous_series <- function(x, ...) {
  k <- x$shape
  s_agg <- x$aggregate_scale
  function(n, ...) stats::rweibull(n, shape = k, scale = s_agg)
}


#' @rdname wei_homogeneous_series
#' @export
mean.wei_homogeneous_series <- function(x, ...) {
  x$aggregate_scale * gamma(1 + 1 / x$shape)
}
