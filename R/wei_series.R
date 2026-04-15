# ==========================================================================
# Closed-form series of Weibull components
# ==========================================================================
#
# For a series of m independent Weibulls with shapes k_j and scales s_j:
#   Cumulative hazard of component j at time t: H_j(t) = (t / s_j)^{k_j}
#   System survival:  S_sys(t) = exp(-sum_j (t / s_j)^{k_j})
#   System hazard:    h_sys(t) = sum_j (k_j / s_j) * (t / s_j)^{k_j - 1}
#   System density:   f_sys(t) = h_sys(t) * S_sys(t)
#
# No closed form for the CDF inverse in general; sampler uses min of
# independently drawn Weibulls.
# ==========================================================================


#' Series of heterogeneous Weibull components (closed form)
#'
#' Constructs a `dist_structure` representing a series system whose
#' components are independent Weibull distributions with possibly
#' different shapes and scales. Closed-form methods are provided for
#' `surv`, `cdf`, `sampler`, and `algebraic.dist::hazard`.
#'
#' @param shapes Positive numeric vector of length `m`: Weibull shape
#'   parameters per component.
#' @param scales Positive numeric vector of length `m` (same length as
#'   `shapes`): Weibull scale parameters per component.
#' @return An object of class
#'   `c("wei_series", "series_dist", "coherent_dist", "dist_structure",
#'   "univariate_dist", "continuous_dist", "dist")`.
#' @examples
#' sys <- wei_series(shapes = c(1, 2, 3), scales = c(1, 2, 3))
#' algebraic.dist::surv(sys)(1)
#' @export
wei_series <- function(shapes, scales) {
  stopifnot(length(shapes) == length(scales),
            all(shapes > 0), all(scales > 0))
  m <- length(shapes)
  components <- lapply(seq_len(m), function(j) {
    algebraic.dist::weibull_dist(shape = shapes[j], scale = scales[j])
  })
  obj <- series_dist(components)
  obj$shapes <- as.numeric(shapes)
  obj$scales <- as.numeric(scales)
  class(obj) <- c("wei_series", class(obj))
  obj
}


#' @rdname wei_series
#' @param x A `wei_series` object.
#' @param ... Ignored.
#' @export
surv.wei_series <- function(x, ...) {
  shapes <- x$shapes
  scales <- x$scales
  function(t, ...) {
    vapply(t, function(ti) {
      exp(-sum((ti / scales)^shapes))
    }, numeric(1L))
  }
}


#' @rdname wei_series
#' @export
cdf.wei_series <- function(x, ...) {
  S <- surv.wei_series(x)
  function(t, ...) 1 - S(t)
}


#' @rdname wei_series
#' @export
sampler.wei_series <- function(x, ...) {
  shapes <- x$shapes
  scales <- x$scales
  m <- length(shapes)
  function(n, ...) {
    mat <- vapply(seq_len(m), function(j) {
      stats::rweibull(n, shape = shapes[j], scale = scales[j])
    }, numeric(n))
    if (!is.matrix(mat)) mat <- matrix(mat, nrow = n)
    apply(mat, 1L, min)
  }
}
