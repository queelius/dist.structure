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
#' @return
#' `wei_series()` returns an object of class
#'   `c("wei_series", "series_dist", "coherent_dist", "dist_structure",
#'   "univariate_dist", "continuous_dist", "dist")`.
#'
#' The associated S3 methods return:
#' - `surv()`, `cdf()`, `hazard()`: a closure `function(t, ...)`.
#' - `sampler()`: a closure `function(n, ...)` returning `n` random
#'   variates from the system lifetime distribution.
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
  # exp(-sum((t/scale_j)^shape_j)) = prod_j exp(-(t/scale_j)^shape_j)
  # = prod_j pweibull(t, shape_j, scale_j, lower.tail = FALSE).
  # Using series_surv_product matches the gamma_series and lognormal_series
  # implementations.
  series_surv_product(stats::pweibull,
                      list(shape = x$shapes, scale = x$scales))
}


#' @rdname wei_series
#' @export
sampler.wei_series <- function(x, ...) {
  samplers <- make_component_samplers(stats::rweibull,
                                      shape = x$shapes, scale = x$scales)
  function(n, ...) {
    apply(sample_component_matrix(samplers, n), 1L, min)
  }
}


# Closed-form Weibull hazard: h_sys(t) = sum_j (k_j / s_j) * (t / s_j)^(k_j - 1).
# Series hazards are additive, and each component's Weibull hazard has a
# direct algebraic form, so the system hazard avoids both numerical
# differentiation and the per-component dispatch cost of the algebraic.dist
# fallback.
#' @rdname wei_series
#' @method hazard wei_series
#' @importFrom algebraic.dist hazard
#' @export
hazard.wei_series <- function(x, ...) {
  shapes <- x$shapes
  scales <- x$scales
  function(t, ...) {
    vapply(t, function(ti) {
      sum((shapes / scales) * (ti / scales)^(shapes - 1))
    }, numeric(1L))
  }
}
