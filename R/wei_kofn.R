# ==========================================================================
# Closed-form k-out-of-n of heterogeneous Weibull components
# ==========================================================================
#
# Same subset-enumeration formula as exp_kofn, with Weibull per-component
# survivals S_j(t) = exp(-(t / scale_j)^shape_j).
# ==========================================================================


#' k-out-of-n system of independent Weibull components (closed form)
#'
#' Constructs a `dist_structure` for a k-out-of-m system whose components
#' are independent (possibly heterogeneous) Weibulls. Closed-form `surv`,
#' `cdf`, and `sampler` via subset enumeration and component order
#' statistics.
#'
#' @param k Minimum functioning components for system operation.
#' @param shapes Positive numeric vector of length `m`.
#' @param scales Positive numeric vector of length `m`.
#' @return An object of class
#'   `c("wei_kofn", "kofn_dist", "coherent_dist", "dist_structure",
#'   "univariate_dist", "continuous_dist", "dist")`.
#' @examples
#' sys <- wei_kofn(k = 2, shapes = c(1, 2, 3), scales = c(1, 2, 3))
#' algebraic.dist::surv(sys)(1)
#' @export
wei_kofn <- function(k, shapes, scales) {
  stopifnot(length(shapes) == length(scales),
            all(shapes > 0), all(scales > 0))
  m <- length(shapes)
  stopifnot(k >= 1L, k <= m)
  components <- lapply(seq_len(m), function(j) {
    algebraic.dist::weibull_dist(shape = shapes[j], scale = scales[j])
  })
  obj <- kofn_dist(k, components)
  obj$shapes <- as.numeric(shapes)
  obj$scales <- as.numeric(scales)
  class(obj) <- c("wei_kofn", class(obj))
  obj
}


#' @rdname wei_kofn
#' @param x A `wei_kofn` object.
#' @param ... Ignored.
#' @export
surv.wei_kofn <- function(x, ...) {
  shapes <- x$shapes
  scales <- x$scales
  k <- x$k
  function(t, ...) {
    vapply(t, function(ti) {
      kofn_surv_probability(exp(-(ti / scales)^shapes), k)
    }, numeric(1L))
  }
}


#' @rdname wei_kofn
#' @export
cdf.wei_kofn <- function(x, ...) {
  S <- surv.wei_kofn(x)
  function(t, ...) 1 - S(t)
}


#' @rdname wei_kofn
#' @export
sampler.wei_kofn <- function(x, ...) {
  order_idx <- length(x$shapes) - x$k + 1L
  samplers <- make_component_samplers(stats::rweibull,
                                      shape = x$shapes, scale = x$scales)
  function(n, ...) {
    apply(sample_component_matrix(samplers, n), 1L,
          function(row) sort(row)[order_idx])
  }
}


#' @rdname wei_kofn
#' @importFrom stats dweibull
#' @export
density.wei_kofn <- function(x, ...) {
  shapes <- x$shapes
  scales <- x$scales
  k <- x$k
  m <- length(shapes)
  function(t, log = FALSE, ...) {
    vals <- vapply(t, function(ti) {
      surv <- exp(-(ti / scales)^shapes)
      dens <- vapply(seq_len(m), function(j) {
        stats::dweibull(ti, shape = shapes[j], scale = scales[j])
      }, numeric(1L))
      kofn_density_value(dens, surv, k)
    }, numeric(1L))
    if (isTRUE(log)) log(vals) else vals
  }
}
