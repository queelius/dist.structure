# ==========================================================================
# Closed-form k-out-of-n of heterogeneous exponential components
# ==========================================================================
#
# For k-of-m of heterogeneous exponentials:
#   S_sys(t) = P(at least k of m exponentials are alive at t)
#            = sum_{A subset of [m], |A|>=k} prod_{j in A} exp(-rate_j t)
#              * prod_{j notin A} (1 - exp(-rate_j t))
# Sampler uses the (m - k + 1)-th order statistic of independently sampled
# component lifetimes.
# ==========================================================================


#' k-out-of-n system of independent exponential components (closed form)
#'
#' Constructs a `dist_structure` for a k-out-of-m system whose components
#' are independent exponentials. Closed-form methods are provided for
#' `surv`, `cdf`, `sampler`, `density`, and `hazard`. `mean` falls back
#' to numerical integration via the `dist_structure` default.
#'
#' @param k Minimum number of functioning components for system operation.
#' @param rates Positive numeric vector of length `m` with `m >= k`.
#' @return
#' `exp_kofn()` returns an object of class
#'   `c("exp_kofn", "kofn_dist", "coherent_dist", "dist_structure",
#'   "univariate_dist", "continuous_dist", "dist")`.
#'
#' The associated S3 methods return:
#' - `surv()`, `density()`, `hazard()`: a closure `function(t, ...)`.
#' - `cdf()` is derived via the `dist_structure` default and returns
#'   a closure `function(t, ...)` equal to `1 - surv(x)(t)`.
#' - `sampler()`: a closure `function(n, ...)` returning `n` random
#'   variates from the system lifetime distribution.
#' @examples
#' sys <- exp_kofn(k = 2, rates = c(1, 2, 3))
#' algebraic.dist::surv(sys)(1)
#' @export
exp_kofn <- function(k, rates) {
  stopifnot(is.numeric(rates), length(rates) >= 1L, all(rates > 0))
  m <- length(rates)
  stopifnot(k >= 1L, k <= m)
  components <- lapply(rates, algebraic.dist::exponential)
  obj <- kofn_dist(k, components)
  obj$rates <- as.numeric(rates)
  class(obj) <- c("exp_kofn", class(obj))
  obj
}


#' @rdname exp_kofn
#' @param x An `exp_kofn` object.
#' @param ... Ignored.
#' @export
surv.exp_kofn <- function(x, ...) {
  rates <- x$rates
  k <- x$k
  function(t, ...) {
    vapply(t, function(ti) {
      kofn_surv_probability(exp(-rates * ti), k)
    }, numeric(1L))
  }
}


#' @rdname exp_kofn
#' @export
sampler.exp_kofn <- function(x, ...) {
  order_idx <- length(x$rates) - x$k + 1L
  samplers <- make_component_samplers(stats::rexp, rate = x$rates)
  function(n, ...) {
    apply(sample_component_matrix(samplers, n), 1L,
          function(row) sort(row)[order_idx])
  }
}


#' @rdname exp_kofn
#' @method hazard exp_kofn
#' @importFrom algebraic.dist hazard
#' @export
hazard.exp_kofn <- function(x, ...) {
  # h_sys(t) = f_sys(t) / S_sys(t); both factors have closed forms here.
  f_fn <- density.exp_kofn(x)
  S_fn <- surv.exp_kofn(x)
  function(t, ...) f_fn(t) / S_fn(t)
}


#' @rdname exp_kofn
#' @importFrom stats density dexp
#' @export
density.exp_kofn <- function(x, ...) {
  rates <- x$rates
  k <- x$k
  function(t, log = FALSE, ...) {
    vals <- vapply(t, function(ti) {
      surv <- exp(-rates * ti)
      dens <- rates * surv
      kofn_density_value(dens, surv, k)
    }, numeric(1L))
    if (isTRUE(log)) log(vals) else vals
  }
}
