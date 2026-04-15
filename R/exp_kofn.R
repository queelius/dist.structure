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
#' are independent exponentials. Closed-form methods for `surv`, `cdf`,
#' and `sampler`; `density`, `hazard`, and `mean` fall back to numerical
#' defaults (or can be overridden by users in specific applications).
#'
#' @param k Minimum number of functioning components for system operation.
#' @param rates Positive numeric vector of length `m` with `m >= k`.
#' @return An object of class
#'   `c("exp_kofn", "kofn_dist", "coherent_dist", "dist_structure",
#'   "univariate_dist", "continuous_dist", "dist")`.
#' @examples
#' sys <- exp_kofn(k = 2, rates = c(1, 2, 3))
#' algebraic.dist::surv(sys)(1)
#' @export
exp_kofn <- function(k, rates) {
  stopifnot(is.numeric(rates), length(rates) >= 1L, all(rates > 0))
  m <- length(rates)
  stopifnot(k >= 1L, k <= m)
  components <- lapply(rates, function(r) algebraic.dist::exponential(r))
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
cdf.exp_kofn <- function(x, ...) {
  S <- surv.exp_kofn(x)
  function(t, ...) 1 - S(t)
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
