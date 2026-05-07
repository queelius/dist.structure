# ==========================================================================
# Closed-form series of exponential components
# ==========================================================================
#
# A series of m independent exponentials with rates lambda_j has itself
# an exponential system lifetime with rate = sum(lambda). This gives
# closed-form survival, sampler, hazard, density, and mean.
# ==========================================================================


#' Series of exponential components (closed form)
#'
#' Constructs a `dist_structure` representing a series system whose
#' components are independent exponentials. The system lifetime is itself
#' an `Exp(sum(rates))` distribution; all dist-level queries have
#' closed-form expressions that bypass the general default methods.
#'
#' @param rates Positive numeric vector of length `m`: per-component
#'   exponential rates.
#' @return
#' `exp_series()` returns an object of class
#'   `c("exp_series", "series_dist", "coherent_dist", "dist_structure",
#'   "univariate_dist", "continuous_dist", "dist")`.
#'
#' The associated S3 methods return:
#' - `surv()`, `cdf()`, `density()`, `hazard()`: a closure `function(t, ...)`
#'   evaluating the named quantity at `t`.
#' - `sampler()`: a closure `function(n, ...)` returning `n` random
#'   variates from the system lifetime distribution.
#' - `mean()`: a numeric scalar (the mean system lifetime).
#' @examples
#' sys <- exp_series(c(0.5, 0.3, 0.2))
#' algebraic.dist::surv(sys)(1)  # equals exp(-sum(rates) * 1)
#' mean(sys)                     # equals 1 / sum(rates)
#' @export
exp_series <- function(rates) {
  stopifnot(is.numeric(rates), length(rates) >= 1L, all(rates > 0))
  components <- lapply(rates, algebraic.dist::exponential)
  obj <- series_dist(components)
  obj$rates <- as.numeric(rates)
  obj$total_rate <- sum(rates)
  class(obj) <- c("exp_series", class(obj))
  obj
}


#' @rdname exp_series
#' @param x An `exp_series` object.
#' @param ... Ignored.
#' @export
surv.exp_series <- function(x, ...) {
  lam <- x$total_rate
  function(t, ...) exp(-lam * t)
}


#' @rdname exp_series
#' @export
sampler.exp_series <- function(x, ...) {
  lam <- x$total_rate
  function(n, ...) stats::rexp(n, rate = lam)
}


#' @rdname exp_series
#' @export
mean.exp_series <- function(x, ...) 1 / x$total_rate


#' @rdname exp_series
#' @importFrom stats density dexp
#' @export
density.exp_series <- function(x, ...) {
  lam <- x$total_rate
  function(t, log = FALSE, ...) stats::dexp(t, rate = lam, log = log)
}


#' @rdname exp_series
#' @method hazard exp_series
#' @importFrom algebraic.dist hazard
#' @export
hazard.exp_series <- function(x, ...) {
  lam <- x$total_rate
  function(t, log.p = FALSE, ...) {
    rep(if (isTRUE(log.p)) log(lam) else lam, length(t))
  }
}
