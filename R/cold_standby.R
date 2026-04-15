# ==========================================================================
# Cold standby system
# ==========================================================================
#
# A cold-standby system has m components arranged sequentially: one is
# active at a time, others wait with zero hazard. When the active one
# fails, a spare takes over (assumed instantaneous, perfect switching).
# System lifetime = T_1 + T_2 + ... + T_m, the sum of independent
# component lifetimes.
#
# Cold standby is NOT a coherent system in the structure-function sense
# (the topology is temporal succession, not order statistics). Therefore
# `cold_standby_dist` does NOT inherit `dist_structure`. It carries
# components and provides standard dist methods (surv, cdf, sampler,
# mean, expectation), but has no phi or min_paths.
#
# Special cases with closed-form aggregate distributions:
#   - iid Exp(rate)        ->  Gamma(shape = m, rate = rate)
#   - heterogeneous Exp    ->  hypoexponential (no simple closed form)
#   - independent normals  ->  Normal(sum mu, sqrt(sum var))
#
# The general default uses Monte Carlo for surv/cdf and exact (sum of
# component means) for `mean`.
# ==========================================================================


#' Cold-standby system distribution
#'
#' Constructs a distribution representing a cold-standby system: one
#' component active at a time with perfect, instantaneous switching to
#' the next spare upon failure. System lifetime equals the sum of
#' independent component lifetimes.
#'
#' Cold standby is not a coherent system in the structure-function sense
#' (its topology is temporal succession, not an order statistic), so the
#' returned object does not inherit `dist_structure`. It IS a `dist`
#' (with `surv`, `cdf`, `sampler`, `mean` available) and exposes
#' [ncomponents()] and [component()] for inspection.
#'
#' Defaults: `sampler` is exact (sample each component independently and
#' sum); `mean` is exact (sum of component means); `surv` and `cdf` use
#' Monte Carlo with a default of `1e5` simulated lifetimes (override via
#' the `mc` argument when calling `surv(x)` or by overriding the method
#' on a specialized subclass with a closed-form aggregate distribution
#' such as Gamma for iid exponential components).
#'
#' @param components List of `dist` objects representing per-stage
#'   component lifetimes.
#' @return An object of class
#'   `c("cold_standby_dist", "univariate_dist", "continuous_dist", "dist")`.
#' @examples
#' # Cold standby of 3 iid Exp(1) components: aggregate is Gamma(3, 1)
#' sys <- cold_standby_dist(replicate(3,
#'   algebraic.dist::exponential(1), simplify = FALSE))
#' mean(sys)  # = 3
#' @export
cold_standby_dist <- function(components) {
  stopifnot(is.list(components), length(components) >= 1L)
  structure(
    list(components = components, m = length(components)),
    class = c("cold_standby_dist", "univariate_dist",
              "continuous_dist", "dist")
  )
}


#' @export
ncomponents.cold_standby_dist <- function(x) x$m


#' @export
component.cold_standby_dist <- function(x, j, ...) {
  stopifnot(j >= 1L, j <= x$m)
  x$components[[j]]
}


#' @rdname cold_standby_dist
#' @param x A `cold_standby_dist` object.
#' @param ... Ignored.
#' @export
sampler.cold_standby_dist <- function(x, ...) {
  m <- x$m
  comp_samplers <- lapply(x$components, algebraic.dist::sampler)
  function(n, ...) {
    rowSums(sample_component_matrix(comp_samplers, n))
  }
}


#' Mean of a cold-standby system: sum of component means
#'
#' @param x A `cold_standby_dist` object.
#' @param ... Ignored.
#' @return Numeric scalar.
#' @export
mean.cold_standby_dist <- function(x, ...) {
  sum(vapply(x$components, mean, numeric(1L)))
}


#' Cold-standby survival via Monte Carlo
#'
#' Estimates `S(t) = P(T_1 + ... + T_m > t)` by simulating `mc` system
#' lifetimes and computing the empirical fraction exceeding `t`.
#' Override on specialized subclasses (e.g., iid exponential -> Gamma)
#' for exact closed forms.
#'
#' @rdname cold_standby_dist
#' @export
surv.cold_standby_dist <- function(x, ...) {
  samp <- sampler.cold_standby_dist(x)
  function(t, mc = 1e5L, ...) {
    samples <- samp(mc)
    vapply(t, function(ti) mean(samples > ti), numeric(1L))
  }
}


#' Cold-standby CDF via Monte Carlo: 1 - surv
#'
#' @rdname cold_standby_dist
#' @export
cdf.cold_standby_dist <- function(x, ...) {
  S <- surv.cold_standby_dist(x)
  function(t, mc = 1e5L, ...) 1 - S(t, mc = mc)
}
