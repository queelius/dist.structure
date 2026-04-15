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
#' sum); `mean` is exact when every component implements `mean()`
#' (otherwise falls back to Monte Carlo via the sampler); `surv` and
#' `cdf` use Monte Carlo with a default of `1e5` simulated lifetimes
#' (override via the `mc` argument). The returned `surv` / `cdf`
#' closures cache their samples after the first call so subsequent
#' evaluations at different `t` values are deterministic given the same
#' `mc`.
#'
#' Methods inherited from `algebraic.dist::univariate_dist` that require
#' `density` or `sup` (notably `expectation`, `vcov`) are NOT supported
#' on cold-standby objects out of the box; specialized subclasses with
#' closed-form aggregate distributions (e.g., iid exponential collapses
#' to `Gamma(m, rate)`) should provide their own methods.
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
#' Computes `sum_j E[T_j]` exactly when every component implements a
#' `mean()` method. Falls back to a Monte Carlo estimate from the
#' sampler when any component lacks an exact mean (e.g., a
#' `dist_structure` component whose `mean` would route through
#' `algebraic.dist::univariate_dist` and require `density` / `sup`,
#' which are not provided on `dist_structure` by default).
#'
#' @param x A `cold_standby_dist` object.
#' @param ... Passed to the Monte Carlo fallback as `mc` (default 1e5).
#' @return Numeric scalar.
#' @export
mean.cold_standby_dist <- function(x, ...) {
  exact <- tryCatch(
    sum(vapply(x$components, mean, numeric(1L))),
    error = function(e) NULL
  )
  if (!is.null(exact)) return(exact)
  args <- list(...)
  mc <- if ("mc" %in% names(args)) args$mc else 1e5L
  mean(sampler.cold_standby_dist(x)(mc))
}


#' Cold-standby survival via Monte Carlo (with cached samples)
#'
#' Estimates `S(t) = P(T_1 + ... + T_m > t)` by simulating `mc` system
#' lifetimes and computing the empirical fraction exceeding `t`. The
#' returned closure caches its sample vector: the first call generates
#' `mc` samples, subsequent calls reuse them as long as `mc` is
#' unchanged (a different `mc` triggers a fresh draw). This makes
#' repeated `S(t)` calls deterministic given the same `mc`.
#'
#' For reproducibility across calls to `surv()` itself (i.e., between
#' separately constructed closures), set the RNG seed externally via
#' `set.seed()` before invoking `surv(x)`. Override the method on a
#' subclass with an exact aggregate distribution (e.g., iid exponential
#' collapses to `Gamma(m, rate)`) when an analytic form is available.
#'
#' @rdname cold_standby_dist
#' @export
surv.cold_standby_dist <- function(x, ...) {
  samp <- sampler.cold_standby_dist(x)
  cache <- new.env(parent = emptyenv())
  cache$samples <- NULL
  cache$mc <- NA_integer_
  function(t, mc = 1e5L, ...) {
    mc <- as.integer(mc)
    if (is.null(cache$samples) || cache$mc != mc) {
      cache$samples <- samp(mc)
      cache$mc <- mc
    }
    samples <- cache$samples
    vapply(t, function(ti) mean(samples > ti), numeric(1L))
  }
}


#' Cold-standby CDF: `1 - surv(x)(t, mc = mc)`
#'
#' Computes the CDF by deferring to a fresh `surv` closure. The CDF
#' closure has its own sample cache (independent of any external `surv`
#' closure), so `cdf(x)(t) + surv(x)(t)` is generally not exactly 1
#' unless callers reuse a single `surv` closure: prefer
#' `S <- surv(x); F_t <- 1 - S(t)` over computing both independently.
#'
#' @rdname cold_standby_dist
#' @export
cdf.cold_standby_dist <- function(x, ...) {
  S <- surv.cold_standby_dist(x)
  function(t, mc = 1e5L, ...) 1 - S(t, mc = mc)
}
