# ==========================================================================
# Closed-form series of Gamma components
# ==========================================================================
#
# For a series of m independent Gammas with shapes k_j and rates r_j:
#   S_sys(t) = prod_j P(Gamma_j > t) = prod_j stats::pgamma(t, k_j, r_j,
#                                                            lower.tail = FALSE)
# No closed-form CDF inverse for the system in general; sampler uses min
# of independently sampled Gammas.
# ==========================================================================


#' Series of independent Gamma components (closed form)
#'
#' Constructs a `dist_structure` for a series system whose components are
#' independent Gammas. Closed-form `surv` is evaluated by the product of
#' per-component upper-tail probabilities; `cdf` is `1 - surv`; `sampler`
#' generates m independent Gammas and takes the min.
#'
#' @param shapes Positive numeric vector of length `m`: per-component
#'   Gamma shape parameters.
#' @param rates Positive numeric vector of length `m`: per-component
#'   Gamma rate parameters.
#' @return An object of class
#'   `c("gamma_series", "series_dist", "coherent_dist", "dist_structure",
#'   "univariate_dist", "continuous_dist", "dist")`.
#' @examples
#' sys <- gamma_series(shapes = c(2, 3), rates = c(1, 2))
#' algebraic.dist::surv(sys)(1)
#' @export
gamma_series <- function(shapes, rates) {
  stopifnot(length(shapes) == length(rates),
            all(shapes > 0), all(rates > 0))
  m <- length(shapes)
  components <- lapply(seq_len(m), function(j) {
    algebraic.dist::gamma_dist(shape = shapes[j], rate = rates[j])
  })
  obj <- series_dist(components)
  obj$shapes <- as.numeric(shapes)
  obj$rates <- as.numeric(rates)
  class(obj) <- c("gamma_series", class(obj))
  obj
}


#' @rdname gamma_series
#' @param x A `gamma_series` object.
#' @param ... Ignored.
#' @export
surv.gamma_series <- function(x, ...) {
  shapes <- x$shapes
  rates <- x$rates
  function(t, ...) {
    vapply(t, function(ti) {
      prod(stats::pgamma(ti, shape = shapes, rate = rates,
                         lower.tail = FALSE))
    }, numeric(1L))
  }
}


#' @rdname gamma_series
#' @export
cdf.gamma_series <- function(x, ...) {
  S <- surv.gamma_series(x)
  function(t, ...) 1 - S(t)
}


#' @rdname gamma_series
#' @export
sampler.gamma_series <- function(x, ...) {
  shapes <- x$shapes
  rates <- x$rates
  m <- length(shapes)
  samplers <- lapply(seq_len(m), function(j) {
    function(n) stats::rgamma(n, shape = shapes[j], rate = rates[j])
  })
  function(n, ...) {
    apply(sample_component_matrix(samplers, n), 1L, min)
  }
}
