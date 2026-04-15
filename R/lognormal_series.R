# ==========================================================================
# Closed-form series of Lognormal components
# ==========================================================================
#
# For a series of m independent Lognormals with parameters mu_j and sd_j:
#   S_sys(t) = prod_j stats::plnorm(t, meanlog = mu_j, sdlog = sd_j,
#                                    lower.tail = FALSE)
# Sampler generates m independent Lognormals and takes the min.
# ==========================================================================


#' Series of independent Lognormal components (closed form)
#'
#' Constructs a `dist_structure` for a series system whose components are
#' independent Lognormals. Closed-form `surv` is the product of
#' per-component upper-tail probabilities; `cdf` is `1 - surv`; `sampler`
#' generates m independent Lognormals and takes the min.
#'
#' @param meanlogs Numeric vector of length `m`: per-component meanlog
#'   parameters.
#' @param sdlogs Positive numeric vector of length `m`: per-component
#'   sdlog parameters.
#' @return An object of class
#'   `c("lognormal_series", "series_dist", "coherent_dist",
#'   "dist_structure", "univariate_dist", "continuous_dist", "dist")`.
#' @examples
#' sys <- lognormal_series(meanlogs = c(0, 1), sdlogs = c(1, 0.5))
#' algebraic.dist::surv(sys)(1)
#' @export
lognormal_series <- function(meanlogs, sdlogs) {
  stopifnot(length(meanlogs) == length(sdlogs), all(sdlogs > 0))
  m <- length(meanlogs)
  components <- lapply(seq_len(m), function(j) {
    algebraic.dist::lognormal(meanlog = meanlogs[j], sdlog = sdlogs[j])
  })
  obj <- series_dist(components)
  obj$meanlogs <- as.numeric(meanlogs)
  obj$sdlogs <- as.numeric(sdlogs)
  class(obj) <- c("lognormal_series", class(obj))
  obj
}


#' @rdname lognormal_series
#' @param x A `lognormal_series` object.
#' @param ... Ignored.
#' @export
surv.lognormal_series <- function(x, ...) {
  series_surv_product(stats::plnorm,
                      list(meanlog = x$meanlogs, sdlog = x$sdlogs))
}


#' @rdname lognormal_series
#' @export
cdf.lognormal_series <- function(x, ...) {
  S <- surv.lognormal_series(x)
  function(t, ...) 1 - S(t)
}


#' @rdname lognormal_series
#' @export
sampler.lognormal_series <- function(x, ...) {
  samplers <- make_component_samplers(stats::rlnorm,
                                      meanlog = x$meanlogs,
                                      sdlog = x$sdlogs)
  function(n, ...) {
    apply(sample_component_matrix(samplers, n), 1L, min)
  }
}
