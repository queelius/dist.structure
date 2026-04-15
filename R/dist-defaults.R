# ==========================================================================
# Distribution defaults on `dist_structure`
# ==========================================================================
#
# System-level distribution methods composed from component distributions
# and topology. Implementors override for closed-form speed.
# ==========================================================================


#' System survival via reliability composition
#'
#' Returns a closure `function(t, ...)` where `surv(x)(t)` equals
#' `reliability(x, surv(component(x, j))(t) for each j)`. This is the
#' classical identity `S_sys(t) = R(S(t))` realized as a composition.
#'
#' @rdname dist_structure
#' @param x A [dist_structure] object.
#' @param ... Ignored.
#' @export
surv.dist_structure <- function(x, ...) {
  m <- ncomponents(x)
  comp_surv_fns <- lapply(seq_len(m), function(j) {
    algebraic.dist::surv(component(x, j))
  })
  function(t, ...) {
    vapply(t, function(ti) {
      ps <- vapply(comp_surv_fns, function(S_j) S_j(ti), numeric(1L))
      reliability(x, ps)
    }, numeric(1L))
  }
}


#' System CDF: `1 - surv(x)(t)`
#'
#' @rdname dist_structure
#' @export
cdf.dist_structure <- function(x, ...) {
  S <- algebraic.dist::surv(x)
  function(t, ...) 1 - S(t, ...)
}


#' System sampler via component samplers and system_lifetime
#'
#' Returns a closure `function(n, ...)` that draws `n` system lifetimes by
#' sampling each component independently and applying [system_lifetime()]
#' to combine.
#'
#' @rdname dist_structure
#' @export
sampler.dist_structure <- function(x, ...) {
  m <- ncomponents(x)
  comp_samplers <- lapply(seq_len(m), function(j) {
    algebraic.dist::sampler(component(x, j))
  })
  function(n, ...) {
    comp_matrix <- sample_component_matrix(comp_samplers, n)
    apply(comp_matrix, 1L, function(times) system_lifetime(x, times))
  }
}
