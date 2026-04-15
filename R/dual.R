# ==========================================================================
# Dual structures
# ==========================================================================


#' Default dual: lazy wrapper
#'
#' Returns a `dual_of_system` object carrying the original structure.
#' `phi_dual(state) = 1 - phi(original, 1 - state)` is evaluated on demand.
#' All other generics fall through to `dist_structure` defaults.
#'
#' @rdname dual
#' @export
dual.dist_structure <- function(x) {
  structure(
    list(original = x, m = ncomponents(x)),
    class = c("dual_of_system", "dist_structure",
              "univariate_dist", "continuous_dist", "dist")
  )
}


#' @export
ncomponents.dual_of_system <- function(x) x$m


#' @export
phi.dual_of_system <- function(x, state) {
  state <- as.integer(state)
  1L - phi(x$original, 1L - state)
}


#' Dual of dual is the original (involution)
#'
#' @rdname dual
#' @export
dual.dual_of_system <- function(x) x$original
