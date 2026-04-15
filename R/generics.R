# ==========================================================================
# S3 generic functions
# ==========================================================================


#' Number of components
#'
#' @param x A [dist_structure] object.
#' @return A positive integer `m`, the number of components.
#' @export
ncomponents <- function(x) UseMethod("ncomponents")


#' The j-th component as a dist
#'
#' Returns the j-th component as an `algebraic.dist::dist` object with
#' parameters baked in (so `sampler`, `surv`, etc. work directly on it).
#'
#' @param x A [dist_structure] object.
#' @param j Integer component index in `1:ncomponents(x)`.
#' @param ... Implementation-specific arguments (e.g., parameter overrides
#'   for lazily parameterized systems).
#' @return An object inheriting from `dist`.
#' @export
component <- function(x, j, ...) UseMethod("component")


#' Structure function
#'
#' Evaluate the coherent structure function `phi: {0, 1}^m -> {0, 1}` at a
#' component state vector `state`. By convention `state[j] = 1` means
#' component `j` is functioning; `phi(state) = 1` means the system is
#' functioning.
#'
#' @param x A [dist_structure] object.
#' @param state Integer or logical vector of length `ncomponents(x)` in
#'   `{0, 1}`.
#' @return Integer scalar, `0` or `1`.
#' @export
phi <- function(x, state) UseMethod("phi")


#' Minimal path sets
#'
#' @param x A [dist_structure] object.
#' @return A list of integer vectors.
#' @export
min_paths <- function(x) UseMethod("min_paths")


#' Minimal cut sets
#'
#' @param x A [dist_structure] object.
#' @return A list of integer vectors.
#' @export
min_cuts <- function(x) UseMethod("min_cuts")


#' System signature
#'
#' Samaniego's signature: `s = (s_1, ..., s_m)` where `s_k` is the
#' probability the system fails at the k-th component failure under
#' i.i.d. absolutely-continuous component lifetimes. Only depends on the
#' structure (phi), not on the component distribution. A default method
#' on `dist_structure` enumerates the `m!` orderings; this is feasible
#' for `m` up to about 8 or 9. Specialized subclasses should override
#' with closed-form expressions.
#'
#' @param x A [dist_structure] object.
#' @return Numeric vector of length `ncomponents(x)` summing to 1.
#' @export
system_signature <- function(x) UseMethod("system_signature")


#' Critical states of a component
#'
#' @param x A [dist_structure] object.
#' @param j Component index.
#' @return Integer matrix with `m - 1` columns; rows are the states of the
#'   other components for which `j` is critical.
#' @export
critical_states <- function(x, j) UseMethod("critical_states")


#' System lifetime from component times
#'
#' @param x A [dist_structure] object.
#' @param times Non-negative numeric vector of length `ncomponents(x)`.
#' @return Scalar system lifetime.
#' @export
system_lifetime <- function(x, times) UseMethod("system_lifetime")


#' Per-component censoring from system observation
#'
#' @param x A [dist_structure] object.
#' @param times Non-negative numeric vector of length `ncomponents(x)`.
#' @return A list with `system_time` (scalar) and `component_status`
#'   (character vector of length `m`, values in `"exact"`, `"left"`,
#'   `"right"`).
#' @export
system_censoring <- function(x, times) UseMethod("system_censoring")


#' Dual structure
#'
#' The dual structure satisfies `phi_dual(state) = 1 - phi(1 - state)`. The
#' dual of a series system is parallel; the dual of k-out-of-n is
#' (n - k + 1)-out-of-n.
#'
#' @param x A [dist_structure] object.
#' @return A [dist_structure] object representing the dual.
#' @export
dual <- function(x) UseMethod("dual")


#' Coherence axiom check
#'
#' @param x A [dist_structure] object.
#' @return `TRUE` if monotone and every component relevant.
#' @export
is_coherent <- function(x) UseMethod("is_coherent")


#' Birnbaum structural importance
#'
#' Fraction of the `2^(m - 1)` states of the other components for which
#' component `j` is critical.
#'
#' @param x A [dist_structure] object.
#' @param j Component index.
#' @return Numeric scalar in `[0, 1]`.
#' @export
structural_importance <- function(x, j) UseMethod("structural_importance")


#' System reliability polynomial
#'
#' `R(p) = E[phi(X)]` where `X_j ~ Bernoulli(p_j)` independent. The
#' multilinear extension of `phi` to component reliabilities.
#'
#' @param x A [dist_structure] object.
#' @param p Numeric vector of length `ncomponents(x)` or a scalar recycled
#'   to all components; values in `[0, 1]`.
#' @return Numeric scalar in `[0, 1]`.
#' @export
reliability <- function(x, p) UseMethod("reliability")
