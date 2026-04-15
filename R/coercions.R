# ==========================================================================
# Coercions
# ==========================================================================
#
# as_dist_structure(x): promote a plain dist (e.g., from algebraic.dist)
# into a 1-component dist_structure (a trivial series of one component).
# Useful for polymorphic functions that accept dist_structure inputs but
# are also handed plain dists.
# ==========================================================================


#' @rdname as_dist_structure
#' @export
as_dist_structure.dist_structure <- function(x, ...) x


#' @rdname as_dist_structure
#' @export
as_dist_structure.dist <- function(x, ...) {
  series_dist(list(x))
}
