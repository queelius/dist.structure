# ==========================================================================
# Virtual base class for structured distributions
# ==========================================================================


#' Virtual base class for structured distributions
#'
#' `dist_structure` is a virtual S3 class: every concrete implementation
#' (`coherent_dist`, `series_dist`, `parallel_dist`, `kofn_dist`,
#' `bridge_dist`, or user-defined subclasses) should include both
#' `"dist_structure"` and the algebraic.dist ancestor `"univariate_dist"`
#' (for system lifetime, which is scalar) plus `"dist"` in its class vector.
#'
#' Concrete implementations provide S3 methods for the generics in this
#' package. The minimum required methods are [ncomponents()], [component()],
#' and one of [phi()] or [min_paths()]; every other generic has a default
#' method on `dist_structure` that composes the primitives.
#'
#' @name dist_structure
NULL


#' Predicate for dist_structure objects
#'
#' @param x Any object.
#' @return `TRUE` if `x` inherits from `"dist_structure"`.
#' @export
is_dist_structure <- function(x) inherits(x, "dist_structure")


#' Format a dist_structure object
#'
#' @param x A [dist_structure] object.
#' @param ... Ignored.
#' @return Character vector suitable for [cat()].
#' @export
format.dist_structure <- function(x, ...) {
  cls <- setdiff(class(x),
                 c("dist_structure", "univariate_dist",
                   "continuous_dist", "dist"))
  head <- if (length(cls)) {
    paste0("<dist_structure: ", cls[[1L]], ">")
  } else {
    "<dist_structure>"
  }
  m <- tryCatch(ncomponents(x), error = function(e) NA_integer_)
  body <- if (is.na(m)) {
    "  components: (ncomponents method not defined)"
  } else {
    paste0("  components: ", m)
  }
  c(head, body)
}


#' Print a dist_structure object
#'
#' @param x A [dist_structure] object.
#' @param ... Passed to [format()].
#' @return `x`, invisibly.
#' @export
print.dist_structure <- function(x, ...) {
  cat(format(x, ...), sep = "\n")
  invisible(x)
}
