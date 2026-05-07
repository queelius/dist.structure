# ==========================================================================
# Virtual base class for structured distributions
# ==========================================================================


#' @title Virtual base class for structured distributions
#'
#' @description
#' `dist_structure` is the virtual S3 class for distributions whose random
#' variable has internal component structure: coherent reliability
#' systems decomposed into components arranged by a structure function.
#' Every concrete implementation (`coherent_dist`, `series_dist`,
#' `parallel_dist`, `kofn_dist`, `bridge_dist`, or user-defined
#' subclasses) should include `"dist_structure"`, the `algebraic.dist`
#' ancestor `"univariate_dist"`, and `"dist"` in its class vector.
#'
#' @details
#' Concrete implementations provide S3 methods for the generics in this
#' package. The minimum required methods are [ncomponents()], [component()],
#' and one of [phi()] or [min_paths()]; every other generic has a default
#' method on `dist_structure` that composes the primitives.
#'
#' If both `phi.<class>()` and `min_paths.<class>()` are provided, the
#' implementor is responsible for keeping them consistent: `phi` derives
#' the in-package generics `reliability`, `critical_states`, and
#' `is_coherent`; `min_paths` derives `min_cuts` and `system_signature`.
#' Inconsistent implementations produce silently inconsistent results.
#'
#' @return
#' This help topic documents the virtual base class together with the
#' three distribution-algebra default methods that compose component
#' distributions through the topology:
#'
#' - [surv.dist_structure()] returns a closure `function(t, ...)` that
#'   evaluates the system survival function via the reliability identity
#'   `S_sys(t) = R(S_1(t), ..., S_m(t))`.
#' - [cdf.dist_structure()] returns a closure `function(t, ...)` equal
#'   to `1 - surv(x)(t)`.
#' - [sampler.dist_structure()] returns a closure `function(n, ...)`
#'   that draws `n` independent system lifetimes by sampling each
#'   component and combining via [system_lifetime()].
#'
#' Concrete subclasses override any of these for closed-form speed; see
#' the closed-form specializations under the See Also entries.
#'
#' @seealso [validate_dist_structure()] for an implementor-side
#'   construction-time validator. [phi()] and [min_paths()] for the
#'   bidirectional protocol primitives. The closed-form families
#'   ([exp_series()], [wei_kofn()], etc.) for reference implementations.
#'
#' @name dist_structure
NULL


#' Predicate for dist_structure objects
#'
#' @param x Any object.
#' @return `TRUE` if `x` inherits from `"dist_structure"`.
#' @export
is_dist_structure <- function(x) inherits(x, "dist_structure")


#' Validate that an object satisfies the dist_structure protocol
#'
#' Checks that `x` declares `"dist_structure"` in its class chain and
#' provides the three required generics: [ncomponents()], [component()],
#' and at least one of [phi()] or [min_paths()].
#'
#' Useful in subclass constructors to fail fast with a clear error
#' message rather than at first method dispatch (where the user sees
#' opaque "no applicable method" errors). A typical pattern:
#'
#' ```r
#' my_subclass <- function(...) {
#'   obj <- structure(list(...), class = c("my_subclass", "dist_structure",
#'                                          "univariate_dist", "continuous_dist", "dist"))
#'   validate_dist_structure(obj)
#'   obj
#' }
#' ```
#'
#' @param x An object claiming to satisfy the `dist_structure` protocol.
#' @return `TRUE` (invisibly) if all checks pass. Stops with an
#'   informative error otherwise.
#' @export
#' @examples
#' # Success: a valid dist_structure (any built-in topology shortcut works).
#' validate_dist_structure(series_dist(replicate(3,
#'   algebraic.dist::exponential(1), simplify = FALSE)))
#'
#' # Failure: an object that declares dist_structure but provides no
#' # primitives. Wrapped in tryCatch for the example's success status.
#' bad <- structure(list(),
#'                  class = c("not_a_real_class", "dist_structure",
#'                            "univariate_dist", "continuous_dist", "dist"))
#' tryCatch(validate_dist_structure(bad),
#'          error = function(e) conditionMessage(e))
validate_dist_structure <- function(x) {
  if (!inherits(x, "dist_structure")) {
    stop("`x` must declare 'dist_structure' in its class chain; ",
         "got class = c(", paste(shQuote(class(x)), collapse = ", "), ").",
         call. = FALSE)
  }
  m <- tryCatch(ncomponents(x),
                error = function(e) {
                  stop("`x` does not implement `ncomponents()`: ", e$message,
                       call. = FALSE)
                })
  if (!is.numeric(m) || length(m) != 1L || m < 1L) {
    stop("`ncomponents(x)` must return a single positive integer; got ",
         paste(format(m), collapse = ", "), ".", call. = FALSE)
  }
  tryCatch(component(x, 1L),
           error = function(e) {
             stop("`x` does not implement `component()`: ", e$message,
                  call. = FALSE)
           })
  # Check via S3 method-table lookup (rather than calling phi/min_paths
  # and catching errors) because the default phi.dist_structure and
  # min_paths.dist_structure are mutually recursive: if a subclass
  # provides neither, calling the dispatched generic would loop. Looking
  # for a non-default registered method on any subclass-specific class is
  # the correct check.
  #
  # Note: `coherent_dist` is intentionally NOT stripped from user_classes.
  # `min_paths.coherent_dist` is registered (returns x$min_paths from the
  # constructor), so any subclass of `coherent_dist` legitimately
  # inherits a working min_paths implementation; validate should accept
  # it. Stripping `coherent_dist` would incorrectly reject those
  # subclasses.
  user_classes <- setdiff(
    class(x),
    c("dist_structure", "univariate_dist", "continuous_dist", "dist")
  )
  has_method <- function(generic) {
    for (cls in user_classes) {
      m <- utils::getS3method(generic, cls, optional = TRUE)
      if (!is.null(m)) return(TRUE)
    }
    FALSE
  }
  if (!has_method("phi") && !has_method("min_paths")) {
    stop("`x` must implement at least one of `phi.<class>()` or ",
         "`min_paths.<class>()`; neither was found on any subclass ",
         "(only the mutually-recursive `dist_structure` defaults).",
         call. = FALSE)
  }
  invisible(TRUE)
}


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
