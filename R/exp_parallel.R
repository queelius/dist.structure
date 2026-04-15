# ==========================================================================
# Closed-form parallel of exponential components
# ==========================================================================
#
# A parallel system of m independent exponentials with rates lambda_j has
#   F_sys(t) = prod_j (1 - exp(-lambda_j t))
#   S_sys(t) = 1 - F_sys(t)
# For iid components, E[max] = (1/lambda)(1 + 1/2 + ... + 1/m).
# For the heterogeneous case, E[max] = sum over nonempty subsets S of
# (-1)^(|S|+1) / sum(rates[S]) via inclusion-exclusion.
# ==========================================================================


#' Parallel of exponential components (closed form)
#'
#' Constructs a `dist_structure` representing a parallel system whose
#' components are independent exponentials. Closed-form methods are
#' provided for `surv`, `cdf`, `sampler`, and `mean` (the last via
#' inclusion-exclusion over the `2^m - 1` non-empty component subsets).
#'
#' @param rates Positive numeric vector of length `m`.
#' @return An object of class
#'   `c("exp_parallel", "parallel_dist", "coherent_dist", "dist_structure",
#'   "univariate_dist", "continuous_dist", "dist")`.
#' @examples
#' sys <- exp_parallel(c(1, 2, 3))
#' algebraic.dist::surv(sys)(1)
#' @export
exp_parallel <- function(rates) {
  stopifnot(is.numeric(rates), length(rates) >= 1L, all(rates > 0))
  components <- lapply(rates, function(r) algebraic.dist::exponential(r))
  obj <- parallel_dist(components)
  obj$rates <- as.numeric(rates)
  class(obj) <- c("exp_parallel", class(obj))
  obj
}


#' @rdname exp_parallel
#' @param x An `exp_parallel` object.
#' @param ... Ignored.
#' @export
cdf.exp_parallel <- function(x, ...) {
  rates <- x$rates
  function(t, ...) {
    vapply(t, function(ti) prod(1 - exp(-rates * ti)), numeric(1L))
  }
}


#' @rdname exp_parallel
#' @export
surv.exp_parallel <- function(x, ...) {
  FF <- cdf.exp_parallel(x)
  function(t, ...) 1 - FF(t)
}


#' @rdname exp_parallel
#' @export
sampler.exp_parallel <- function(x, ...) {
  rates <- x$rates
  samplers <- lapply(rates, function(r) function(n) stats::rexp(n, rate = r))
  function(n, ...) {
    apply(sample_component_matrix(samplers, n), 1L, max)
  }
}


#' @rdname exp_parallel
#' @export
mean.exp_parallel <- function(x, ...) {
  rates <- x$rates
  m <- length(rates)
  # Fast path: iid components. E[max of m iid Exp(lambda)] = H_m / lambda.
  if (length(unique(rates)) == 1L) {
    return(sum(1 / seq_len(m)) / rates[[1L]])
  }
  # General heterogeneous case: inclusion-exclusion over 2^m - 1 subsets.
  total <- 0
  for (sz in seq_len(m)) {
    subsets <- utils::combn(m, sz, simplify = FALSE)
    for (S in subsets) {
      total <- total + (-1)^(sz + 1L) / sum(rates[S])
    }
  }
  total
}
