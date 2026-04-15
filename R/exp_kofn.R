# ==========================================================================
# Closed-form k-out-of-n of heterogeneous exponential components
# ==========================================================================
#
# For k-of-m of heterogeneous exponentials:
#   S_sys(t) = P(at least k of m exponentials are alive at t)
#            = sum_{A subset of [m], |A|>=k} prod_{j in A} exp(-rate_j t)
#              * prod_{j notin A} (1 - exp(-rate_j t))
# Sampler uses the (m - k + 1)-th order statistic of independently sampled
# component lifetimes.
# ==========================================================================


#' k-out-of-n system of independent exponential components (closed form)
#'
#' Constructs a `dist_structure` for a k-out-of-m system whose components
#' are independent exponentials. Closed-form methods for `surv`, `cdf`,
#' and `sampler`; `density`, `hazard`, and `mean` fall back to numerical
#' defaults (or can be overridden by users in specific applications).
#'
#' @param k Minimum number of functioning components for system operation.
#' @param rates Positive numeric vector of length `m` with `m >= k`.
#' @return An object of class
#'   `c("exp_kofn", "kofn_dist", "coherent_dist", "dist_structure",
#'   "univariate_dist", "continuous_dist", "dist")`.
#' @examples
#' sys <- exp_kofn(k = 2, rates = c(1, 2, 3))
#' algebraic.dist::surv(sys)(1)
#' @export
exp_kofn <- function(k, rates) {
  stopifnot(is.numeric(rates), length(rates) >= 1L, all(rates > 0))
  m <- length(rates)
  stopifnot(k >= 1L, k <= m)
  components <- lapply(rates, function(r) algebraic.dist::exponential(r))
  obj <- kofn_dist(k, components)
  obj$rates <- as.numeric(rates)
  class(obj) <- c("exp_kofn", class(obj))
  obj
}


#' @rdname exp_kofn
#' @param x An `exp_kofn` object.
#' @param ... Ignored.
#' @export
surv.exp_kofn <- function(x, ...) {
  rates <- x$rates
  k <- x$k
  m <- length(rates)
  function(t, ...) {
    vapply(t, function(ti) {
      comp_surv <- exp(-rates * ti)
      comp_fail <- 1 - comp_surv
      total <- 0
      for (sz in k:m) {
        subsets <- utils::combn(m, sz, simplify = FALSE)
        for (A in subsets) {
          alive <- prod(comp_surv[A])
          failed <- setdiff(seq_len(m), A)
          dead <- if (length(failed) == 0L) 1 else prod(comp_fail[failed])
          total <- total + alive * dead
        }
      }
      total
    }, numeric(1L))
  }
}


#' @rdname exp_kofn
#' @export
cdf.exp_kofn <- function(x, ...) {
  S <- surv.exp_kofn(x)
  function(t, ...) 1 - S(t)
}


#' @rdname exp_kofn
#' @export
sampler.exp_kofn <- function(x, ...) {
  rates <- x$rates
  k <- x$k
  m <- length(rates)
  order_idx <- m - k + 1L
  function(n, ...) {
    mat <- vapply(seq_len(m), function(j) {
      stats::rexp(n, rate = rates[j])
    }, numeric(n))
    if (!is.matrix(mat)) mat <- matrix(mat, nrow = n)
    apply(mat, 1L, function(row) sort(row)[order_idx])
  }
}
