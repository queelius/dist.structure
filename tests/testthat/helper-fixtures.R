# ==========================================================================
# Test fixtures
# ==========================================================================

iid_exp_components <- function(m, rate = 1) {
  replicate(m, algebraic.dist::exponential(rate), simplify = FALSE)
}
