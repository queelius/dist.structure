# ==========================================================================
# Test fixtures
# ==========================================================================

exp_component <- function(rate = 1) algebraic.dist::exponential(rate)

iid_exp_components <- function(m, rate = 1) {
  replicate(m, exp_component(rate), simplify = FALSE)
}
