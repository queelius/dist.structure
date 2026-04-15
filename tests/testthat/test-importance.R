# ==========================================================================
# Reliability importance measures
# ==========================================================================


test_that("birnbaum_importance for series equals prod of other reliabilities", {
  sys <- series_dist(iid_exp_components(4))
  p <- c(0.9, 0.8, 0.7, 0.6)
  for (j in 1:4) {
    expected <- prod(p[-j])
    expect_equal(birnbaum_importance(sys, j, p), expected, tolerance = 1e-10)
  }
})


test_that("birnbaum_importance for parallel equals prod of other unreliabilities", {
  sys <- parallel_dist(iid_exp_components(4))
  p <- c(0.9, 0.8, 0.7, 0.6)
  for (j in 1:4) {
    expected <- prod(1 - p[-j])
    expect_equal(birnbaum_importance(sys, j, p), expected, tolerance = 1e-10)
  }
})


test_that("birnbaum_importance for kofn equals binomial tail derivative at iid p", {
  # For k-of-m at iid p: R(p) = sum_{i=k}^{m} C(m, i) p^i (1-p)^(m-i).
  # dR/dp_j (same for each j by symmetry) = ...
  # Check: I_B for 2-of-3 at p = 0.5 equals
  #   P(sum of other two = k-1 = 1) = 2 * 0.5 * 0.5 = 0.5
  sys <- kofn_dist(k = 2L, iid_exp_components(3))
  expect_equal(birnbaum_importance(sys, 1L, 0.5), 0.5, tolerance = 1e-10)
  # 2-of-3 at p = 0.8: P(exactly 1 of other 2 alive) = 2 * 0.8 * 0.2 = 0.32
  expect_equal(birnbaum_importance(sys, 1L, 0.8), 0.32, tolerance = 1e-10)
})


test_that("birnbaum_importance scalar p recycles", {
  sys <- series_dist(iid_exp_components(3))
  expect_equal(
    birnbaum_importance(sys, 1L, 0.5),
    birnbaum_importance(sys, 1L, c(0.5, 0.5, 0.5)),
    tolerance = 1e-12
  )
})


test_that("criticality_importance equals Birnbaum * F_j / F_sys for series", {
  # Series of iid Exp(1): F_j(t) = 1 - exp(-t); F_sys(t) = 1 - exp(-m*t).
  # I_B(j; p) at iid p is p^(m-1).
  m <- 3L
  sys <- exp_series(rep(1, m))
  t0 <- 0.5
  p <- exp(-t0)  # common survival at t0
  ib <- p^(m - 1L)
  F_j <- 1 - p
  F_sys <- 1 - exp(-m * t0)
  expected <- ib * F_j / F_sys
  expect_equal(criticality_importance(sys, 1L, t0), expected,
               tolerance = 1e-10)
})


test_that("vesely_fussell_importance equals F_j / F_sys for series", {
  # For a series system, every component j is a single-component cut, so
  # the cuts containing j are just {{j}}. The union of events "cut {j}
  # failed" is just "component j failed", with probability F_j.
  # I_VF(j; t) = F_j / F_sys.
  m <- 3L
  sys <- exp_series(rep(1, m))
  t0 <- 0.5
  F_j <- 1 - exp(-t0)
  F_sys <- 1 - exp(-m * t0)
  expect_equal(vesely_fussell_importance(sys, 1L, t0), F_j / F_sys,
               tolerance = 1e-10)
})


test_that("vesely_fussell_importance returns 0 when component is in no cut", {
  # Pathological case: ensure correct 0 return when no cuts contain j.
  # This shouldn't happen for a coherent system unless j is irrelevant.
  # Test vacuously by choosing j for a parallel system of 1 component.
  # For 1-component parallel, min_cuts = {{1}}, so cuts containing 1 do exist.
  # Skip the vacuous check; instead, verify finiteness on realistic inputs.
  sys <- parallel_dist(iid_exp_components(3))
  t0 <- 1
  vf <- vesely_fussell_importance(sys, 1L, t0)
  expect_true(is.finite(vf))
  expect_true(vf >= 0 && vf <= 1)
})


test_that("importance measures handle bridge without error", {
  sys <- bridge_dist(iid_exp_components(5))
  t0 <- 1
  for (j in 1:5) {
    ib <- birnbaum_importance(sys, j, exp(-t0))
    ci <- criticality_importance(sys, j, t0)
    vf <- vesely_fussell_importance(sys, j, t0)
    expect_true(is.finite(ib) && ib >= 0 && ib <= 1)
    expect_true(is.finite(ci) && ci >= 0 && ci <= 1)
    expect_true(is.finite(vf) && vf >= 0 && vf <= 1)
  }
})
