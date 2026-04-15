# ==========================================================================
# Dual system tests
# ==========================================================================


test_that("dual of series is parallel-like via phi", {
  sys <- series_dist(iid_exp_components(3))
  dsys <- dual(sys)
  # dual phi(state) = 1 - phi(1 - state); series.phi(1-state) = 1 iff all 1 flipped, i.e., state = 0
  expect_equal(phi(dsys, c(0, 0, 0)), 0L)
  expect_equal(phi(dsys, c(1, 0, 0)), 1L)
  expect_equal(phi(dsys, c(1, 1, 1)), 1L)
})


test_that("dual of k-of-m satisfies phi identity", {
  sys <- kofn_dist(k = 2L, iid_exp_components(4))
  dsys <- dual(sys)
  for (a in 0:1) for (b in 0:1) for (c in 0:1) for (d in 0:1) {
    state <- c(a, b, c, d)
    expect_equal(phi(dsys, state), 1L - phi(sys, 1L - state))
  }
})


test_that("dual of dual equals original (involution)", {
  sys <- kofn_dist(k = 2L, iid_exp_components(3))
  dd <- dual(dual(sys))
  # Either returns the same object (coherent_dist path via min_cuts) or a
  # double-wrapped equivalent. Check phi-equivalence.
  for (a in 0:1) for (b in 0:1) for (c in 0:1) {
    state <- c(a, b, c)
    expect_equal(phi(dd, state), phi(sys, state))
  }
})


test_that("dual_of_system inherits dist_structure defaults", {
  sys <- series_dist(iid_exp_components(3))
  dsys <- dual(sys)  # falls through to coherent_dist method producing a new coherent_dist
  expect_true(is_dist_structure(dsys))
  expect_equal(ncomponents(dsys), 3L)
})


test_that("coherent_dist dual swaps paths and cuts", {
  sys <- coherent_dist(
    min_paths = list(c(1L, 2L), c(3L)),
    components = iid_exp_components(3)
  )
  dsys <- dual(sys)
  # dual min_paths should equal original min_cuts
  expect_setequal(
    lapply(min_paths(dsys), sort),
    lapply(min_cuts(sys), sort)
  )
})
