# ==========================================================================
# Topology generics and defaults
# ==========================================================================


test_that("series_dist has correct phi and min_paths", {
  sys <- series_dist(iid_exp_components(3))
  expect_equal(ncomponents(sys), 3L)
  expect_equal(phi(sys, c(1, 1, 1)), 1L)
  expect_equal(phi(sys, c(1, 1, 0)), 0L)
  expect_equal(min_paths(sys), list(1:3))
})


test_that("parallel_dist has correct phi and min_paths", {
  sys <- parallel_dist(iid_exp_components(3))
  expect_equal(phi(sys, c(0, 0, 0)), 0L)
  expect_equal(phi(sys, c(1, 0, 0)), 1L)
  expect_equal(length(min_paths(sys)), 3L)
})


test_that("kofn_dist has correct phi", {
  sys <- kofn_dist(k = 2L, iid_exp_components(3))
  expect_equal(phi(sys, c(1, 1, 0)), 1L)
  expect_equal(phi(sys, c(1, 0, 0)), 0L)
  expect_equal(phi(sys, c(1, 1, 1)), 1L)
})


test_that("bridge_dist has expected 4 minimal path sets", {
  sys <- bridge_dist(iid_exp_components(5))
  expect_equal(length(min_paths(sys)), 4L)
  expect_equal(ncomponents(sys), 5L)
})


test_that("min_cuts defaults via Berge transversal from min_paths", {
  # series of 3: min_paths = {1,2,3}, min_cuts = {1},{2},{3}
  series_sys <- series_dist(iid_exp_components(3))
  cuts <- min_cuts(series_sys)
  expect_equal(length(cuts), 3L)
  expect_setequal(unlist(cuts), c(1L, 2L, 3L))

  # parallel of 3: min_paths = {1},{2},{3}, min_cuts = {1,2,3}
  parallel_sys <- parallel_dist(iid_exp_components(3))
  cuts <- min_cuts(parallel_sys)
  expect_equal(length(cuts), 1L)
  expect_setequal(cuts[[1]], c(1L, 2L, 3L))
})


test_that("min_cuts returns distinct sets (no duplicates) for kofn", {
  # k-of-n where k > 1 and m > k + 1 used to produce duplicates because
  # the Berge-transversal extension generates equal-length transversals
  # through independent paths. 2-of-3 has exactly 3 min cuts: all pairs.
  sys <- kofn_dist(k = 2L, iid_exp_components(3))
  cuts <- min_cuts(sys)
  expect_equal(length(cuts), 3L)
  # Verify distinct
  expect_equal(length(unique(cuts)), length(cuts))
  # Each cut is a pair (size m - k + 1 = 2)
  expect_true(all(vapply(cuts, length, integer(1L)) == 2L))
})


test_that("min_cuts of bridge has exactly 4 distinct sets", {
  # Classical bridge: cuts = {1,2}, {4,5}, {1,3,5}, {2,3,4}
  sys <- bridge_dist(iid_exp_components(5))
  cuts <- min_cuts(sys)
  expect_equal(length(cuts), 4L)
  expect_equal(length(unique(cuts)), 4L)
})


test_that("critical_states default returns expected count for k-of-m", {
  # For k-of-m, component j is critical iff exactly (k-1) of the other
  # (m-1) components are functioning.
  for (m in 3:4) for (k in 1:m) {
    sys <- kofn_dist(k, iid_exp_components(m))
    expect_equal(nrow(critical_states(sys, 1L)),
                 choose(m - 1L, k - 1L))
  }
})


test_that("system_lifetime of k-of-m is the (m - k + 1)-th order statistic", {
  sys <- kofn_dist(k = 2L, iid_exp_components(3))
  # 2-of-3 fails when m - k + 1 = 2 components have failed.
  expect_equal(system_lifetime(sys, c(5, 1, 3)), 3)
  # Series (k = m): fails at min
  expect_equal(system_lifetime(series_dist(iid_exp_components(4)),
                               c(2, 5, 1, 3)), 1)
  # Parallel (k = 1): fails at max
  expect_equal(system_lifetime(parallel_dist(iid_exp_components(4)),
                               c(2, 5, 1, 3)), 5)
})


test_that("system_censoring classifies components correctly", {
  sys <- kofn_dist(k = 2L, iid_exp_components(3))
  cens <- system_censoring(sys, c(1, 3, 5))
  expect_equal(cens$system_time, 3)
  expect_equal(cens$component_status, c("left", "exact", "right"))
})


test_that("is_coherent returns TRUE for valid systems", {
  expect_true(is_coherent(series_dist(iid_exp_components(3))))
  expect_true(is_coherent(parallel_dist(iid_exp_components(3))))
  expect_true(is_coherent(kofn_dist(2L, iid_exp_components(4))))
  expect_true(is_coherent(bridge_dist(iid_exp_components(5))))
})


test_that("structural_importance is symmetric for k-of-m and matches formula", {
  sys <- kofn_dist(k = 2L, iid_exp_components(4))
  imps <- vapply(1:4, function(j) structural_importance(sys, j), 0)
  expect_equal(imps, rep(imps[[1]], 4L))
  # For k-of-m: choose(m - 1, k - 1) / 2^(m - 1)
  expect_equal(imps[[1]], choose(3L, 1L) / 2^3)
})


test_that("reliability at iid p = 0.5 matches expected for k-of-m", {
  # For k-of-m at iid p: R(p) = P(Binomial(m, p) >= k)
  sys <- kofn_dist(k = 2L, iid_exp_components(3))
  for (p in c(0.2, 0.5, 0.7)) {
    expect_equal(
      reliability(sys, p),
      sum(stats::dbinom(2:3, size = 3, prob = p)),
      tolerance = 1e-12
    )
  }
})


test_that("reliability of series at iid p is p^m", {
  sys <- series_dist(iid_exp_components(4))
  for (p in c(0.1, 0.5, 0.9)) {
    expect_equal(reliability(sys, p), p^4, tolerance = 1e-12)
  }
})


test_that("is_dist_structure recognizes subclasses", {
  expect_true(is_dist_structure(series_dist(iid_exp_components(2))))
  expect_true(is_dist_structure(parallel_dist(iid_exp_components(2))))
  expect_true(is_dist_structure(kofn_dist(1L, iid_exp_components(2))))
  expect_false(is_dist_structure(list()))
  expect_false(is_dist_structure(algebraic.dist::exponential(1)))
})
