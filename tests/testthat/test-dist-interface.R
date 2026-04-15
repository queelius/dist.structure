# ==========================================================================
# Distribution interface: dist_structure IS-A dist
# ==========================================================================


test_that("series_dist objects are dists", {
  sys <- series_dist(iid_exp_components(3))
  expect_true(algebraic.dist::is_dist(sys))
  expect_s3_class(sys, "dist_structure")
  expect_s3_class(sys, "dist")
})


test_that("surv(series_dist) equals product of component survs", {
  # Series of 3 iid Exp(1): S_sys(t) = exp(-3t)
  sys <- series_dist(iid_exp_components(3, rate = 1))
  S <- algebraic.dist::surv(sys)
  for (ti in c(0.1, 0.5, 1.0, 2.0)) {
    expect_equal(S(ti), exp(-3 * ti), tolerance = 1e-10)
  }
})


test_that("surv(parallel_dist) equals 1 - product of component cdfs", {
  # Parallel of 3 iid Exp(1): F_sys(t) = (1 - exp(-t))^3, S = 1 - F
  sys <- parallel_dist(iid_exp_components(3, rate = 1))
  S <- algebraic.dist::surv(sys)
  for (ti in c(0.1, 0.5, 1.0, 2.0)) {
    expect_equal(S(ti), 1 - (1 - exp(-ti))^3, tolerance = 1e-10)
  }
})


test_that("cdf equals 1 - surv", {
  sys <- series_dist(iid_exp_components(2))
  S <- algebraic.dist::surv(sys)
  F_fn <- algebraic.dist::cdf(sys)
  for (ti in c(0.3, 0.7, 1.5)) {
    expect_equal(F_fn(ti), 1 - S(ti), tolerance = 1e-12)
  }
})


test_that("sampler returns the right number of samples", {
  sys <- series_dist(iid_exp_components(3, rate = 1))
  samp <- algebraic.dist::sampler(sys)
  set.seed(1)
  x <- samp(100)
  expect_length(x, 100L)
  expect_true(all(x >= 0))
})


test_that("sampler of series iid Exp(rate) has approximately the right mean", {
  # Series of m iid Exp(rate) is Exp(m * rate). Mean = 1/(m*rate).
  m <- 3L
  rate <- 2
  sys <- series_dist(iid_exp_components(m, rate = rate))
  set.seed(42)
  x <- algebraic.dist::sampler(sys)(10000)
  expected_mean <- 1 / (m * rate)
  expect_equal(mean(x), expected_mean, tolerance = 0.05)
})


test_that("sampler of parallel iid has mean close to harmonic-series formula", {
  # Max of m iid Exp(1): E[max] = 1 + 1/2 + ... + 1/m
  m <- 4L
  sys <- parallel_dist(iid_exp_components(m, rate = 1))
  set.seed(42)
  x <- algebraic.dist::sampler(sys)(10000)
  expected_mean <- sum(1 / (1:m))
  expect_equal(mean(x), expected_mean, tolerance = 0.05)
})


test_that("surv consistent with reliability identity S_sys(t) = R(S(t))", {
  sys <- bridge_dist(iid_exp_components(5, rate = 1))
  S_sys <- algebraic.dist::surv(sys)
  for (ti in c(0.2, 0.8, 1.5)) {
    # Manual computation
    p <- exp(-ti)
    expected <- reliability(sys, p)
    expect_equal(S_sys(ti), expected, tolerance = 1e-10)
  }
})
