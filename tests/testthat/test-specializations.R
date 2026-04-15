# ==========================================================================
# Closed-form specializations: exp_series, wei_series, wei_homogeneous_series
# ==========================================================================
#
# Each specialization is tested against:
#   (a) a known analytic identity
#   (b) the general dist_structure default (correctness check)
# ==========================================================================


test_that("exp_series has class hierarchy including series_dist", {
  sys <- exp_series(c(0.5, 0.3, 0.2))
  expect_s3_class(sys, "exp_series")
  expect_s3_class(sys, "series_dist")
  expect_s3_class(sys, "dist_structure")
  expect_s3_class(sys, "dist")
  expect_equal(ncomponents(sys), 3L)
})


test_that("exp_series surv matches Exp(sum(rates))", {
  rates <- c(0.5, 0.3, 0.2)
  sys <- exp_series(rates)
  S <- algebraic.dist::surv(sys)
  lam <- sum(rates)
  for (ti in c(0.1, 0.5, 1, 2, 5)) {
    expect_equal(S(ti), exp(-lam * ti), tolerance = 1e-12)
  }
})


test_that("exp_series sampler mean matches 1 / sum(rates)", {
  rates <- c(0.5, 0.3, 0.2)
  sys <- exp_series(rates)
  set.seed(42)
  samples <- algebraic.dist::sampler(sys)(10000L)
  expect_equal(mean(samples), 1 / sum(rates), tolerance = 0.05)
})


test_that("exp_series mean is 1 / sum(rates) exactly", {
  rates <- c(1, 2, 3)
  sys <- exp_series(rates)
  expect_equal(mean(sys), 1 / 6, tolerance = 1e-12)
})


test_that("exp_series agrees with general series_dist on surv", {
  rates <- c(0.5, 0.3, 0.2)
  sys_special <- exp_series(rates)
  sys_general <- series_dist(lapply(rates, function(r)
    algebraic.dist::exponential(r)))
  S1 <- algebraic.dist::surv(sys_special)
  S2 <- algebraic.dist::surv(sys_general)
  for (ti in c(0.2, 1, 3)) {
    expect_equal(S1(ti), S2(ti), tolerance = 1e-10)
  }
})


test_that("wei_series has correct class hierarchy", {
  sys <- wei_series(shapes = c(1, 2), scales = c(1, 2))
  expect_s3_class(sys, "wei_series")
  expect_s3_class(sys, "series_dist")
  expect_s3_class(sys, "dist_structure")
})


test_that("wei_series surv matches closed-form cumulative hazard", {
  shapes <- c(1, 2, 3)
  scales <- c(1, 2, 3)
  sys <- wei_series(shapes, scales)
  S <- algebraic.dist::surv(sys)
  for (ti in c(0.3, 1, 2)) {
    expected <- exp(-sum((ti / scales)^shapes))
    expect_equal(S(ti), expected, tolerance = 1e-12)
  }
})


test_that("wei_series sampler produces values consistent with Monte Carlo surv", {
  shapes <- c(2, 2, 2)
  scales <- c(1, 2, 3)
  sys <- wei_series(shapes, scales)
  set.seed(1)
  samples <- algebraic.dist::sampler(sys)(10000L)
  t0 <- 0.5
  empirical_surv <- mean(samples > t0)
  expected_surv <- algebraic.dist::surv(sys)(t0)
  expect_equal(empirical_surv, expected_surv, tolerance = 0.02)
})


test_that("wei_series reduces to Exp when all shapes are 1", {
  # Series of Weibull(1, 1/rate) is Exp(rate); series of them is Exp(sum rates).
  scales <- c(1, 2, 4)  # corresponds to rates 1, 0.5, 0.25
  rates <- 1 / scales
  sys <- wei_series(shapes = rep(1, 3), scales = scales)
  S <- algebraic.dist::surv(sys)
  for (ti in c(0.3, 1, 2)) {
    expect_equal(S(ti), exp(-sum(rates) * ti), tolerance = 1e-12)
  }
})


test_that("wei_homogeneous_series aggregates into a single Weibull", {
  shape <- 2
  scales <- c(1, 2, 3)
  sys <- wei_homogeneous_series(shape, scales)
  agg_scale <- (sum(1 / scales^shape))^(-1 / shape)
  S <- algebraic.dist::surv(sys)
  for (ti in c(0.2, 0.7, 1.5)) {
    expected <- exp(-(ti / agg_scale)^shape)
    expect_equal(S(ti), expected, tolerance = 1e-12)
  }
})


test_that("wei_homogeneous_series mean matches Weibull formula", {
  shape <- 3
  scales <- c(1, 1.5, 2)
  sys <- wei_homogeneous_series(shape, scales)
  agg_scale <- (sum(1 / scales^shape))^(-1 / shape)
  expected <- agg_scale * gamma(1 + 1 / shape)
  expect_equal(mean(sys), expected, tolerance = 1e-12)
})


test_that("wei_homogeneous_series is a wei_series subclass and inherits methods", {
  sys <- wei_homogeneous_series(shape = 2, scales = c(1, 2))
  expect_s3_class(sys, "wei_homogeneous_series")
  expect_s3_class(sys, "wei_series")
  expect_s3_class(sys, "series_dist")
  # component() still works
  c1 <- component(sys, 1L)
  expect_s3_class(c1, "weibull_dist")
  expect_equal(c1$scale, 1)
  expect_equal(c1$shape, 2)
})


test_that("wei_homogeneous_series agrees with wei_series on survival", {
  shape <- 2
  scales <- c(1, 2, 3)
  sys_hom <- wei_homogeneous_series(shape, scales)
  sys_het <- wei_series(rep(shape, 3), scales)
  S1 <- algebraic.dist::surv(sys_hom)
  S2 <- algebraic.dist::surv(sys_het)
  for (ti in c(0.2, 1, 2)) {
    expect_equal(S1(ti), S2(ti), tolerance = 1e-12)
  }
})


test_that("exp_series is faster than general series_dist on surv (smoke test)", {
  # Not a strict benchmark, but verifies both paths produce matching values.
  rates <- runif(5, 0.1, 2)
  sys_special <- exp_series(rates)
  sys_general <- series_dist(lapply(rates, algebraic.dist::exponential))
  ts <- c(0.1, 0.5, 1, 2)
  S1 <- algebraic.dist::surv(sys_special)(ts)
  S2 <- algebraic.dist::surv(sys_general)(ts)
  expect_equal(S1, S2, tolerance = 1e-10)
})
