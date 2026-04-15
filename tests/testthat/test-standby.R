# ==========================================================================
# Cold-standby system
# ==========================================================================


test_that("cold_standby_dist mean equals sum of component means", {
  comps <- list(
    algebraic.dist::exponential(1),     # mean 1
    algebraic.dist::exponential(2),     # mean 0.5
    algebraic.dist::gamma_dist(2, 1)    # mean 2
  )
  sys <- cold_standby_dist(comps)
  expect_equal(mean(sys), 1 + 0.5 + 2, tolerance = 1e-10)
})


test_that("cold_standby_dist of m iid Exp(rate) sampler ~ Gamma(m, rate)", {
  m <- 3L
  rate <- 2
  sys <- cold_standby_dist(replicate(m,
    algebraic.dist::exponential(rate), simplify = FALSE))
  set.seed(1)
  x <- algebraic.dist::sampler(sys)(10000L)
  # E[Gamma(m, rate)] = m / rate
  expect_equal(mean(x), m / rate, tolerance = 0.05)
  # Var = m / rate^2
  expect_equal(var(x), m / rate^2, tolerance = 0.1)
})


test_that("cold_standby_dist surv via Monte Carlo approximates Gamma surv", {
  m <- 4L
  rate <- 1
  sys <- cold_standby_dist(replicate(m,
    algebraic.dist::exponential(rate), simplify = FALSE))
  set.seed(1)
  S_mc <- algebraic.dist::surv(sys)
  for (ti in c(2, 4, 6)) {
    mc_estimate <- S_mc(ti, mc = 50000L)
    exact <- stats::pgamma(ti, shape = m, rate = rate, lower.tail = FALSE)
    expect_equal(mc_estimate, exact, tolerance = 0.02)
  }
})


test_that("cold_standby_dist surv closure caches samples (deterministic across calls)", {
  sys <- cold_standby_dist(replicate(3,
    algebraic.dist::exponential(1), simplify = FALSE))
  set.seed(42)
  S <- algebraic.dist::surv(sys)
  v1 <- S(2, mc = 1000L)
  v2 <- S(2, mc = 1000L)  # same mc -> uses cached samples
  expect_identical(v1, v2)
  v3 <- S(c(1, 2, 3), mc = 1000L)
  expect_identical(v3[2], v2)
})


test_that("cold_standby_dist surv regenerates when mc changes", {
  sys <- cold_standby_dist(replicate(3,
    algebraic.dist::exponential(1), simplify = FALSE))
  set.seed(42)
  S <- algebraic.dist::surv(sys)
  v_small <- S(2, mc = 1000L)
  v_large <- S(2, mc = 5000L)
  # different sample counts means different cache; values may differ
  expect_true(is.numeric(v_small) && is.numeric(v_large))
})


test_that("cold_standby_dist mean falls back to MC for components without exact mean", {
  # A series_dist of 2 exponentials does not provide a `mean()` method;
  # mean() routes through univariate_dist::mean -> expectation -> sup,
  # which is not implemented for dist_structure. cold_standby_dist's
  # mean.cold_standby_dist must catch this and fall back to MC.
  inner <- series_dist(replicate(2,
    algebraic.dist::exponential(1), simplify = FALSE))
  sys <- cold_standby_dist(list(inner, algebraic.dist::exponential(1)))
  set.seed(1)
  m <- mean(sys)
  # Inner is min of 2 iid Exp(1) = Exp(2), mean 0.5; outer Exp(1) mean 1.
  # Total expected ~ 1.5; MC tolerance ~5%.
  expect_equal(m, 1.5, tolerance = 0.1)
})


test_that("cold_standby_dist class does NOT include dist_structure", {
  sys <- cold_standby_dist(list(algebraic.dist::exponential(1)))
  expect_false(is_dist_structure(sys))
  expect_s3_class(sys, "cold_standby_dist")
  expect_s3_class(sys, "dist")
})


test_that("cold_standby_dist exposes ncomponents and component", {
  comps <- list(algebraic.dist::exponential(1),
                algebraic.dist::gamma_dist(2, 1))
  sys <- cold_standby_dist(comps)
  expect_equal(ncomponents(sys), 2L)
  expect_s3_class(component(sys, 1L), "exponential")
  expect_s3_class(component(sys, 2L), "gamma_dist")
})


test_that("cold_standby_dist sampler returns positive values of correct length", {
  sys <- cold_standby_dist(list(
    algebraic.dist::exponential(1),
    algebraic.dist::exponential(1)
  ))
  set.seed(1)
  x <- algebraic.dist::sampler(sys)(50L)
  expect_length(x, 50L)
  expect_true(all(x > 0))
})


test_that("criticality_importance validates j is in range", {
  sys <- series_dist(replicate(3,
    algebraic.dist::exponential(1), simplify = FALSE))
  expect_error(criticality_importance(sys, 0L, 0.5))
  expect_error(criticality_importance(sys, 4L, 0.5))
  expect_error(criticality_importance(sys, c(1L, 2L), 0.5))
})
