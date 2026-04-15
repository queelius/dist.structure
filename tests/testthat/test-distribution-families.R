# ==========================================================================
# Additional distribution-family series specializations
# ==========================================================================


test_that("gamma_series surv equals product of per-component upper tails", {
  shapes <- c(2, 3, 1.5)
  rates <- c(1, 0.5, 2)
  sys <- gamma_series(shapes, rates)
  S <- algebraic.dist::surv(sys)
  for (ti in c(0.3, 1, 3)) {
    expected <- prod(stats::pgamma(ti, shape = shapes, rate = rates,
                                   lower.tail = FALSE))
    expect_equal(S(ti), expected, tolerance = 1e-12)
  }
})


test_that("gamma_series with shape = 1 reduces to exp_series", {
  rates <- c(0.5, 1, 2)
  sys_gamma <- gamma_series(shapes = rep(1, 3), rates = rates)
  sys_exp <- exp_series(rates)
  S1 <- algebraic.dist::surv(sys_gamma)
  S2 <- algebraic.dist::surv(sys_exp)
  for (ti in c(0.3, 1, 2)) {
    expect_equal(S1(ti), S2(ti), tolerance = 1e-10)
  }
})


test_that("gamma_series sampler empirical surv matches closed form", {
  shapes <- c(2, 3, 1.5)
  rates <- c(1, 0.5, 2)
  sys <- gamma_series(shapes, rates)
  set.seed(42)
  x <- algebraic.dist::sampler(sys)(10000L)
  for (t0 in c(0.3, 1)) {
    expect_equal(mean(x > t0), algebraic.dist::surv(sys)(t0),
                 tolerance = 0.02)
  }
})


test_that("gamma_series class hierarchy is correct", {
  sys <- gamma_series(shapes = c(2, 3), rates = c(1, 1))
  expect_s3_class(sys, "gamma_series")
  expect_s3_class(sys, "series_dist")
  expect_s3_class(sys, "dist_structure")
  expect_s3_class(sys, "dist")
})


test_that("lognormal_series surv equals product of per-component upper tails", {
  mu <- c(0, 1, -0.5)
  sd <- c(1, 0.5, 0.8)
  sys <- lognormal_series(mu, sd)
  S <- algebraic.dist::surv(sys)
  for (ti in c(0.5, 1, 3)) {
    expected <- prod(stats::plnorm(ti, meanlog = mu, sdlog = sd,
                                   lower.tail = FALSE))
    expect_equal(S(ti), expected, tolerance = 1e-12)
  }
})


test_that("lognormal_series sampler empirical surv matches closed form", {
  mu <- c(0, 1)
  sd <- c(1, 0.5)
  sys <- lognormal_series(mu, sd)
  set.seed(42)
  x <- algebraic.dist::sampler(sys)(20000L)
  for (t0 in c(0.5, 1.5)) {
    expect_equal(mean(x > t0), algebraic.dist::surv(sys)(t0),
                 tolerance = 0.03)
  }
})


test_that("lognormal_series class hierarchy is correct", {
  sys <- lognormal_series(meanlogs = c(0, 1), sdlogs = c(1, 0.5))
  expect_s3_class(sys, "lognormal_series")
  expect_s3_class(sys, "series_dist")
  expect_s3_class(sys, "dist_structure")
  expect_s3_class(sys, "dist")
})


test_that("gamma_series component(j) returns gamma_dist", {
  sys <- gamma_series(shapes = c(2, 3), rates = c(1, 2))
  c1 <- component(sys, 1L)
  expect_s3_class(c1, "gamma_dist")
  expect_equal(c1$shape, 2)
  expect_equal(c1$rate, 1)
})


test_that("lognormal_series component(j) returns lognormal", {
  sys <- lognormal_series(meanlogs = c(0, 1), sdlogs = c(1, 0.5))
  c2 <- component(sys, 2L)
  expect_s3_class(c2, "lognormal")
})
