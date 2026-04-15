# ==========================================================================
# Parallel and k-of-n closed-form specializations
# ==========================================================================


test_that("exp_parallel surv matches the classical product-complement formula", {
  rates <- c(1, 2, 3)
  sys <- exp_parallel(rates)
  S <- algebraic.dist::surv(sys)
  for (ti in c(0.2, 0.7, 1.5, 3)) {
    expected <- 1 - prod(1 - exp(-rates * ti))
    expect_equal(S(ti), expected, tolerance = 1e-12)
  }
})


test_that("exp_parallel mean matches harmonic sum for iid", {
  # iid case with m components and common rate lambda:
  # E[max] = (1/lambda) * (1 + 1/2 + ... + 1/m)
  for (m in 2:5) {
    lambda <- 2
    sys <- exp_parallel(rep(lambda, m))
    expected <- sum(1 / seq_len(m)) / lambda
    expect_equal(mean(sys), expected, tolerance = 1e-10)
  }
})


test_that("exp_parallel mean matches sampler empirical mean", {
  rates <- c(0.5, 1, 2)
  sys <- exp_parallel(rates)
  set.seed(1)
  x <- algebraic.dist::sampler(sys)(10000)
  expect_equal(mean(x), mean(sys), tolerance = 0.05)
})


test_that("exp_parallel class hierarchy is correct", {
  sys <- exp_parallel(c(1, 2))
  expect_s3_class(sys, "exp_parallel")
  expect_s3_class(sys, "parallel_dist")
  expect_s3_class(sys, "dist_structure")
  expect_s3_class(sys, "dist")
})


test_that("exp_kofn reduces to exp_series for k = m", {
  rates <- c(0.5, 0.3, 0.2)
  sys_kofn <- exp_kofn(k = 3L, rates)
  sys_series <- exp_series(rates)
  S1 <- algebraic.dist::surv(sys_kofn)
  S2 <- algebraic.dist::surv(sys_series)
  for (ti in c(0.3, 1, 2)) {
    expect_equal(S1(ti), S2(ti), tolerance = 1e-10)
  }
})


test_that("exp_kofn reduces to exp_parallel for k = 1", {
  rates <- c(1, 2, 3)
  sys_kofn <- exp_kofn(k = 1L, rates)
  sys_par <- exp_parallel(rates)
  S1 <- algebraic.dist::surv(sys_kofn)
  S2 <- algebraic.dist::surv(sys_par)
  for (ti in c(0.2, 0.8, 2)) {
    expect_equal(S1(ti), S2(ti), tolerance = 1e-10)
  }
})


test_that("exp_kofn sampler empirical mean matches k-th order statistic (iid)", {
  # For iid Exp(1), the i-th smallest value has E[T_(i)] = sum_{j=1}^{i} 1/(m-j+1).
  # A k-of-m system fails at the (m-k+1)-th smallest time (= the i=m-k+1 order stat).
  m <- 5L; k <- 3L
  sys <- exp_kofn(k, rep(1, m))
  set.seed(1)
  x <- algebraic.dist::sampler(sys)(10000)
  i <- m - k + 1L
  expected_mean <- sum(1 / ((m - seq_len(i) + 1L)))
  expect_equal(mean(x), expected_mean, tolerance = 0.05)
})


test_that("exp_kofn class hierarchy is correct", {
  sys <- exp_kofn(k = 2L, rates = c(1, 2, 3))
  expect_s3_class(sys, "exp_kofn")
  expect_s3_class(sys, "kofn_dist")
  expect_s3_class(sys, "dist_structure")
  expect_s3_class(sys, "dist")
})


test_that("wei_kofn reduces to wei_series for k = m", {
  shapes <- c(2, 2, 2); scales <- c(1, 2, 3)
  sys_kofn <- wei_kofn(3L, shapes, scales)
  sys_series <- wei_series(shapes, scales)
  S1 <- algebraic.dist::surv(sys_kofn)
  S2 <- algebraic.dist::surv(sys_series)
  for (ti in c(0.3, 1, 2)) {
    expect_equal(S1(ti), S2(ti), tolerance = 1e-10)
  }
})


test_that("wei_kofn sampler survival probability matches closed-form surv", {
  shapes <- c(2, 2, 2); scales <- c(1, 2, 3)
  sys <- wei_kofn(2L, shapes, scales)
  set.seed(1)
  x <- algebraic.dist::sampler(sys)(10000)
  for (t0 in c(0.3, 0.8, 1.5)) {
    emp_surv <- mean(x > t0)
    analytical_surv <- algebraic.dist::surv(sys)(t0)
    expect_equal(emp_surv, analytical_surv, tolerance = 0.02)
  }
})


test_that("wei_kofn class hierarchy is correct", {
  sys <- wei_kofn(k = 2L, shapes = c(1, 2), scales = c(1, 2))
  expect_s3_class(sys, "wei_kofn")
  expect_s3_class(sys, "kofn_dist")
  expect_s3_class(sys, "dist_structure")
  expect_s3_class(sys, "dist")
})


test_that("exp_parallel and exp_kofn agree on surv at k = 1", {
  rates <- runif(4, 0.1, 2)
  sys_par <- exp_parallel(rates)
  sys_kofn <- exp_kofn(1L, rates)
  for (ti in c(0.2, 1, 3)) {
    expect_equal(
      algebraic.dist::surv(sys_par)(ti),
      algebraic.dist::surv(sys_kofn)(ti),
      tolerance = 1e-10
    )
  }
})


test_that("exp_kofn surv matches the general default (small m)", {
  rates <- c(1, 2)
  sys_special <- exp_kofn(1L, rates)
  # Compare against the general kofn_dist default
  sys_general <- kofn_dist(1L,
    lapply(rates, function(r) algebraic.dist::exponential(r)))
  for (ti in c(0.3, 1, 2)) {
    expect_equal(algebraic.dist::surv(sys_special)(ti),
                 algebraic.dist::surv(sys_general)(ti),
                 tolerance = 1e-10)
  }
})
