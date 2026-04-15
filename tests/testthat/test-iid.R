# ==========================================================================
# IID convenience constructors
# ==========================================================================


test_that("min_iid equals series_dist of m copies", {
  d <- algebraic.dist::exponential(2)
  sys <- min_iid(d, m = 4L)
  expect_s3_class(sys, "series_dist")
  expect_equal(ncomponents(sys), 4L)
  # All components are the same exponential
  for (j in 1:4) {
    expect_s3_class(component(sys, j), "exponential")
  }
})


test_that("max_iid equals parallel_dist of m copies", {
  d <- algebraic.dist::exponential(1)
  sys <- max_iid(d, m = 3L)
  expect_s3_class(sys, "parallel_dist")
  expect_equal(ncomponents(sys), 3L)
})


test_that("order_statistic constructs correct kofn", {
  d <- algebraic.dist::exponential(1)
  # Median of 5: T_(3). System fails at 3rd order statistic.
  sys <- order_statistic(d, k = 3L, m = 5L)
  expect_s3_class(sys, "kofn_dist")
  expect_equal(ncomponents(sys), 5L)
  # System lifetime = T_(k) with k = 3 means 3rd smallest time
  times <- c(0.5, 0.1, 0.3, 0.7, 0.9)
  # Sorted: 0.1, 0.3, 0.5, 0.7, 0.9. T_(3) = 0.5
  expect_equal(system_lifetime(sys, times), 0.5)
})


test_that("min_iid survival matches m-th power of component survival", {
  # min of m iid Exp(1) is Exp(m)
  d <- algebraic.dist::exponential(1)
  sys <- min_iid(d, m = 3L)
  S_sys <- algebraic.dist::surv(sys)
  for (ti in c(0.2, 0.5, 1.0)) {
    expect_equal(S_sys(ti), exp(-3 * ti), tolerance = 1e-10)
  }
})
