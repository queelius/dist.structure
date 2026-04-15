# ==========================================================================
# Coercions
# ==========================================================================


test_that("as_dist_structure on a plain dist returns a 1-component series", {
  d <- algebraic.dist::exponential(2)
  ds <- as_dist_structure(d)
  expect_true(is_dist_structure(ds))
  expect_s3_class(ds, "series_dist")
  expect_equal(ncomponents(ds), 1L)
  c1 <- component(ds, 1L)
  expect_s3_class(c1, "exponential")
  expect_equal(c1$rate, 2)
})


test_that("as_dist_structure on a dist_structure returns it unchanged", {
  sys <- series_dist(iid_exp_components(3))
  result <- as_dist_structure(sys)
  expect_identical(result, sys)
})


test_that("as_dist_structure preserves dist methods", {
  d <- algebraic.dist::exponential(1)
  ds <- as_dist_structure(d)
  S_d <- algebraic.dist::surv(d)
  S_ds <- algebraic.dist::surv(ds)
  for (ti in c(0.5, 1, 2)) {
    expect_equal(S_d(ti), S_ds(ti), tolerance = 1e-12)
  }
})
