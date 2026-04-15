# ==========================================================================
# system_signature default
# ==========================================================================


test_that("series signature is (1, 0, ..., 0)", {
  sys <- series_dist(iid_exp_components(4))
  sig <- system_signature(sys)
  expect_equal(sig, c(1, 0, 0, 0), tolerance = 1e-12)
})


test_that("parallel signature is (0, ..., 0, 1)", {
  sys <- parallel_dist(iid_exp_components(4))
  sig <- system_signature(sys)
  expect_equal(sig, c(0, 0, 0, 1), tolerance = 1e-12)
})


test_that("signature sums to 1 for k-of-n", {
  for (m in 3:5) for (k in 1:m) {
    sys <- kofn_dist(k, iid_exp_components(m))
    sig <- system_signature(sys)
    expect_length(sig, m)
    expect_equal(sum(sig), 1, tolerance = 1e-12)
    # k-of-n fails at the (m - k + 1)-th failure, so all mass is there
    expected <- integer(m)
    expected[m - k + 1L] <- 1L
    expect_equal(sig, expected, tolerance = 1e-12)
  }
})


test_that("bridge signature is the standard (0, 1/5, 3/5, 1/5, 0)", {
  # Classical result for the bridge network.
  sys <- bridge_dist(iid_exp_components(5))
  sig <- system_signature(sys)
  expect_equal(sig, c(0, 1/5, 3/5, 1/5, 0), tolerance = 1e-12)
})
