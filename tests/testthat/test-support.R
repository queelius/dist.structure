test_that("closed-form system survivals respect lifetime support", {
  systems <- list(
    exp_series(c(1, 2)), exp_parallel(c(1, 2)),
    exp_kofn(2, c(1, 2, 3)),
    wei_series(c(0.5, 2), c(1, 2)),
    wei_homogeneous_series(0.5, c(1, 2)),
    wei_kofn(2, c(0.5, 1.5, 2), c(1, 2, 3)),
    gamma_series(c(1, 2), c(1, 2)),
    lognormal_series(c(0, 1), c(1, 2))
  )
  times <- c(-Inf, -2, -0.1, 0, Inf)
  for (sys in systems) {
    expect_equal(algebraic.dist::surv(sys)(times), c(1, 1, 1, 1, 0),
                 info = class(sys)[1])
    expect_equal(algebraic.dist::cdf(sys)(times), c(0, 0, 0, 0, 1),
                 info = class(sys)[1])
  }
})

test_that("closed-form densities and hazards are zero before time zero", {
  systems <- list(exp_series(c(1, 2)), exp_kofn(2, c(1, 2, 3)),
                  wei_series(c(0.5, 2), c(1, 2)),
                  wei_kofn(2, c(0.5, 1.5, 2), c(1, 2, 3)))
  for (sys in systems) {
    expect_equal(algebraic.dist::hazard(sys)(c(-Inf, -1)), c(0, 0),
                 info = class(sys)[1])
  }
  for (sys in systems[c(1, 2, 4)]) {
    expect_equal(stats::density(sys)(c(-Inf, -1)), c(0, 0),
                 info = class(sys)[1])
    expect_equal(stats::density(sys)(c(-Inf, -1), log = TRUE), c(-Inf, -Inf),
                 info = class(sys)[1])
  }
})
