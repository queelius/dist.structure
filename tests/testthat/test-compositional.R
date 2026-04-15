# ==========================================================================
# Compositional operations: substitute_component, compose_systems
# ==========================================================================


test_that("substitute_component preserves topology", {
  sys <- series_dist(iid_exp_components(3, rate = 1))
  new_comp <- algebraic.dist::exponential(5)
  sys2 <- substitute_component(sys, 2L, new_comp)
  expect_equal(ncomponents(sys2), 3L)
  # phi unchanged
  expect_equal(phi(sys2, c(1, 1, 1)), 1L)
  expect_equal(phi(sys2, c(1, 0, 1)), 0L)
  # component 2 is the replacement
  c2 <- component(sys2, 2L)
  expect_s3_class(c2, "exponential")
  expect_equal(c2$rate, 5)
})


test_that("substitute_component changes system survival", {
  sys <- exp_series(c(1, 1, 1))
  # Original: Exp(3). Replace component 2 with Exp(2). New: Exp(1+2+1) = Exp(4).
  sys2 <- substitute_component(sys, 2L, algebraic.dist::exponential(2))
  S2 <- algebraic.dist::surv(sys2)
  t0 <- 0.5
  expected <- exp(-4 * t0)
  expect_equal(S2(t0), expected, tolerance = 1e-10)
})


test_that("compose_systems of series-of-series flattens to a bigger series", {
  # A series of two series sub-systems (each with 2 components) has 4
  # components in series, min_paths = {{1,2,3,4}}.
  inner_a <- series_dist(iid_exp_components(2, rate = 1))
  inner_b <- series_dist(iid_exp_components(2, rate = 1))
  outer <- series_dist(list(inner_a, inner_b))
  composed <- compose_systems(outer, list(inner_a, inner_b))
  expect_equal(ncomponents(composed), 4L)
  paths <- min_paths(composed)
  expect_equal(length(paths), 1L)
  expect_setequal(paths[[1]], 1:4)
})


test_that("compose_systems of series-of-parallel has the right min_paths", {
  # series_dist(list(parallel_of_2, parallel_of_2)):
  #   outer min_paths = {{1, 2}} (series of 2 outer components)
  #   inner_a min_paths = {{1}, {2}}
  #   inner_b min_paths = {{1}, {2}}
  # Composed min_paths: pick one inner path for each outer component in
  # the outer path, union. Four combinations:
  #   {1} + {3} = {1, 3}
  #   {1} + {4} = {1, 4}
  #   {2} + {3} = {2, 3}
  #   {2} + {4} = {2, 4}
  inner_a <- parallel_dist(iid_exp_components(2, rate = 1))
  inner_b <- parallel_dist(iid_exp_components(2, rate = 1))
  # outer structure: series of 2 things, each a "component" representing
  # one inner
  outer_comps <- lapply(1:2, function(i) algebraic.dist::exponential(1))
  outer <- series_dist(outer_comps)
  composed <- compose_systems(outer, list(inner_a, inner_b))
  expect_equal(ncomponents(composed), 4L)
  paths <- min_paths(composed)
  expect_equal(length(paths), 4L)
  path_sorted <- lapply(paths, sort)
  expect_true(list(c(1L, 3L)) %in% lapply(path_sorted, as.integer) ||
              any(vapply(path_sorted, function(P)
                setequal(P, c(1L, 3L)), logical(1L))))
  # Check each expected path appears
  expected_sets <- list(c(1L, 3L), c(1L, 4L), c(2L, 3L), c(2L, 4L))
  for (E in expected_sets) {
    expect_true(any(vapply(path_sorted, function(P)
      setequal(P, E), logical(1L))))
  }
})


test_that("compose_systems of parallel-of-series gives one path with all components", {
  # parallel of 2 series sub-systems (each with 2 components):
  #   outer min_paths = {{1}, {2}}
  #   inner_a min_paths = {{1, 2}}
  #   inner_b min_paths = {{1, 2}}  (shifted: {{3, 4}})
  # Composed paths: for each outer path {k}, union of inner paths for k.
  #   {1}: inner_a paths = {{1, 2}}, one choice, gives {1, 2}
  #   {2}: inner_b paths = {{3, 4}}, one choice, gives {3, 4}
  inner_a <- series_dist(iid_exp_components(2, rate = 1))
  inner_b <- series_dist(iid_exp_components(2, rate = 1))
  outer_comps <- lapply(1:2, function(i) algebraic.dist::exponential(1))
  outer <- parallel_dist(outer_comps)
  composed <- compose_systems(outer, list(inner_a, inner_b))
  expect_equal(ncomponents(composed), 4L)
  paths <- min_paths(composed)
  expect_equal(length(paths), 2L)
  path_sorted <- lapply(paths, sort)
  expected_sets <- list(c(1L, 2L), c(3L, 4L))
  for (E in expected_sets) {
    expect_true(any(vapply(path_sorted, function(P)
      setequal(P, E), logical(1L))))
  }
})


test_that("compose_systems with plain dist as inner treats it as single component", {
  outer <- series_dist(lapply(1:2, function(i) algebraic.dist::exponential(1)))
  composed <- compose_systems(outer,
                              list(algebraic.dist::exponential(1),
                                   algebraic.dist::exponential(2)))
  # Single-component inner: min_paths = {{1}}; combined with series outer
  # path {{1, 2}}, composed min_paths = {{1, 2}}.
  expect_equal(ncomponents(composed), 2L)
  paths <- min_paths(composed)
  expect_equal(length(paths), 1L)
  expect_equal(sort(paths[[1]]), 1:2)
})
