#' dist.structure: Structured Random Variables for Reliability System
#' Distributions
#'
#' Extends the algebraic.dist distribution algebra to random variables with
#' internal structure: coherent reliability systems decomposed into components
#' arranged by a structure function (series, parallel, k-out-of-n, bridge,
#' arbitrary topologies). Every `dist_structure` object is also a `dist`, so
#' the full distribution algebra works automatically via default methods that
#' compose component-level distributions through the topology.
#'
#' @section Two algebras, one hierarchy:
#'
#' `algebraic.dist` provides the base algebra of random variables:
#' `density`, `surv`, `hazard`, `sampler`, `mean`, `vcov`, `expectation`,
#' and arithmetic on RVs (transformations).
#'
#' `dist.structure` extends it with structure-preserving operations for
#' reliability systems: component-wise decomposition, structure function
#' evaluation, path and cut sets, system signature, critical states, dual
#' structures, and Birnbaum importance measures. Because every
#' `dist_structure` is a `dist`, the base algebra inherits automatically.
#'
#' @section Primitives implementors provide:
#'
#' - [ncomponents()]: number of components `m`
#' - [component()]: return the j-th component as a `dist`
#' - One of [phi()] or [min_paths()]: the other defaults automatically
#'
#' @section Topology queries (default methods on `dist_structure`):
#'
#' [phi()], [min_paths()], [min_cuts()], [critical_states()],
#' [system_signature()], [system_lifetime()], [system_censoring()],
#' [dual()], [is_coherent()], [structural_importance()], [reliability()]
#'
#' @section Distribution queries (default methods via component composition):
#'
#' [surv()], [cdf()], [sampler()]. Further base dist methods (`mean`,
#' `vcov`, `expectation`) inherit from `univariate_dist` via numerical
#' integration of `surv`. Implementors override any default for speed.
#'
#' @section Reference implementations:
#'
#' - [coherent_dist()]: general coherent system given minimal path sets +
#'   components
#' - Topology shortcuts: [series_dist()], [parallel_dist()], [kofn_dist()],
#'   [bridge_dist()], [consecutive_k_dist()]
#' - IID convenience: [min_iid()], [max_iid()], [order_statistic()]
#'
#' Closed-form specializations (exp_series, wei_series,
#' wei_homogeneous_series) and richer importance measures
#' (birnbaum_importance, criticality_importance) are deferred to later
#' releases.
#'
#' @keywords internal
"_PACKAGE"
