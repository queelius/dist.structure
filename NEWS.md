# dist.structure 0.5.0

## Documentation

* Added six vignettes covering core workflows:
  - Getting started: five-minute tour of the package.
  - Coherent systems: phi, min_paths, min_cuts, signature, reliability
    polynomial, dual, critical states; series, parallel, k-of-n,
    bridge, and custom coherent systems.
  - Distribution interface: dist_structure inherits from dist; default
    composition of surv/cdf/sampler via the reliability polynomial;
    when to choose specializations over the general path.
  - Importance measures: structural, Birnbaum reliability, criticality,
    Vesely-Fussell; worked through the bridge network.
  - Composition and substitution: substitute_component and
    compose_systems for hierarchical construction.
  - Non-coherent systems: cold_standby_dist and its deliberate
    exclusion from the dist_structure protocol.
* Added `README.md` with quick-tour and ecosystem context.
* Added pkgdown configuration (`_pkgdown.yml`).

## Fixes

* `hazard.exp_series` is now properly registered as an S3 method
  rather than being exported as a regular function.

# dist.structure 0.4.2

* Added closed-form `density.exp_kofn` and `density.wei_kofn` methods
  via the critical-state subset enumeration formula
  `f_sys(t) = sum_j f_j(t) * P(j critical at t)`.
  Required for downstream packages (kofn) that want to compute
  exact loglik contributions without duplicating the math.

# dist.structure 0.4.1

* Added closed-form `density.exp_series` and `hazard.exp_series`
  methods so downstream packages can delegate via
  `algebraic.dist::density()` / `hazard()` rather than
  reimplementing.

# dist.structure 0.4.0

## New families

* `gamma_series(shapes, rates)`: closed-form series of Gamma components.
* `lognormal_series(meanlogs, sdlogs)`: closed-form series of
  Lognormal components.
* `cold_standby_dist(components)`: non-coherent system; deliberately
  does NOT inherit `dist_structure`. Provides exact mean (sum of
  component means) and sampler; Monte Carlo surv/cdf with per-closure
  sample caching for determinism.

## New operation

* `as_dist_structure(x)`: coercion. `dist` -> 1-component `series_dist`;
  `dist_structure` returned unchanged.

## Internal

* Added `make_component_samplers` and `series_surv_product` helpers
  (now shared across the exp/wei series and k-of-n specializations).

# dist.structure 0.3.0

## New features

* Importance measures: `birnbaum_importance`, `criticality_importance`,
  `vesely_fussell_importance`.
* Compositional operations: `substitute_component`, `compose_systems`.

## Fixes

* Critical bug in `minimize_sets`: equal-length duplicate sets were not
  being removed, causing `min_cuts` to return duplicated cut sets for
  k-of-n and bridge topologies. Fixed by deduplicating input before
  subset-minimality checks.
* `criticality_importance` now validates `j` bounds (matches other
  importance measures).

## Defaults

* Added `system_signature.dist_structure` default via `m!` permutation
  enumeration (warns for `m > 9`).

# dist.structure 0.2.0

* Added closed-form `exp_parallel`, `exp_kofn`, `wei_kofn`
  specializations.

# dist.structure 0.1.0

* Initial release. Protocol (generics + virtual base class),
  topology defaults (phi, min_cuts via Berge transversal,
  critical_states, system_lifetime, system_censoring, is_coherent,
  structural_importance, reliability), dist defaults (surv, cdf,
  sampler via component composition), reference implementations
  (coherent_dist, series/parallel/kofn/bridge/consecutive_k_dist), iid
  convenience constructors (min_iid, max_iid, order_statistic), and
  closed-form series specializations (exp_series, wei_series,
  wei_homogeneous_series).
