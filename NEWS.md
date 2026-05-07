# dist.structure 0.5.0

Initial CRAN submission.

## New exports

* `validate_dist_structure(x)`: helper for subclass implementors.
  Checks the class chain, presence of `ncomponents()` / `component()`,
  and (via `getS3method` lookup) whether at least one of `phi.X` or
  `min_paths.X` is registered. Fails fast at construction time
  instead of at first method dispatch.
* `min_paths.dist_structure`: default derived by enumeration from
  `phi`. Closes a protocol gap: the contract said "provide phi or
  min_paths" but only the phi-from-paths direction had a default;
  users providing only `phi` previously hit "no applicable method"
  on `min_paths()` and any default that depended on it
  (`min_cuts`, `system_signature`, `reliability`, ...).
* `surv.series_dist`: product of component survival functions.
  Avoids the `2^m` reliability-polynomial enumeration in the
  `dist_structure` default for series systems whose components are
  not in a closed-form family.
* `surv.parallel_dist`: `1 - prod(component CDFs)`. Same rationale
  as `surv.series_dist`.
* `hazard.wei_series`: closed-form additive Weibull hazard
  `sum_j (k_j/s_j) * (t/s_j)^(k_j - 1)`.
* `hazard.wei_kofn`, `hazard.exp_kofn`: composition of the existing
  closed-form `density()` and `surv()`.
* `component.dual_of_system`: delegates to the underlying original.
  User-defined `dist_structure` subclasses that route `dual()`
  through the lazy-wrapper path now compose with `surv()`,
  `sampler()`, and other defaults.

## Documentation

* Added seven vignettes covering core workflows:
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
  - Implementing a dist_structure subclass: 4-step recipe with a
    worked alarm-system example.
* Added `README.md` with quick-tour and ecosystem context.
* Added pkgdown configuration (`_pkgdown.yml`).
* Per-method `@return` blocks across `dist_structure` and the
  closed-form family topics.
* `kofn_dist`: explicit `:G` vs `:F` convention note with
  `@seealso` to `order_statistic`.
* `bridge_dist`: Barlow-Proschan reference and component-role
  description.

## Fixes

* `hazard.exp_series` is now properly registered as an S3 method
  rather than being exported as a regular function.
* `system_signature` guards against always-zero `phi` (signature is
  undefined for non-functioning systems; previously silently
  returned `(1, 0, ..., 0)`).
* `mean.cold_standby_dist` falls back to Monte Carlo when component
  means silently propagate `NA` instead of returning `NA`.
* `coherent_dist` (and inherited by all topology shortcuts) rejects
  `NULL` or non-`dist` components with a clear message instead of
  letting the error surface deep in default-method dispatch.

## Internal cleanup

* Split `R/defaults.R`: shared numeric helpers move to
  `R/internal-utils.R`; `defaults.R` now contains only the S3
  default methods (165 lines down from 318).
* Removed seven redundant `cdf` overrides; `cdf.dist_structure`
  default produces an identical closure via dispatch.
* `exp_parallel`: flipped so `surv` is primary, matching the
  family convention.
* Removed redundant `min_paths.kofn_dist` (the parent
  `coherent_dist` already caches `combn` output as `$min_paths`).
* Test seeds migrated to `withr::local_seed()` (auto-restore on
  `test_that` exit); `withr` added to `Suggests`.

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
