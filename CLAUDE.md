# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

`dist.structure` extends the `algebraic.dist` distribution algebra to
random variables with internal structure (coherent reliability systems
decomposed into components via a structure function). Every
`dist_structure` object is a `dist`, so the base algebra inherits
automatically.

Contains (v0.5.0):
- Protocol: S3 generics + virtual base class `dist_structure` (inherits
  `dist`)
- Topology defaults on `dist_structure` (phi, min_cuts, critical_states,
  reliability, dual, is_coherent, structural_importance, system_lifetime,
  system_censoring)
- Dist defaults on `dist_structure` (sampler, surv, cdf via component +
  topology composition)
- Reference implementations: `coherent_dist`, topology shortcuts
  (`series_dist`, `parallel_dist`, `kofn_dist`, `bridge_dist`,
  `consecutive_k_dist`), iid constructors (`min_iid`, `max_iid`,
  `order_statistic`)
- Closed-form specializations: `exp_series`, `wei_series`,
  `wei_homogeneous_series`, `gamma_series`, `lognormal_series`,
  `exp_parallel`, `exp_kofn`, `wei_kofn`
- Importance measures: `structural_importance`, `birnbaum_importance`,
  `criticality_importance`, `vesely_fussell_importance`
- Compositional operations: `substitute_component`, `compose_systems`
- Coercions: `as_dist_structure`
- Non-coherent systems: `cold_standby_dist` (sum-of-lifetimes; not a
  `dist_structure`, but shares `ncomponents`/`component`)

## Protocol contract

Every implementation of a reliability system distribution:

1. Declares class containing `"dist_structure"`, `"univariate_dist"`,
   `"dist"` (plus any intermediate subclass like `"coherent_dist"`).
2. Provides `ncomponents.<cls>(x)` returning the count `m`.
3. Provides `component.<cls>(x, j, ...)` returning an `algebraic.dist::dist`.
4. Provides at least one of:
   - `phi.<cls>(x, state)` returning 0 or 1
   - `min_paths.<cls>(x)` returning a list of integer vectors
5. Optionally overrides any default for performance.

## Class hierarchy

```
dist (algebraic.dist)
  └── univariate_dist
        ├── dist_structure (virtual)
        │     └── coherent_dist
        │           ├── series_dist     ── exp_series, wei_series, wei_homogeneous_series, gamma_series, lognormal_series
        │           ├── parallel_dist   ── exp_parallel
        │           ├── kofn_dist       ── exp_kofn, wei_kofn
        │           ├── bridge_dist
        │           └── consecutive_k_dist
        │     └── dual_of_system (lazy wrapper)
        └── cold_standby_dist (NOT a dist_structure; sum-of-lifetimes)
```

`cold_standby_dist` deliberately sits outside `dist_structure` because
it has no static structure function (the active component is determined
dynamically by the failure history), but it implements `ncomponents` and
`component` so user code can iterate uniformly across system types.

## Dependency

```
algebraic.dist  →  dist.structure
```

`algebraic.dist` is the only runtime Imports. No flexhaz or
likelihood.model dependency here; those live elsewhere and can depend
on `dist.structure`.

## Commands

```bash
Rscript -e 'devtools::document()'   # regenerate NAMESPACE and man/
Rscript -e 'devtools::test()'       # run testthat
Rscript -e 'devtools::check()'      # full R CMD check
Rscript -e 'covr::package_coverage()'
```

Run a single test file:

```bash
Rscript -e 'testthat::test_file("tests/testthat/test-topology.R")'
```

## File layout

- `R/dist.structure-package.R`: package-level docs
- `R/dist_structure-class.R`: virtual base class, predicate, `format`,
  `print`
- `R/generics.R`: S3 `UseMethod` stubs with roxygen
- `R/defaults.R`: topology defaults (phi from min_paths, min_cuts via
  Berge transversal, critical_states, system_lifetime, system_censoring,
  is_coherent, structural_importance, reliability) plus shared numeric
  helpers (kofn_surv_probability, series_surv_product,
  make_component_samplers, sample_component_matrix, binary_grid,
  permutations, minimize_sets)
- `R/dist-defaults.R`: `surv`, `cdf`, `sampler` defaults composing
  component-level distributions through the topology
- `R/dual.R`: default `dual` lazy wrapper and `dual_of_system` subclass
  methods
- `R/coherent_dist.R`: general coherent system constructor plus topology
  shortcut constructors (`series_dist`, `parallel_dist`, `kofn_dist`,
  `bridge_dist`, `consecutive_k_dist`)
- `R/iid_constructors.R`: `min_iid`, `max_iid`, `order_statistic`
- `R/exp_series.R`, `R/wei_series.R`, `R/wei_homogeneous_series.R`,
  `R/gamma_series.R`, `R/lognormal_series.R`: closed-form series
  specializations
- `R/exp_parallel.R`: closed-form parallel for exponential components
- `R/exp_kofn.R`, `R/wei_kofn.R`: closed-form k-out-of-n
- `R/importance.R`: Birnbaum reliability, criticality, Vesely-Fussell
- `R/compositional.R`: `substitute_component`, `compose_systems`
- `R/cold_standby.R`: non-coherent cold-standby spare arrangement
- `R/coercions.R`: `as_dist_structure`

## Key design notes

- `dist_structure` is a virtual class inheriting `univariate_dist`:
  system lifetime is a scalar RV.
- `component(x, j, ...)` returns a fully-parameterized `dist`. For v0.1,
  components carry their own parameter values; lazy parameterization via
  extra `...` arguments is left to implementors.
- The system survival default uses the classical identity
  `S_sys(t) = R(S_1(t), ..., S_m(t))` where R is the reliability
  polynomial (multilinear extension of phi).
- `min_iid`, `max_iid`, `order_statistic` parallel the base-algebra
  operators (`min`, `max`, order statistics) but preserve topology so
  structural queries remain available. Use them when you want
  structure-aware objects; use plain `min`/`max` when you only want the
  distribution.

## Testing conventions

Test fixtures in `tests/testthat/helper-fixtures.R` define
`iid_exp_components(m)` for building quick iid-exponential component
lists used throughout the test suite.

Tests target default-method code paths (phi, min_cuts, critical_states,
etc.) as well as the dist interface (surv, sampler, cdf) and verify
statistical identities (e.g., min of iid Exp is Exp with rate sum).

## Notes

- No em-dashes in code comments or docs (repository hook enforcement).
- Package name uses a dot (`dist.structure`), matching the ecosystem
  convention (`algebraic.dist`, `likelihood.model`, `likelihood.contr`).
