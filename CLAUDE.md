# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

`dist.structure` extends the `algebraic.dist` distribution algebra to
random variables with internal structure (coherent reliability systems
decomposed into components via a structure function). Every
`dist_structure` object is a `dist`, so the base algebra inherits
automatically.

Contains:
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

Deferred to future versions:
- v0.2: closed-form specializations (`exp_series`, `wei_series`,
  `wei_homogeneous_series`)
- v0.3: richer importance measures (Birnbaum reliability importance,
  criticality, Vesely-Fussell), compositional operations
  (compose_systems, substitute_component), coercions to algebraic.dist
  base operations

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
        └── dist_structure (virtual)
              └── coherent_dist
                    ├── series_dist
                    ├── parallel_dist
                    ├── kofn_dist
                    ├── bridge_dist
                    └── consecutive_k_dist
              └── dual_of_system (lazy wrapper)
```

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
  is_coherent, structural_importance, reliability)
- `R/dist-defaults.R`: `surv`, `cdf`, `sampler` defaults composing
  component-level distributions through the topology
- `R/dual.R`: default `dual` lazy wrapper and `dual_of_system` subclass
  methods
- `R/coherent_dist.R`: general coherent system constructor plus topology
  shortcut constructors (`series_dist`, `parallel_dist`, `kofn_dist`,
  `bridge_dist`, `consecutive_k_dist`)
- `R/iid_constructors.R`: `min_iid`, `max_iid`, `order_statistic`

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
