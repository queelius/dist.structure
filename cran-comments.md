# dist.structure 0.5.0 (initial CRAN submission)

## Summary

This is the **initial CRAN submission** of `dist.structure`.

`dist.structure` is a protocol package that extends the `algebraic.dist`
distribution algebra to structured random variables: coherent reliability
systems (series, parallel, k-out-of-n, bridge, and arbitrary topologies
specified via minimal path sets) decomposed into components arranged by
a structure function. Every `dist_structure` object is also a `dist`, so
the full distribution algebra (`mean`, `surv`, `cdf`, `sampler`, etc.)
works automatically via default methods that compose component-level
distributions through the topology. The package adds structural queries
(structure function evaluation, minimal paths and cuts, system signature,
critical states, dual structures, Birnbaum and reliability importance,
system reliability) and ships closed-form specializations for common
parametric families (`exp_series`, `wei_series`, `wei_homogeneous_series`,
`gamma_series`, `lognormal_series`, `exp_parallel`, `exp_kofn`, `wei_kofn`).

The package is the shared protocol layer underpinning a small collection
of reliability and survival-analysis packages: `serieshaz` 0.2.0,
`maskedcauses` 0.10.0, and `kofn` 0.4.0 (all developed by the same
author, currently distributed via r-universe) will each be resubmitted
to CRAN once `dist.structure` is accepted, and will then `Imports:
dist.structure`. The package's distinctive contribution relative to
existing CRAN packages is providing **one S3 protocol** (with virtual
base class, generics, defaults, and reference implementations) that
those downstream packages and any user-defined subclass can implement,
rather than each reimplementing topology and structure-function
machinery independently.

## Test environments

- Local: Ubuntu Linux 6.17, R 4.3.3
- Planned pre-submission checks: win-builder (R-devel and R-release),
  R-hub (`linux`, `macOS`, `windows-x86_64-devel`).

## R CMD check results

```
0 errors | 0 warnings | 1 note (the persistent
"checking for future file timestamps ... unable to verify current time"
NOTE caused by the local check environment being unable to reach the
network worldclock service used to verify timestamps; not a package
issue, and does not appear on win-builder).
```

The full test suite (365 testthat tests) passes; vignettes build; URLs
are clean (`urlchecker::url_check()`); spelling is clean
(`spelling::spell_check_package()`).

## Downstream dependencies

There are currently **no reverse dependencies** on CRAN; this is a new
package. Three packages by the same author depend on it locally and are
queued for resubmission once `dist.structure` is accepted:

- `serieshaz` (CRAN: 0.1.1; local: 0.2.0)
- `maskedcauses` (CRAN: 0.9.3; local: 0.10.0)
- `kofn` (not yet on CRAN; local: 0.4.0)

The maintainer will not submit any of those updates until
`dist.structure` is accepted, and will run reverse-dependency checks
locally before each submission.

## Notes for the reviewer

- The package depends only on `algebraic.dist` (already on CRAN,
  version 1.0.0) and base R (`stats`, `utils`).
- All exported functions have `\value{}` documentation. Examples
  complete in well under the 5s budget; none use the network or
  filesystem.
- The package is MIT-licensed; both `LICENSE` (the CRAN-required
  template) and `LICENSE.md` (the human-readable form, listed in
  `.Rbuildignore`) are included.
- Tests use `withr::local_seed()` to avoid mutating the user's
  global RNG state.
