# dist.structure

<!-- badges: start -->
[![R-universe status](https://queelius.r-universe.dev/badges/dist.structure)](https://queelius.r-universe.dev/dist.structure)
<!-- badges: end -->

Structured random variables for reliability system distributions.
`dist.structure` extends the `algebraic.dist` distribution algebra to
random variables with internal structure: coherent reliability systems
decomposed into components arranged by a structure function (series,
parallel, k-out-of-n, bridge, consecutive-k, and arbitrary topologies
via minimal path sets).

Every `dist_structure` object is a `dist`, so the full distribution
algebra (`mean`, `vcov`, `sampler`, `surv`, `cdf`) works automatically
via default methods that compose component-level distributions through
the topology. `dist.structure` adds structural queries (`phi`,
`min_paths`, `min_cuts`, `system_signature`, `critical_states`, `dual`,
`is_coherent`), importance measures (Birnbaum structural, Birnbaum
reliability, criticality, Vesely-Fussell), and compositional operations
(`substitute_component`, `compose_systems`).

## Installation

From r-universe:

```r
install.packages("dist.structure",
  repos = c("https://queelius.r-universe.dev", "https://cran.r-project.org"))
```

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("queelius/dist.structure")
```

## Quick tour

```r
library(dist.structure)

# A 2-out-of-3 system of independent exponentials.
sys <- exp_kofn(k = 2, rates = c(1, 2, 3))

# Distribution queries (inherited from algebraic.dist).
algebraic.dist::surv(sys)(1)           # system survival at t = 1
algebraic.dist::sampler(sys)(10)       # 10 system-lifetime samples

# Topology queries.
phi(sys, state = c(1, 1, 0))           # system functions? -> 1
min_paths(sys)                         # minimal path sets
system_signature(sys)                  # Samaniego signature

# Importance.
structural_importance(sys, j = 1)      # fraction of pivotal states
birnbaum_importance(sys, j = 1, p = 0.8)  # dR/dp_j at p = 0.8
reliability(sys, p = 0.9)              # system reliability

# Composition.
sys2 <- substitute_component(sys, j = 2,
  new_component = algebraic.dist::weibull_dist(shape = 2, scale = 1))
```

See `vignette("getting-started", package = "dist.structure")` for a
longer tour, and the Articles tab on the [pkgdown
site](https://queelius.github.io/dist.structure/) for deeper topics.

## Ecosystem

`dist.structure` sits at the bottom of the reliability ecosystem. It
depends only on `algebraic.dist`. Downstream packages implement the
protocol on their own subclasses:

- **serieshaz** registers dist.structure methods on `dfr_dist_series`
  so series systems built from arbitrary dynamic failure rate
  components get the protocol for free.
- **kofn** provides maximum likelihood estimation for k-out-of-n
  systems, using dist.structure as its DGP and topology layer.
- **maskedcauses** wraps its exponential-series API over
  `dist.structure::exp_series`.

## License

MIT. See `LICENSE`.
