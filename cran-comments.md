# dist.structure 0.5.1

## Changes

This patch corrects evaluation outside the non-negative lifetime support in
closed-form exponential and Weibull systems. Survival is now one, and CDF,
density and hazard are zero, for negative times. Component probabilities use
base R distribution functions. Regression tests cover eight system families.

## Test environment and results

Local Ubuntu 24.04.3 LTS, R 4.3.3, x86_64 Linux, checked 2026-10-28:

R CMD check --as-cran: 0 errors, 0 warnings, 2 notes.

* Unable to verify current time (timestamp check).
* HTML validation skipped because the system command 'tidy' is unavailable.

CRAN incoming/dependency checks completed. Package tests, examples, vignette
rebuilds and the PDF manual passed. Full development tests also passed, with
84.86% line coverage. Integration checks loaded the modified dist.structure
and flexhaz source packages together, then exercised serieshaz fitting and
interoperability and maskedhaz fitting and cross-validation.

These are local Linux results; no new cross-platform checks have been run
for this patch.
