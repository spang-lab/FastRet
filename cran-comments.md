## Submission

This is an update of FastRet from version 1.3.0 (currently on CRAN) to version 1.5.2.

- Fixes the R-devel NOTE "Found calls to structure() using deprecated special names" (`.Names` in `R/server.R` was replaced by `stats::setNames()`).

- Fixes a bug where model training with cross-validation required the suggested package toscutil.

- Adds an `inst/CITATION` file and a reference to the accompanying publication (Fadil et al., 2026, <doi:10.1021/acs.jcim.6c01344>) in the Description field.

- Includes further features, bug fixes and updated example data; see NEWS.md for details.

## Test environments

- TODO: local OS, R version

- TODO: GitHub Actions (macOS, Windows, Ubuntu; R devel, release, oldrel-1, 4.1)

- TODO: win-builder (devel, release)

## R CMD check results

TODO: 0 errors | 0 warnings | 0 notes

## Reverse dependencies

The package has no reverse dependencies.
