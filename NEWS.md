# FastRet 1.1.4 <!-- Commit Date: 2025-02-08 -->

- Added a cache cleanup handler that gets registered via
  `reg.finalizer()` upon package loading to ensure that the cache
  directory is removed if it doesn't contain any files that should
  persist between R sessions.

- Added an article about installation details incl. a
  troubleshooting section

- Improved function docs

- Improved examples by removing `donttest` blocks

- Improved examples & tests by using smaller example datasets to
  reduce runtime

# FastRet 1.1.3 <!-- Commit Date: 2024-06-24 -->

- Moved `patch.R` from the `R` folder to `misc/scripts`, which is
  excluded from the package build using `.Rbuildignore`. The file
  is conditionally sourced by the private function
  `start_gui_in_devmode()` if available, allowing its use during
  development without including it in the package.

- Added `\value` tags to the mentioned `.Rd` files describing the
  functions' return values.

- Added *Bonini et al. (2020) <doi:10.1021/acs.analchem.9b05765>*
  as reference to the description part of the DESCRIPTION file,
  listing  it as *Related work*. This reference is used in the
  documentation for `read_retip_hilic_data()` and `ram_cache`. No
  additional references are used in the package documentation.

- Added Fadi Fadil as a contributor. Fadi measured the example
  datasets shipped with FastRet.

- Added ORCID IDs for contributors as described in [CRAN's
  checklist for submissions].

[CRAN's checklist for submissions]:
    https://cran.r-project.org/web/packages/submission_checklist.html

# FastRet 1.1.2 <!-- Commit Date: 2024-06-18 -->

- Wrapped examples of `read_rp_xlsx()` and `read_rpadj_xlsx()`
  into `donttest` to prevent note "Examples with CPU time > 2.5
  times elapsed time: ...". By now I'm pretty sure the culprit is
  the `xlsx` package, which uses a java process for reading the
  file. Maybe we should switch to openxlsx or readxl in the
  future.

# FastRet 1.1.1 <!-- Commit Date: 2024-06-18 -->

- Improved examples of `preprocess_data()` to prevent note
  "Examples with CPU time > 2.5 times elapsed time:
  preprocess_data (CPU=2.772, elapsed=0.788)".

# FastRet 1.1.0 <!-- Commit Date: 2024-06-17 -->

- Added RAM caching to `getCDs()`

# FastRet 1.0.3 <!-- Commit Date: 2024-06-13 -->

- Added examples to `start_gui()`, `fastret_app()`, `getCDsFor1Molecule()`,
  `analyzeCDNames()`, `check_lm_suitabilitym()`, `plot_lm_suitability()`,
  `extendedTask()`, `selective_measuring()`, `train_frm()`, `adjust_frm()`,
  `get_predictors()`
- Improved lots of existing examples
- Added additionalo logging messages at various places
- Submitted to CRAN, but rejected because the following examples
  caused at least one of the following notes on the CRAN testing
  machines: (1) "CPU time > 5s", (2) "CPU time > 2.5 times elapsed
  time". In this context, "CPU time" is calculated as the sum of
  the measured "user" and "system" times.

  | function             | user  | system | elapsed | ratio |
  | -------------------- | ------| ------ | ------- | ----- |
  | check_lm_suitability | 5.667 | 0.248  | 2.211   | 2.675 |
  | predict.frm          | 2.477 | 0.112  | 0.763   | 3.393 |
  | getCDs               | 2.745 | 0.089  | 0.961   | 2.949 |

# FastRet 1.0.2 <!-- Commit Date: 2024-06-11 -->

- Improved github actions

# FastRet 1.0.1 <!-- Commit Date: 2024-06-07 -->

- Improved github actions

# FastRet 1.0.0 <!-- Commit Date: 2024-06-07 -->

Completely refactored source code, e.g.:

- Added a test suite covering all important functions

- The UI now uses Extended Tasks for background processing,
  allowing GUI usage by multiple users at the same time

- The clustering now uses Partitioning Around Medoids (PAM)
  instead of k-means, which is faster and much better suited for
  our use case

- The training of the Lasso and/or XGBoost models is no longer
  done using `caret` but using `glmnet` and `xgboost` directly.
  The new implementation is much faster and allows for full
  control over the number of workers started.

- Function `getCDs` now caches the results on Disk, making the
  retrieval of chemical descriptors much faster

- The GUI now has a console element, showing the progress of the
  background tasks like clustering and model training

- The GUI has a cleaner interface, because lots of the options are
  now hidden in the "Advanced" tab by default and are only
  displayed upon user request

# FastRet 0.99.7 <!-- Commit Date: 2023-11-30 -->

- Fixed R CMD check findings
- Fixed Github Actions
- Fixed Dockerfile

# FastRet 0.99.6 <!-- Commit Date: 2023-11-30 -->

- Added contact, privacy policy
- Fixed Dockerfile

# FastRet 0.99.4 <!-- Commit Date: 2023-11-29 -->

- Improved Github Actions and Formatting of source code

# FastRet 0.99.3 <!-- Commit Date: 2023-11-29 -->

- Reduce required R version in DESCRIPTION from 4.2 to 4.1
- Added Dockerfile
- Fixed R CMD check warnings
- Fixed R CMD check action

# FastRet 0.99.2 <!-- Commit Date: 2023-11-27 -->

- Added documentation website at:
  https://spang-lab.github.io/FastRet/

# FastRet 0.99.1 <!-- Commit Date: 2023-11-27 -->

- Initial version.

  Copy of commit cd243aa82a56df405df8060b84535633cf06b692 of
  [Christian Amesöders
  Repository](https://github.com/ChristianAmes/FastRet.git).
  (Christian wrote this initial version of FastRet as part of his
  master thesis at the Institute of functional Genomics,
  University of Regensburg).
