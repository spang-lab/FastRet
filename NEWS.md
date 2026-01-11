
# FastRet 1.3.2 <!-- Commit Date: 2026-01-11 -->

Bugfix:

1. `adjust_frm()` automatically switched to "lm" adjustment if `predictors`
   contained only one predictor, regardless of whether `add_cds` was TRUE or
   FALSE. This has now been fixed, so that "lasso", "ridge" and "gbtree"
   adjustment is now possible if either `add_cds` is TRUE, `add_cds` is
   `NULL` or `predictors` contains more than one predictor.

# FastRet 1.3.1 <!-- Commit Date: 2026-01-11 -->

API Improvements:

1. Added arguments `match_rts` and `match_keys` to `adjust_frm()`:
   - If `match_rts=TRUE` (default), RTs are obtained by matching rows in
     `new_data` to rows in `frm$df` based on `match_keys`.
   - If `match_rts=FALSE`, RTs are obtained by applying the base model to
     `new_data`.
   - `match_keys` can be any combination of "INCHIKEY", "SMILES" and "NAME". If
     left at default NULL, SMILES+INCHIKEY is used if both columns are present
     in the adjusted and the original training data. Otherwise, SMILES+NAME is
     used.

# FastRet 1.3.0 <!-- Commit Date: 2025-11-12 -->

API Improvements:

1. `getCDs()`:
   - Calculation of CDs is skipped, if all CDs are already present in the input
     dataframe.
   - In addition to a dataframe with column "SMILES", a plain character vector
     of SMILES strings can now be provided as input.
   - Improved progress output.
   - More Smiles are now pre-cached internally to speed up retrieval of CDs for
     large datasets.
2. `plot_frm()`:
   - Now available as exported function (the function existed before, but was
    only exposed via the Graphical User Interface). Now users can call it
     directly from R scripts.
   - Added semi-transparent background to legends for better readability
   - Changed the point character from unfilled circles with colored borders to
     filled circles with black borders for better visibility.
   - Constant predictions (causing correlation to be NA) do no longer cause
     a warning.
3. `preprocess_data()`:
   - Added argument `add_cds` to control whether chemical descriptors should be
     added to the input data using `getCDs()`.
   - Added argument `rm_ucs` to control whether unsupported columns (i.e.
     columns that are neither mandatory nor optional) should be removed from the
     input data.
   - Added argument `rt_terms` to control whether transformations of the RT column
     (square, cube, log, exp, sqrt) should be added to the input data.
   - In case of missing mandatory columns (SMILES, RT, NAME) an error is now raised.
   - INCHIKEY and Chemical Descriptors listed in `CDFeatures` are allowed as
     optional columns.
   - Columns that are neither mandatory nor optional are automatically removed.
   - Improved runtime for generation of polynomial features and/or interaction
     terms.
   - Removal of NA and/or near-zero-variance predictors is now done after adding
     polynomial features and/or interaction terms.
   - Improved progress output.
4. `train_frm()`:
   - Added argument `do_cv` to control whether cross-validation should be
     performed for performance estimation. Default is TRUE.
   - Removal of near-zero-variance predictors and/or removal of NA values is now
     done as part of the internal model training, i.e. it happens separately for
     each fold during cross-validation. This prevents data leakage from the
     training set to the validation set. The corresponding hint about
     "overoptimistic cross-validation results" has consequently been removed
     from the documentation.
   - Argument `method` now accepts two values for training models with xgbtree
     base: "gbtreeDefault" (train xgboost with default params) and "gbtreeRP"
     (train xgboost with parameters optimized for the RP dataset). The old value
     "gbtree" still works and is now an alias for "gbtreeDefault".
   - Improved documentation of the return value (i.e. `frm` objects are fully
     specified now).
   - Added type checking for each user input.
   - Performance estimation via cross-validation now uses the new clipping
     mechanism provided by `clip_predictions()`. Of course, the clipping is
     always based on the RT range of training folds, not the whole original
     training data.
   - Calculating proportion of variance explained (R²) no longer throws a
     warning for constant predictions (causing correlation to be NA). Instead,
     R² is set to 0 in such cases.
5. `predict.frm()`:
   - Data transformations applied to the training data (adding polynomial
     features and/or adding interaction terms) are now automatically applied to
     new data as well. This was not the case before, leading to errors when
     trying to predict RTs using models trained with `degree_polynomial>1`
     and/or `interaction_terms=TRUE`, unless the transformations were manually
     applied to the new data beforehand.
   - Added argument `clip` to allow clipping of predictions to be within
     the RT range of the training data. Works for both adjusted and unadjusted
     models.
   - Predictions are now clipped to be within a sensible range by default. To
     produce unclipped predictions, set `clip=FALSE`. See `clip_predictions()`
     for details.
   - If a chemical descriptor is NA in the new data, but required by the model,
     the NA value is now replaced by the mean of that descriptor in the training
     data. Previously, the prediction for such entries was NA. The old behavior
     can be restored by setting `impute=FALSE` in `predict.frm()`.
6. `selective_measuring()`:
   - Added argument `rt_coef`, allowing user to control the influence of RT on
     the clustering. A value of 0 means that RT is ignored, a value of
     "max_ridge_coefficient" means that RT has the same weight as the most
     important chemical descriptor and a value of 1 means no scaling at all
     (except standardization to z-scores, which is applied before to the whole
     dataset before the ridge regression is trained).
7. `adjust_frm()`:
   - Added argument `seed` to allow reproducible results.
   - Added argument `do_cv` to control whether cross-validation should be
     performed for performance estimation. Default is TRUE.
   - Added argument `adj_type` to control which model should be trained for
     adjustment: supported options are "lm", "lasso", "ridge", or "gbtree".
     Previously, only "lm" was supported. To stay backwards compatible, the default
     is "lm".
   - Added argument `add_cds` to control whether chemical descriptors should be
     added to the input data using `getCDs()`. Only recommended for adj_type other
     than "lm".
   - Added support for mapping by SMILES+INCHIKEY in addition to SMILES+NAME.
     SMILES+INCHIKEY is used by default if both columns are present in the
     adjusted and the original training data. Otherwise SMILES+NAME is used as
     before.
   - Improved error handling. Previously, unmappable entries in the new data had
     been ignored silently. Now, an error is raised in such cases.
   - Function arguments are now stored in the returned frm object for better
     reproducibility.
   - Mapping is now performed by matching each new entry to the average RT of
     all original training entries with the same key (SMILES+INCHIKEY or
     SMILES+NAME). Example: if the new dataset contains a key twice, and the
     original training data contains the key three times, both new entries are
     mapped to the average RT of the three original entries.
   - Performance estimation via cross-validation now uses the new clipping
     mechanism provided by `clip_predictions()`. Of course, the clipping is
     always based on the RT range of training folds, not the whole original
     training data.
   - Calculating proportion of variance explained (R²) no longer throws a
     warning for constant predictions (causing correlation to be NA). Instead,
     R² is set to 0 in such cases.
8. `print.frm()`:
   - frm objects can now be printed directly to the console in a user-friendly
     format.
9. `clip_predictions()`:
   - New utility function for clipping predicted RTs to be within a sensible
     range. Used internally by `train_frm()`, `predict.frm()` and
     `adjust_frm()`.
10. `get_predictors()`:
    - Added arguments `base` and `adjust` to control whether predictors for the
      base model, the adjustment model or both should be returned.

Bugfixes:

1. Interaction terms generated by `preprocess_data()` are now generated
   correctly as product of the involved features instead of a division. This
   follows common practice in regression modeling and avoids division by zero
   issues. Passing older models, trained with division-based interaction terms,
   to downstream functions like `predict.frm()` or `adjust_frm()` will now lead
   to an error. (This is not a breaking change, as `predict.frm()` and friends
   have in fact never been able to handle such models).
2. `plot_frm()` with type "scatter.cv.adj" or "scatter.train.adj" now correctly
   shows retention times from the new data (used for model adjustment) as x-axis
   values instead of the original training retention times.
3. `catf()` now only emits escape codes (i.e. colored output), it the output is
   directed to a terminal. If the output is redirected to a file or a pipe, no
   escape codes are emitted anymore. Since `catf()` is used throughout the
   package for logging, this fixes the output for the whole package.

Internal Improvements:

1. Added or improved unit tests for:
   - `adjust_frm()`
   - `fit_gbtree()`
   - `fit_glmnet()`
   - `get_param_grid()`
   - `get_predictors()`
   - `getCDs()`
   - `plot_frm()`
   - `predict_frm()`
   - `preprocess_data()`
   - `selective_measuring()`
   - `train_frm()`
   - `validate_inputdata()`
2. Removed `caret` dependency by adding custom implementations for:
   - `createFolds()`
   - `nearZeroVar()`
3. Extract mapping and merging part of `adjust_frm()` into a private function
   `merge_dfs()`.
4. Replaced `fit_glmnet()`, `fit_lasso()` and `fit_ridge()` with a single
   function `fit_glmnet()`, that takes the method ("lasso" or "ridge") as
   parameter. Instead of a dataframe `df` that has to contain only predictors
   plus the RT column (as reponse), the function now takes a matrix of
   predictors `X` and a vector of responses `y`. This makes the function more
   flexible and easier to test.
5. Replaced `fit_gbtree_grid()` with a much simpler function
   `find_params_best()`. Instead of allowing the specification of every grid
   parameter, the new function instead accepts a keyword `searchspace` for
   specifying predefined grids to choose from.
6. Improved `fit_gbtree` by exposing lots of hardcoded internal xgboost
   parameters as function parameters with sensible defaults. In particular, the
   user can now set `xpar` to "default", "rpopt" or a predefined grid-size to
   train the model with different hyperparameter settings. Furthermore, the
   function is now written in a way that works with both, version 1.7.9.1 and
   the new 3.1.2.1 version published on 2025/12/03 (yes, version 2.x was skipped
   completely).
7. Added helper function `get_param_grid()` for returning predefined
   hyperparameter grids for xgboost model training based on keywords like
   "tiny", "small" or "large".
8. Added function `benchmark_find_params()` to benchmark runtime of
   `find_params_best()` for different numbers of cores and/or threads. As it
   turns out, choosing a higher number of cores is usually more efficient (at
   the cost of worse progress output).
9. Added utility functions `named()`, `as_str()`, `is_valid_smiles()` and
   `as_canonical()`

# FastRet 1.2.2 <!-- Commit Date: 2025-11-05 -->

- Improved `selective_measuring()` by aligning glmnet coefficients to columns by
  name (more stable) and by including RT, scaled by `max(abs(coefs))`,  in PAM
  clustering.
- Added `libwebp-dev` as dependency to Dockerfile.

# FastRet 1.2.1 <!-- Branch: 2025-09-23 -->

- Add updated Measurements `Measurements_v8.xlsx` to `inst/extdata/`. The new
  list contains corrections to the old `RP` dataset plus 1660 new measurements
  measured on a total of 18 different chromatographic environments.
- Reintroduced RAM caching (although hugely simplified).

# FastRet 1.2.0 <!-- Commit Date: 2025-09-22 -->

- Added `seed` parameter to `selective_measuring()` function for reproducible
  clustering results
- Enhanced documentation for `train_frm()` function
- Removed `digest` and `shinybusy` dependencies
- Major refactoring of caching system and related functions
- Removed mock files from `inst/mockdata/`
- Removed objects: `getCDsFor1Molecule()`, `get_cache_dir()`, `ram_cache` (these
  were exported, but declared as internal)
- Added private function `parLapply2`
- Added comprehensive GitHub Copilot instructions file
- Improved code organization and documentation across multiple R files

# FastRet 1.1.5 <!-- Commit Date: 2025-06-26 -->

- Improved `read_retip_hilic_data()`:
  the dataset is now only downloaded from GitHub if the package is not installed.
  If it is installed, the dataset is loaded directly.

- Internal Changes:
  - Removed `TODOS.md`
  - Bumped version to 1.1.5
  - Moved all data related functions from `util.R` to `data.R`
  - Added a README to `misc/datasets`
  - Added functions `load_all()` and `document()` to `util.R`
  - Replaced `xlsx` and `readxl` packages with `openxlsx`

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
- Added additional logging messages at various places
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

  Copy of commit `cd243aa82a56df405df8060b84535633cf06b692` of [Christian
  Amesöders Repository](https://github.com/ChristianAmes/FastRet).
  (Christian wrote this initial version of FastRet as part of his master thesis
  at the Institute of functional Genomics, University of Regensburg).

