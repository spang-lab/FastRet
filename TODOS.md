# Todos

- [Open](#open)
  - [Make sure tests run on Linux and MacOS](#make-sure-tests-run-on-linux-and-macos)
  - [Update included Measurements](#update-included-measurements)
  - [Add tests for the shiny GUI](#add-tests-for-the-shiny-gui)
  - [Remove RP](#remove-rp)
  - [Ensure test coverage](#ensure-test-coverage)
  - [Ensure backwards compatibility](#ensure-backwards-compatibility)
  - [Resubmit to CRAN](#resubmit-to-cran)
- [Done](#done)
  - [Allow disabling of cross-validation](#allow-disabling-of-cross-validation)
  - [Fix failing tests in Github CI](#fix-failing-tests-in-github-ci)
  - [Allow Adjustment based on chemical descriptors](#allow-adjustment-based-on-chemical-descriptors)
  - [Add gbtree and glmnet as adjustment models](#add-gbtree-and-glmnet-as-adjustment-models)

## Open

### Make sure tests run on Linux and MacOS

Currently failing tests:

```log
══ Failed tests ════════════════════════════════════════════════════════════════
── Error ('test-fit_gbtree.R:47:5'): fit_gbtree works with xgboost 1.7 ─────────
<callr_status_error/callr_error/rlib_error_3_0/rlib_error/error/condition>
Error: Error: ! in callr subprocess.
Caused by error:
! Expected `packageVersion("xgboost", lib.loc = new_lib) == package_version("1.7.9.1")` to be TRUE.
Differences:
`actual`:   FALSE
`expected`: TRUE

Installed xgboost version: 1.7.11.1

[ FAIL 1 | WARN 0 | SKIP 1 | PASS 94 ]
Error:
! Test failures.
Execution halted
```


### Update included Measurements

- Remove `inst/extdata/Measurements_v8.xlsx`
- Add `data/datasets.rda`, containing datasets described in
  [../freda/README.md](../freda/README.md), i.e., HILIC_1, HILIC_2,
  HILIC_Retip_1, HILIC_Retip_2, RP_1, RP_2, RP_AXMM_1, RP_AXMM_2, RP_B_1,
  RP_B_2, RP_Flat, RP_FR25_Flat, RP_Steep, RP_T25_Flat, RP_T25_FR25_Flat,
  RP_T25_FR25_Steep, RP_Val, RP_Val_Flat, RP_Val_FR25_Flat, RP_Val_Steep,
  RP_Val_T25_Flat, RP_Val_T25_FR25_Flat, RP_Val_T25_FR25_Steep.

### Add tests for the shiny GUI

Do some research about best-practices for testing shiny applications.
If no good options exists, define a minimal set of tests that ensure that the
shiny application works as expected. And then execute each test manually.

### Remove RP

Replace all mentions of the following objects with corresponding objects from
the new `datasets.rda` object. See issue 'Update included Measurements'.

- `RP.rda`
- `RP.xlsx`
- `RP_adj.xlsx`
- `RP_lass_model.rds`
- `Measurements_v8.xlsx`

Then mark the all old objects as deprecated in the documentation

### Ensure test coverage

Ensure we have tests for the following scenarios:

1. Train a ridge model, predict with it, adjust with lm/RT, predict with adjusted model
2. Train a lasso model, predict with it, adjust with gbtree/RTTCD, predict with adjusted model
3. Train a gbtree model, predict with it, adjust with lasso/RTCD, predict with adjusted model

### Ensure backwards compatibility

Include an adjusted lasso and gbtree model that was trained and adjusted with
FastRet version 1.2.2 or earlier. Write testcase to ensure that all *.frm
functions like adjust_frm, predict.frm, get_predictors, coef.frm, plot.frm work
as expected with those.

### Resubmit to CRAN

Resubmit FastRet v1.3.0 to CRAN.

## Done

### Allow disabling of cross-validation

In `adjust_frm` allow `docv = FALSE` (do-cross-validation). In `train_frm` allow
`docv = FALSE` (do-cross-validation). The corresponding `cv` element of the
returned object should be `NULL` in those cases. The docs and tests must be
adjusted accordingly.

### Fix failing tests in Github CI

The following tests:

```txt
── Error ('test-train_frm-gbtree.R:5:5'): train_frm works if `method == "GBTree"` ──
<subscriptOutOfBoundsError/error/condition>
Error in `FUN(X[[i]], ...)`: subscript out of bounds
Backtrace:
1. └─FastRet::train_frm(...) at test-train_frm-gbtree.R:5:5
2.   └─base::lapply(tmp, "[[", 2)

── Error ('test-fit_gbtree.R:8:5'): fit.gbtrees works as expected ──────────────
Error in `begin_iteration:end_iteration`: argument of length 0
Backtrace:
1. └─FastRet:::fit_gbtree(df, verbose = 0) at test-fit_gbtree.R:8:5
2.   └─FastRet:::fit_gbtree_grid(...)
3.     └─xgboost::xgb.train(...)

── Error ('test-fit_gbtree.R:16:5'): fit.gbtrees works for data from reverse phase column ──
Error in `begin_iteration:end_iteration`: argument of length 0
Backtrace:
1. └─FastRet:::fit_gbtree(df, verbose = 0) at test-fit_gbtree.R:16:5
2.   └─FastRet:::fit_gbtree_grid(...)
3.     └─xgboost::xgb.train(...)
```

fail in:

- https://www.r-project.org/nosvn/R.check/r-devel-linux-x86_64-fedora-clang/FastRet-00check.html
- https://www.r-project.org/nosvn/R.check/r-devel-linux-x86_64-debian-gcc/FastRet-00check.html
- https://www.r-project.org/nosvn/R.check/r-devel-linux-x86_64-fedora-gcc/FastRet-00check.html
- https://www.r-project.org/nosvn/R.check/r-release-linux-x86_64/FastRet-00check.html
- https://www.r-project.org/nosvn/R.check/r-oldrel-windows-x86_64/FastRet-00check.html

I strongly assume, it has to do with the update of xgboost from version 1.7.11
to 3.1.2.1 that was published on 2025-12-03. I.e., we need to check if we can
find a way that works with both xgboost versions (1.x.x) and (3.x.x). If we
cannot, we should introduce a version check and call different code paths
depending on the installed xgboost version. In this case, we should add tests in
.github/workflows for both xgboost versions to make sure that both code paths are
tested in the future. E.g. doing something linke this:

- {os: macos-latest,    r: 'release',  xgboost: 'lastest'}
- {os: windows-latest,  r: 'release',  xgboost: 'lastest'}
- {os: ubuntu-latest,   r: 'devel',    xgboost: 'lastest', http-user-agent: 'release'}
- {os: ubuntu-latest,   r: 'release',  xgboost: 'lastest'}
- {os: ubuntu-latest,   r: 'release',  xgboost: '1.7.11'}
- {os: ubuntu-latest,   r: 'oldrel-1', xgboost: 'lastest'}
- {os: ubuntu-latest,   r: '4.1',      xgboost: 'lastest'}

Instead of

- {os: macos-latest,    r: 'release'}
- {os: windows-latest,  r: 'release'}
- {os: ubuntu-latest,   r: 'devel', http-user-agent: 'release'}
- {os: ubuntu-latest,   r: 'release'}
- {os: ubuntu-latest,   r: 'release'}
- {os: ubuntu-latest,   r: 'oldrel-1'}
- {os: ubuntu-latest,   r: '4.1'}

in `.github/workflows/r-cmd-check.yaml`.

### Allow Adjustment based on chemical descriptors

In `adjust_frm` add argument `add_cds`. If `add_cds = TRUE` chemical descriptors
should be calculated using `preprocess_data()` (as done in `train_frm()`) but
always without interaction terms and without polynomial terms. The resulting CDs
should be included in the adjustment model.

### Add gbtree and glmnet as adjustment models

Currently, adjustment models are only fitted using `lm`. In addition, `glmnet`
and `gbtree` models (as returned by `fit_glmnet()` and `fit_gbtree()`) should be
supported as well. Introduce a new argument `adj_type` to `adjust_frm()` that
can take values `"lm"` (default), `"glmnet"`, and `"gbtree"` and based on this
value decide which adjustment model to fit.
