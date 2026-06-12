# FastRet — R package

Retention-time (RT) prediction in liquid chromatography via QSRR
(Quantitative Structure–Retention Relationship). This is the scientific
software being published; it is on CRAN and has a Shiny web GUI
(<https://fastret.spang-lab.de>). Docs:
<https://spang-lab.github.io/FastRet/>.

This repo is one of three in the workspace — see
[../CLAUDE.md](https://spang-lab.github.io/CLAUDE.md) for how it relates
to **freda** (which depends on it) and **fastret-overleaf** (the
manuscript).

- **Current branch: `1.3.0-dev`** (newest; ahead of `main`). Version
  `1.3.5` (DESCRIPTION).
- License GPL-3. Maintainer: Tobias Schmidt.

## What it does (user-facing capabilities)

1.  **Train** a custom RT model for a chromatography column from
    `(NAME, SMILES, RT)` data (`train_frm`). Methods: `"lasso"`,
    `"ridge"` (glmnet), `"gbtreeDefault"`/`"gbtreeRP"` (xgboost).
2.  **Predict** RTs for new molecules with a trained model
    (`predict.frm`).
3.  **Adjust** an existing model to altered conditions (temperature, pH,
    column age, …) by overlaying an adjustment model (`adjust_frm`).
4.  **Selective measuring** — pick which molecules to re-measure on a
    new column to train an adjustment model cheaply, via ridge-weighted
    PAM clustering (`selective_measuring`).
5.  **Web GUI** (Shiny) and a **CLI** (R console) for all of the above.

## The core object: `frm` (FastRet model)

S3 object of class `frm`. Key fields: `$model` (glmnet or xgb.Booster),
`$df` (training data: `NAME, SMILES, RT` + chemical descriptors), `$cv`
(cross-validation folds/models/stats/preds), `$seed`, `$version`,
`$args`, and `$adj` (the adjustment model, same sub-structure, present
only after `adjust_frm`). Save/load with `saveRDS`/`readRDS`.
Polynomial/interaction terms are **not** stored in `$df` — they’re
regenerated on the fly during prediction, so `predict.frm` must reapply
the training-time transforms (it does this automatically).

## Code map (`R/`)

| File | Role |
|----|----|
| `train.R` | Core: `train_frm`, `predict.frm`, `adjust_frm`, fit helpers (glmnet/gbtree/lm), CV, `get_predictors`, `clip_predictions`, xgboost grid search |
| `getcds.R` | Chemical-descriptor (CD) computation via **rcdk** (`getCDs`); disk+RAM caching; `CDFeatures`/`CDNames` constants; `as_canonical` |
| `prepro.R` | `preprocess_data`: column validation, CD addition, polynomial/interaction/RT terms, NA & near-zero-variance removal (done **per-fold** in CV to avoid leakage) |
| `sm.R` | `selective_measuring`: ridge weighting + PAM clustering; `rt_coef` tuning (`0`, `"max"`, `"inf"`) |
| `plot.R` | `plot_frm`: measured-vs-predicted scatter (CV/train/adjust variants, optional log2 axes) |
| `data.R` | `read_rp_xlsx`, `read_rpadj_xlsx`, `read_retip_hilic_data`, `read_rp_lasso_model_rds`; RP dataset docs |
| `app.R`, `ui.R`, `server.R` | Shiny GUI (see below) |
| `util.R` | Logging (`catf`), context managers (`withTimeout`/`withSink`/…), `collect`, `pkg_file`, `named` |
| `benchmark.R` | `benchmark_find_params`: runtime measurement of the param grid search |

Exports are roxygen-generated; grouped in
[\_pkgdown.yml](https://spang-lab.github.io/FastRet/_pkgdown.yml) by
`@keywords` `public` (user API) / `internal` (exported for parallel
workers, not user-facing) / `dataset`.

## Web GUI (Shiny)

Launch with
[`FastRet::start_gui()`](https://spang-lab.github.io/FastRet/reference/start_gui.md)
(interactive R) → serves on `http://localhost:8080`.
[`fastret_app()`](https://spang-lab.github.io/FastRet/reference/fastret_app.md)
returns the app object without running. Uses **Shiny ≥ 1.8.1
ExtendedTask** to run long jobs (training, selective measuring) on
background
[`future::multisession`](https://future.futureverse.org/reference/multisession.html)
workers (`nw` workers, each up to `nsw` sub-workers) so the UI stays
responsive. Tabs: Train / Predict / Adjust / Selective Measuring. Each
session is isolated. `app.R` is **excluded from the built package**
(`.Rbuildignore`) — it’s a dev launcher.

## Dev workflow

``` r

devtools::load_all()     # develop
devtools::document()     # regenerate NAMESPACE + man/*.Rd from roxygen (never hand-edit these)
devtools::test()         # testthat (edition 3, parallel)
devtools::check()        # R CMD check
lintr::lint_package()    # optional; .lintr is permissive (not enforced in CI)
```

Tests live in `tests/testthat/`. Some are named to run first
(descriptor/cache warm-up, e.g. `getCDs`, `train_frm-*`). `misc/` holds
dev helpers (`scripts/check-version.R` for the CI version check,
`scripts/patch-shiny.R` for dev-mode GUI reload) and is excluded from
the build. This package has **no Dockerfiles** — the deployment
container lives in the workspace’s `docker-images/` repo (see
[../CLAUDE.md](https://spang-lab.github.io/CLAUDE.md)).

**CI** (`.github/workflows/`): `R-CMD-check` (macOS/Windows/Ubuntu ×
several R + xgboost versions), `test-coverage` (covr → codecov),
`pkgdown` (→ GitHub Pages, on `main`), plus rhub and a version-increment
check. The pkgdown site builds from `vignettes/` (GUI-Usage, CLI-Usage,
Package-Internals, Installation, Contributing).

## Git / commits

- **NEVER add a `Co-Authored-By: Claude ...` trailer (or any
  Claude/Anthropic co-author) to commit messages or PR bodies.** Do not
  add yourself as co-author.

## Dependencies & system requirements

- **Java SDK is required** for `rcdk` (descriptor computation). Without
  it, descriptor and canonicalization code fails. Install a JDK via your
  OS package manager and run `R CMD javareconf` if rcdk fails to load
  (see the Installation vignette).
- Key Imports: `rcdk`, `glmnet`, `xgboost`, `data.table`, `future`,
  `cluster`, `shiny` (≥1.8.1), `shinyjs`/`shinyhelper`/`bslib`/`DT`
  (GUI), `openxlsx`, `withr`.
- R ≥ 4.1.0.

## Data

`data/RP.rda` (lazy-loaded RP dataset, ~442 metabolites) and
`inst/extdata/` (RP.xlsx, RP_adj.xlsx, a pre-trained
`RP_lasso_model.rds`, larger measurement workbooks).
`inst/cachedata/CDs.rds` is a precomputed descriptor cache for ~1000
SMILES, loaded into RAM on first `getCDs` call to avoid recomputation.
HILIC data is pulled from the Retip package (CC BY 4.0) via
`read_retip_hilic_data`.

## Gotchas

- **rcdk/Java** is the most common setup failure. Check
  `Sys.which("java")`.
- **Descriptor computation is the slow step** (rcdk, ~0.5–2 s per unique
  SMILES); xgboost training + CV is also slow. Lasso/ridge are fast.
  Caching (disk `CDs.rds` + RAM option) matters — don’t defeat it.
- **CV avoids data leakage** by removing near-zero-var/NA features and
  generating polynomial/interaction terms *inside each fold* (since
  v1.3.0). Preserve this when touching `prepro.R`/`train.R`.
- **Predictions are clipped** by default (`clip=TRUE`) to a log-normal
  tail range of observed RTs; missing descriptors are **mean-imputed**
  from training data (`impute=TRUE`). Both are user-visible behaviors
  reviewers asked about (see the manuscript’s `mean-imputation-bias`
  revision item).
- Reproducibility relies on `frm$seed`; same input + seed ⇒ same model.
