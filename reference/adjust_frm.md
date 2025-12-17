# Adjust an existing FastRet model for use with a new column

The goal of this function is to train a model that predicts RT_ADJ
(retention time measured on a new, adjusted column) from RT (retention
time measured on the original column) and to attach this adjustment
model to an existing FastRet model.

## Usage

``` r
adjust_frm(
  frm,
  new_data,
  predictors = 1:6,
  nfolds = 5,
  verbose = 1,
  seed = NULL,
  do_cv = TRUE,
  adj_type = "lm",
  add_cds = NULL
)
```

## Arguments

- frm:

  An object of class `frm` as returned by
  [`train_frm()`](https://spang-lab.github.io/FastRet/reference/train_frm.md).

- new_data:

  Data frame with required columns "RT", "NAME", "SMILES"; optional
  "INCHIKEY". "RT" must be the retention time measured on the adjusted
  column. Each row must match at least one row in `frm$df`. The exact
  matching behavior is described in 'Details'.

- predictors:

  Numeric vector specifying which transformations to include in the
  model. Available options are: 1=RT, 2=RT^2, 3=RT^3, 4=log(RT),
  5=exp(RT), 6=sqrt(RT). Note that predictor 1 (RT) is always included,
  even if not specified explicitly.

- nfolds:

  The number of folds for cross validation.

- verbose:

  Show progress messages?

- seed:

  An integer value to set the seed for random number generation to allow
  for reproducible results.

- do_cv:

  A logical value indicating whether to perform cross-validation. If
  FALSE, the `cv` element in the returned adjustment object will be
  NULL.

- adj_type:

  A string representing the adjustment model type. Either "lm", "lasso",
  "ridge", or "gbtree".

- add_cds:

  A logical value indicating whether to add chemical descriptors as
  predictors to new data. Default is TRUE if `adj_type` is "lasso",
  "ridge" or "gbtree" and FALSE if `adj_type` is "lm".

## Value

An object of class `frm`, as returned by
[`train_frm()`](https://spang-lab.github.io/FastRet/reference/train_frm.md),
but with an additional element `adj` containing the adjustment model.
Components of `adj` are:

- `model`: The fitted adjustment model. Class depends on `adj_type` and
  is one of `lm`, `glmnet`, or `xgb.Booster`.

- `df`: The data frame used for training the adjustment model. Including
  columns "NAME", "SMILES", "RT", "RT_ADJ" and optionally "INCHIKEY", as
  well as any additional predictors specified via the `predictors`
  argument.

- `cv`: A named list containing the cross validation results (see
  'Details'), or NULL if `do_cv = FALSE`. When not NULL, elements are:

  - `folds`: A list of integer vectors specifying the samples in each
    fold.

  - `models`: A list of adjustment models trained on each fold.

  - `stats`: A list of vectors with RMSE, Rsquared, MAE, pBelow1Min per
    fold. Added with v1.3.0.

  - `preds`: Retention time predictions obtained during CV by applying
    the adjustment model to the hold-out data.

  - `preds_adjonly`: Removed (i.e. NULL) since v1.3.0.

- `args`: Function arguments used for adjustment (excluding `frm`,
  `new_data` and `verbose`). Added with v1.3.0.

- `version`: The version of the FastRet package used to train the
  adjustment model. Added with v1.3.0.

## Details

Matching is done via "SMILES"+"INCHIKEY" if both datasets have
non-missing INCHIKEYs for all rows; otherwise via "SMILES"+"NAME". If
multiple rows in `frm$df` match the same row in `new_data`, their RT
values are averaged first, and this average is used for training the
adjustment model.

Example: if `frm$df` equals data.frame OLD shown below and `new_data`
equals data.frame NEW, then the resulting, paired data.frame will look
like PAIRED.

    OLD <- data.frame(
        NAME   = c("A", "B",  "B",  "C"  ),
        SMILES = c("C", "CC", "CC", "CCC"),
        RT     = c(5.0,  8.0,  8.2,  9.0 )
    )
    NEW <- data.frame(
        NAME   = c("A", "B",  "B",  "B"),
        SMILES = c("C", "CC", "CC", "CC"),
        RT     = c(2.5,  5.5,  5.7,  5.6)
    )
    PAIRED <- data.frame(
        NAME   = c("A", "B",  "B",  "B"),
        SMILES = c("C", "CC", "CC", "CC"),
        RT     = c(5.0,  8.1,  8.1,  8.1), # Average of OLD$RT[2:3]
        RT_ADJ = c(2.5,  5.5,  5.7,  5.6)  # Taken from NEW
    )

If `do_cv` is TRUE, the adjustment procedure is evaluated in
cross-validation. However, care must be taken when interpreting the CV
results, as the model performance depends on both the adjustment layer
and the original model, which was trained on the full base dataset.
Therefore, the observed CV metrics should be read as "expected
performance when predicting RTs for molecules that were part of the
base-model training but not part of the adjustment set" instead of
"expected performance when predicting RTs for completely new molecules".

## Examples

``` r
frm <- read_rp_lasso_model_rds()
new_data <- read_rpadj_xlsx()
frm_adj <- adjust_frm(frm, new_data, verbose = 0)
```
