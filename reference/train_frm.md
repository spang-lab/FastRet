# Train a new FastRet model (FRM) for retention time prediction

Trains a new model from molecule SMILES to predict retention times (RT)
using the specified method.

## Usage

``` r
train_frm(
  df,
  method = "lasso",
  verbose = 1,
  nfolds = 5,
  nw = 1,
  degree_polynomial = 1,
  interaction_terms = FALSE,
  rm_near_zero_var = TRUE,
  rm_na = TRUE,
  rm_ns = FALSE,
  seed = NULL,
  do_cv = TRUE
)
```

## Arguments

- df:

  A dataframe with columns "NAME", "RT", "SMILES" and optionally a set
  of chemical descriptors. If no chemical descriptors are provided, they
  are calculated using the function
  [`preprocess_data()`](https://spang-lab.github.io/FastRet/reference/preprocess_data.md).

- method:

  A string representing the prediction algorithm. Either "lasso",
  "ridge", "gbtree", "gbtreeDefault" or "gbtreeRP". Method "gbtree" is
  an alias for "gbtreeDefault".

- verbose:

  A logical value indicating whether to print progress messages.

- nfolds:

  An integer representing the number of folds for cross validation.

- nw:

  An integer representing the number of workers for parallel processing.

- degree_polynomial:

  An integer representing the degree of the polynomial. Polynomials up
  to the specified degree are included in the model.

- interaction_terms:

  A logical value indicating whether to include interaction terms in the
  model.

- rm_near_zero_var:

  A logical value indicating whether to remove near zero variance
  predictors.

- rm_na:

  A logical value indicating whether to remove NA values before
  training. Highly recommended to avoid issues during model fitting.
  Setting this to FALSE with `method = "lasso"` will most likely lead to
  errors.

- rm_ns:

  A logical value indicating whether to remove chemical descriptors that
  were considered as not suitable for linear regression based on a
  previous analysis of an independent dataset. Currently not used.

- seed:

  An integer value to set the seed for random number generation to allow
  for reproducible results.

- do_cv:

  A logical value indicating whether to perform cross-validation. If
  FALSE, the `cv` element in the returned object will be NULL.

## Value

A 'FastRet Model', i.e., an object of class `frm`. Components are:

- `model`: The fitted base model. This can be an object of class
  `glmnet` (for Lasso or Ridge regression) or `xgb.Booster` (for GBTree
  models).

- `df`: The data frame used for training the model. The data frame
  contains all user-provided columns (including mandatory columns RT,
  SMILES and NAME) as well the calculated chemical descriptors. (But no
  interaction terms or polynomial features, as these can be recreated
  within a few milliseconds).

- `cv`: A named list containing the cross validation results, or NULL if
  `do_cv = FALSE`. When not NULL, elements are:

  - `folds`: A list of integer vectors specifying the samples in each
    fold.

  - `models`: A list of models trained on each fold.

  - `stats`: A list of vectors with RMSE, Rsquared, MAE, pBelow1Min per
    fold.

  - `preds`: Retention time predictions obtained in CV as numeric
    vector.

- `seed`: The seed used for random number generation.

- `version`: The version of the FastRet package used to train the model.

- `args`: The value of function arguments besides `df` as named list.

## Examples

``` r
m <- train_frm(df = RP[1:40, ], method = "lasso", nfolds = 2, verbose = 0)
# For the sake of a short runtime, only the first 40 rows of the RP dataset
# are used in this example. In practice, you should always use the entire
# training dataset for model training.
```
