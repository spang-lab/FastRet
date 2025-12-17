# Predict retention times using a FastRet Model

Predict retention times for new data using a FastRet Model (FRM).

## Usage

``` r
# S3 method for class 'frm'
predict(
  object = train_frm(),
  df = object$df,
  adjust = NULL,
  verbose = 0,
  clip = TRUE,
  impute = TRUE,
  ...
)
```

## Arguments

- object:

  An object of class `frm` as returned by
  [`train_frm()`](https://spang-lab.github.io/FastRet/reference/train_frm.md).

- df:

  A data.frame with the same columns as the training data.

- adjust:

  If `object` was adjusted using
  [`adjust_frm()`](https://spang-lab.github.io/FastRet/reference/adjust_frm.md),
  it will contain a property `object$adj`. If `adjust` is TRUE,
  `object$adj` will be used to adjust predictions obtained from
  `object$model`. If FALSE `object$adj` will be ignored. If NULL,
  `object$model` will be used, if available.

- verbose:

  A logical value indicating whether to print progress messages.

- clip:

  Clip predictions to be within RT range of training data?

- impute:

  Impute missing predictor values using column means of training data?

- ...:

  Not used. Required to match the generic signature of
  [`predict()`](https://rdrr.io/r/stats/predict.html).

## Value

A numeric vector with the predicted retention times.

## See also

[`train_frm()`](https://spang-lab.github.io/FastRet/reference/train_frm.md),
[`adjust_frm()`](https://spang-lab.github.io/FastRet/reference/adjust_frm.md)

## Examples

``` r
object <- read_rp_lasso_model_rds()
df <- head(RP)
yhat <- predict(object, df)
```
