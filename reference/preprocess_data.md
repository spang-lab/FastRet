# Preprocess data

Preprocess data so they can be used as input for
[`train_frm()`](https://spang-lab.github.io/FastRet/reference/train_frm.md).

## Usage

``` r
preprocess_data(
  data,
  degree_polynomial = 1,
  interaction_terms = FALSE,
  verbose = 1,
  nw = 1,
  rm_near_zero_var = TRUE,
  rm_na = TRUE,
  add_cds = TRUE,
  rm_ucs = TRUE,
  rt_terms = 1,
  mandatory = c("NAME", "RT", "SMILES")
)
```

## Arguments

- data:

  Dataframe with following columns:

  - Mandatory: NAME, RT and SMILES.

  - Recommmended: INCHIKEY.

  - Optional: Any of the chemical descriptors listed in
    [CDFeatures](https://spang-lab.github.io/FastRet/reference/CDFeatures.md).
    All other columns will be removed. See 'Details'.

- degree_polynomial:

  Add predictors with polynomial terms up to the specified degree, e.g.
  2 means "add squares", 3 means "add squares and cubes". Set to 1 to
  leave descriptors unchanged.

- interaction_terms:

  Add interaction terms? Polynomial terms are not included in the
  generation of interaction terms.

- verbose:

  0: no output, 1: show progress, 2: progress and warnings.

- nw:

  Number of workers to use for parallel processing.

- rm_near_zero_var:

  Remove near zero variance predictors?

- rm_na:

  Remove NA values?

- add_cds:

  Add chemical descriptors using
  [`getCDs()`](https://spang-lab.github.io/FastRet/reference/getCDs.md)?
  See 'Details'.

- rm_ucs:

  Remove unsupported columns?

- rt_terms:

  Which retention-time transformations to append as extra predictors.
  Supply a numeric vector referencing predefined rt_terms (1=RT,
  2=I(RT^2), 3=I(RT^3), 4=log(RT), 5=exp(RT), 6=sqrt(RT)) or a character
  vector with the explicit transformation terms. Character values are
  passed to [`model.frame()`](https://rdrr.io/r/stats/model.frame.html),
  so they must use valid formula syntax (e.g. "I(RT^2)" rather than
  "RT^2").

- mandatory:

  Character vector of mandatory columns that must be present in `data`.
  If any of these columns are missing, an error is raised.

## Value

A dataframe with the preprocessed data.

## Details

If `add_cds = TRUE`, chemical descriptors are added using
[`getCDs()`](https://spang-lab.github.io/FastRet/reference/getCDs.md).
If **all** chemical descriptors listed in
[CDFeatures](https://spang-lab.github.io/FastRet/reference/CDFeatures.md)
are already present in the input `data` object,
[`getCDs()`](https://spang-lab.github.io/FastRet/reference/getCDs.md)
will leave them unchanged. If one or more chemical descriptors are
missing, **all** chemical descriptors will be recalculated and existing
ones will be overwritten.

## Examples

``` r
data <- head(RP, 3)
pre <- preprocess_data(data, verbose = 0)
```
