# Done

## Remove caret from FastRet

Currently, there are only three functions used from caret:
`caret::createFolds()` and `caret::nearZeroVar()`. These functions should
be replaced with custom implementations to reduce package dependencies.

Done 2025-11-07, 11am, with version 1.2.3.

# Open

## Ensure predict.frm handles polynomial/interaction features

Problem
- Models can be trained with polynomial features (degree > 1) and/or interaction terms via `preprocess_data()`, but `predict.frm()` only ensures base chemical descriptors are present (via `getCDs()`). It does not reconstruct polynomial or interaction columns. If a trained model's predictors include names like `Fsp3^2` or `Fsp3*apol`, prediction will fail with "predictors are missing in df".

Impact
- Users training with `degree_polynomial >= 2` or `interaction_terms = TRUE` cannot deploy the model to new data unless they manually compute matching columns.

Proposed fix
- In `predict.frm()` grep the feature names of the `frm` base model (available
  via `get_predictors()`) for `^` and `*` to identify polynomial and interaction
  features. If any are found, verify that these features are indeed composed
  from elements of `CDFeatures`. If yes, call `preprocess_data` with
  `degree_polynomial=...`, `interaction_terms=...`, `rm_near_zero_var=FALSE` and
  `rm_na=FALSE)` to reconstruct polynomial and interaction features without
  re-filtering by variance/NA across the full dataset.
- Then subset to `get_predictors(object)` and predict.
- Nice addition: modify `preprocess_data()` to store function arguments optionally as attributes of the returned data frame, so that `predict.frm()` can read `degree_polynomial` and `interaction_terms` from `attributes(frm$df)` instead of having to grep column names. (But we still want to have the grepping functionality as fallback).

Acceptance criteria
- A model trained with polynomials and/or interactions can be used to predict on new data containing only NAME/SMILES/RT or NAME/SMILES without manual preprocessing.
- Add unit tests covering polynomial+interaction training and subsequent prediction on new data (with and without pre-computed CDs).

## Document frm object structure

Document the structure of `frm` objects created by `train_frm()`, including
all its components and their meanings. This documentation should be added to a
dedicated help file. The corresponding roxygen2 documentation should be found in
`R/train.R`, near the `train_frm()` function definition.

After the docs exist, search all functions for places where `frm` objects are
mentioned and replace the word with an actual link to the `frm` documentation.

## Do not preprocess before CV

Currently, we preprocess the data, then train our model and then estimate its
performance in CV. This is incorrect, as preprocessing should be part of the
training, i.e., we should implement two seperate functions `train_frm()` and
`cv_train_frm()`, with the latter calling the former in each CV fold and
`train_frm()` doing the preprocessing internally.

Plan:

1. Make calculate of chemical descriptors in `preprocess_data()` optional
2. Implement a function `cross_validate`, that takes a fit function plus other
   required args and then performs a cross validation where the fit function
   is evaluated on each fold.
3. Change the `train_frm` as follows:
   - Instead of doing the whole preprocessing in the beginning, only calculate
     chemical descriptors if they are missing. Do not remove zero variance
     features or add interaction terms or other data transformations. I.e.,
     generated `df_cdk` from `df_raw`.
   - Call `fit_func` to fit the model on `df_raw`. All data transformations
     should be done inside `fit_func`.
   - Call `cross_validate()` with `fit_func` to evaluate `fit_func` via CV.
     `cross_validate()` must do nothing else than splitting the data and calling
     `fit_func` (with all preprocessing args) on each fold.


## Make scaling of RT in SM optional

In selective measuring, we scale RTs by the largest ridge coefficent. This
should be made a parameter of the function. Instead of max scaling, we should
allows the user to choose between "max" and "none" (we use strings in case we
want to add more options later on).

## Accept InCHIKEYs in FastRet

Change the implementation of `train_frm`, `adjust_frm`,
`selective_measuring`, etc. to accept an INCHIKEY column in addition to
the NAME, SMILES and RT columns. If that column is present, INCHIKEYs are
used for mapping, instead of NAME/SMILES combinations.

## Test FastRet with rcdk v3.3

Make sure we test once with an old version of rcdk and corresponding rcdlibs.
Install these in the github workflow.
Set the smallest test version as dependency in the DESCRIPTION file.

## Make bounding in predict_frm optional

Make bounding of predicted RTs to the min/max +- x% of training RTs optional
via a function parameter in `predict_frm()`. Default should be 33%.

## Check that all parallel calls work on Windows

Check all places where we use parallel processing and make sure they work
on Windows as well.
