# Done

## Remove caret from FastRet

Currently, there are only three functions used from caret:
`caret::createFolds()` and `caret::nearZeroVar()`. These functions should
be replaced with custom implementations to reduce package dependencies.

Done 2025-11-07, 11am, with version 1.2.3.

# Open

## Do not preprocess before CV

Currently, we preprocess the data, then train our model and then estimate its
performance in CV. This is incorrect, as preprocessing should be part of the
training, i.e., we should implement two seperate functions `train_frm()` and
`cv_train_frm()`, with the latter calling the former in each CV fold and
`train_frm()` doing the preprocessing internally.

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

