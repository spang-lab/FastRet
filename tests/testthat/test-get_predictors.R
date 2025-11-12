library(testthat)

test_that("get_predictors returns predictor names for glmnet model", {
    set.seed(1)
    frm <- train_frm(df = RP[1:40, ], method = "lasso", nfolds = 2, nw = 1, verbose = 0)
    preds <- get_predictors(frm)
    expect_true(is.character(preds))
    expect_true(length(preds) > 10) # should have many descriptors
})

test_that("get_predictors returns predictor names for gbtree model", {
    skip_on_cran()
    set.seed(1)
    frm <- train_frm(df = RP[1:30, ], method = "gbtree", nfolds = 2, nw = 1, verbose = 0)
    preds <- get_predictors(frm)
    expect_true(is.character(preds))
    expect_true(length(preds) > 10)
})
