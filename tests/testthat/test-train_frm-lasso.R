library(testthat)

test_that("train_frm works if method == 'lasso'", {
    model1 <- train_frm(df=RP[1:20, ], method="lasso", nfolds=2, nw=1, verbose=0, seed=42)
    model2 <- train_frm(df=RP[1:20, ], method="lasso", nfolds=2, nw=2, verbose=0, seed=42)
    expect_equal(names(model1), c("model", "df", "cv", "seed", "version", "args"))
    model2$args$nw <- 1 # Set to same value as in model1 before comparison

    expect_equal(model1, model2)
})

test_that("train_frm works with do_cv = FALSE", {
    model1 <- train_frm(df=RP[1:20, ], method="lasso", nfolds=2, verbose=0, seed=42, do_cv=FALSE)
    model2 <- train_frm(df=RP[1:20, ], method="lasso", nfolds=2, verbose=0, seed=42, do_cv=TRUE)
    
    # Model with do_cv=FALSE should have NULL cv element
    expect_null(model1$cv)
    expect_equal(names(model1), c("model", "df", "cv", "seed", "version", "args"))
    
    # Model with do_cv=TRUE should have cv element
    expect_false(is.null(model2$cv))
    expect_true(is.list(model2$cv))
    expect_equal(names(model2$cv), c("folds", "models", "stats", "preds"))
    
    # The actual models should be the same (only cv differs)
    expect_equal(model1$model$beta, model2$model$beta)
    expect_equal(model1$df, model2$df)
})
