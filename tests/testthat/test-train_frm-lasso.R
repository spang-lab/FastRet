library(testthat)

test_that("train_frm works if method == 'lasso'", {
    model1 <- train_frm(df=RP[1:20, ], method="lasso", nfolds=2, nw=1, verbose=0, seed=42)
    model2 <- train_frm(df=RP[1:20, ], method="lasso", nfolds=2, nw=2, verbose=0, seed=42)
    expect_equal(names(model1), c("model", "df", "cv", "seed", "version", "args"))
    model2$args$nw <- 1 # Set to same value as in model1 before comparison

    expect_equal(model1, model2)
})
