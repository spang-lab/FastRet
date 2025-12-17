library(testthat)

test_that("train_frm works if `method == \"GBTree\"`", {
    model1 <- train_frm(df=RP[1:20, ], method="gbtree", nfolds=2, nw=1, verbose=0, seed=42)
    model2 <- train_frm(df=RP[1:20, ], method="gbtree", nfolds=2, nw=2, verbose=0, seed=42)
    expect_equal(names(model1), c("model", "df", "cv", "seed", "version", "args"))
    model2$args$nw <- 1 # Set to same value as in model1 before comparison
    expect_true(all.equal(model1, model2)) # (1)
    # (1) we can't use expect_equal because the xgboost objects contain pointers
    # to the underlying C structs which causes expect_equal to fail
})
