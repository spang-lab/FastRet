library(testthat)

test_that("train_frm works if `method == \'lasso\'`", {
    set.seed(1)
    model <- train_frm(
        df = RP[1:88, ], # Use only 20% of the data to speed up execution time
        method = "lasso",
        nfolds = 2,
        nw = 2,
        verbose = 0
    )
    expect_equal(names(model), c("model", "df", "cv", "seed", "version"))
})
