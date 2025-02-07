library(testthat)

test_that("train_frm works if `method == \"GBTree\"`", {
    set.seed(1)
    model <- train_frm(
        df = RP[1:44, ], # Use only 10% of the data to speed up execution time
        method = "gbtree",
        nfolds = 2,
        nw = 2,
        verbose = 0
    )
    expect_equal(names(model), c("model", "df", "cv", "seed", "version"))
})
