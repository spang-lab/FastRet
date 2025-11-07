test_that("createFolds works", {
    y <- rnorm(10, 5, 1) # e.g. 10 retention time scattered around 5 minutes
    folds <- createFolds(y, 3)
    expect_equal(unname(sort(lengths(folds))), c(3, 3, 4))
    expect_equal(names(folds), c("Fold1", "Fold2", "Fold3"))
})
