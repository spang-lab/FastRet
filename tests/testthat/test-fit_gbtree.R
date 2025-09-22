library(testthat)

test_that("fit.gbtrees works as expected", {
    set.seed(123) # for reproducibility
    n <- 100 # number of observations
    Fsp3 <- rnorm(n); nSmallRings <- rnorm(n); nAromRings <- rnorm(n); noise <- rnorm(n); RT <- Fsp3^2 + sin(nSmallRings) + noise
    df <- data.frame(RT, Fsp3, nSmallRings, nAromRings)
    result <- fit_gbtree(df, verbose = 0)
    expect_true(inherits(result, "xgb.Booster"))
})

test_that("fit.gbtrees works for data from reverse phase column", {
    df <- preprocess_data(data = RP[1:5, ], verbose = 0)[1:5, 1:10] # use a smaller dataset to speed up the test
    result <- fit_gbtree(df, verbose = 0)
    expect_true(inherits(result, "xgb.Booster"))
})
