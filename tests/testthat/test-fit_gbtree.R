library(testthat)

test_that("fit.gbtrees works as expected", {
    set.seed(123) # for reproducibility
    n <- 100 # number of observations
    Fsp3 <- rnorm(n); nSmallRings <- rnorm(n); nAromRings <- rnorm(n); noise <- rnorm(n); RT <- Fsp3^2 + sin(nSmallRings) + noise
    df <- data.frame(RT, Fsp3, nSmallRings, nAromRings)
    X <- df[, c("Fsp3", "nSmallRings", "nAromRings")]
    y <- df$RT
    result <- fit_gbtree(X, y, verbose = 0)
    expect_true(inherits(result, "xgb.Booster"))
})

test_that("fit.gbtrees works for data from reverse phase column", {
    df_full <- preprocess_data(data = RP[1:20, ], verbose = 0, rm_near_zero_var = FALSE)
    meta <- which(colnames(df_full) %in% c("NAME", "SMILES", "RT", "INCHIKEY"))
    X <- df_full[1:20, - meta]
    y <- df_full[1:20, "RT"]
    result <- fit_gbtree(X, y, verbose = 0)
    expect_true(inherits(result, "xgb.Booster"))
})
