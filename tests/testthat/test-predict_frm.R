library(testthat)

test_that("predict.frm returns numeric vector (unadjusted)", {
    frm <- read_rp_lasso_model_rds()
    df <- frm$df[1:5, ]
    yhat <- predict(frm, df, adjust = FALSE, verbose = 0)
    expect_type(yhat, "double")
    expect_length(yhat, nrow(df))
})

test_that("predict.frm applies adjustment when available", {
    skip_on_cran()
    frm <- read_rp_lasso_model_rds()
    new_data <- read_rpadj_xlsx()
    frm2 <- adjust_frm(frm, new_data, predictors = 1:2, nfolds = 2, verbose = 0)
    df <- frm$df[1:5, ]
    yhat_adj <- predict(frm2, df, adjust = TRUE, verbose = 0)
    expect_type(yhat_adj, "double")
    expect_length(yhat_adj, nrow(df))
})
