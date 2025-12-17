library(testthat)

test_that("validate_inputdata errors on missing required columns", {
    df <- data.frame(RT = 1:3, NAME = c("a","b","c"), stringsAsFactors = FALSE)
    expect_error(validate_inputdata(df), "missing columns")
})

test_that("validate_inputdata passes with minimal valid data", {
    frm <- read_rp_lasso_model_rds()
    df <- frm$df[1:2, c("RT", "SMILES", "NAME", CDFeatures[1])]
    expect_invisible(validate_inputdata(df, min_cds = 1))
})

test_that("validate_inputdata errors on unknown columns and min_cds", {
    frm <- read_rp_lasso_model_rds()
    df <- frm$df[1:2, c("RT", "SMILES", "NAME", CDFeatures[1])]
    df$FOO <- 1
    expect_error(validate_inputdata(df), "Unknown columns")
    df2 <- frm$df[1:2, c("RT", "SMILES", "NAME", CDFeatures[1])]
    expect_error(validate_inputdata(df2, min_cds = 5), "At least")
})

test_that("validate_inputmodel detects missing elements", {
    good <- list(model = 1, df = data.frame(), cv = list())
    expect_invisible(validate_inputmodel(good))
    bad <- list(model = 1)
    expect_error(validate_inputmodel(bad), "Model object is")
})
