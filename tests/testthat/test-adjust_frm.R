library(testthat)

test_that("adjust_frm works", {
    frm = readRDS(pkg_file("extdata/RP_lasso_model.rds"))
    new_data = read_rpadj_xlsx()
    m2 <- adjust_frm(frm, new_data, predictors = 1:2)
    m6 <- adjust_frm(frm, new_data, predictors = 1:6)
    expect_identical(
        object = names(coef(m2$adj$model)),
        expected = c("(Intercept)", "RT", "I(RT^2)")
    )
    expect_identical(
        object = names(coef(m6$adj$model)),
        expected = c("(Intercept)", "RT", "I(RT^2)", "I(RT^3)", "log(RT)", "exp(RT)", "sqrt(RT)")
    )
    expect_error(
        object = adjust_frm(frm, new_data, predictors = c()),
        regexp = ".*Invalid predictors.*"
    )
})
