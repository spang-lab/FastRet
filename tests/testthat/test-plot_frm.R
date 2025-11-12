library(testthat)

test_that("plot_frm works for not-adjusted models", {
    df <- RP[1:20, ]
    frm <- train_frm(df, "lasso", verbose = 0, nfolds = 2, nw = 1, seed = 123)
    testthat::expect_no_error(object = {
        plot_frm(frm, type = "scatter.cv")
        plot_frm(frm, type = "scatter.train")
    })
    expect_error(
        object = plot_frm(frm = frm, type = "scatter.cv.adj"),
        regexp = "the model has not been adjusted yet"
    )
    expect_error(
        object = plot_frm(frm = frm, type = "scatter.train.adj"),
        regexp = "the model has not been adjusted yet"
    )
})

test_that("plot_frm works for adjusted models", {
    frm <- readRDS(pkg_file("extdata/RP_lasso_model.rds"))
    frmadj <- adjust_frm(frm, new_data = read_rpadj_xlsx(), predictors = 1:2, verbose = 0)
    testthat::expect_no_error(object = {
        plot_frm(frm = frmadj, type = "scatter.cv")
        plot_frm(frm = frmadj, type = "scatter.train")
        plot_frm(frm = frmadj, type = "scatter.cv.adj")
        plot_frm(frm = frmadj, type = "scatter.train.adj")
    })
})