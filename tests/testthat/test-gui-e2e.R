# End-to-end GUI tests driving the real Shiny app via headless Chrome.
#
# These exercise the actual button-press workflows for all four modes and every
# model type: Train (Lasso + XGBoost), Selective Measuring, Predict, and Adjust
# (Lasso + LM + XGBoost). For each, we set the inputs, click the task button,
# wait for the asynchronous ExtendedTask to resolve, and assert that the
# corresponding output was populated -- i.e. true button-press end-to-end, with
# no UI snapshots.
#
# Requirements: shinytest2 + chromote + a Chrome/Chromium binary + Java (rcdk).
# None are guaranteed on CRAN, hence the skips below.

skip_on_cran()
skip_if_not_installed("shinytest2")
skip_if_not_installed("chromote")
skip_if(is.null(chromote::find_chrome()), "No Chrome/Chromium binary found")
skip_if_not(nzchar(Sys.which("java")), "Java not available (required by rcdk)")

# Build a small xlsx fixture from a subset of the bundled RP dataset so that
# XGBoost training + cross-validation stays fast (all descriptors are cached).
make_subset_xlsx <- function(n = 20, require = c("RT", "NAME", "SMILES")) {
    src <- system.file("extdata", "RP.xlsx", package = "FastRet")
    df <- openxlsx::read.xlsx(src, sheet = 1)
    df <- df[seq_len(min(n, nrow(df))), require, drop = FALSE]
    path <- tempfile(fileext = ".xlsx")
    openxlsx::write.xlsx(df, path, rowNames = FALSE)
    path
}

model_rds <- system.file("extdata", "RP_lasso_model.rds", package = "FastRet")
adj_xlsx  <- system.file("extdata", "RP_adj.xlsx", package = "FastRet")

# A single app instance is reused across the modes to avoid paying the
# Chrome-startup cost four times.
app <- shinytest2::AppDriver$new(
    app_dir = test_path("apps", "fastret"),
    name = "fastret-gui",
    timeout = 30 * 1000,
    load_timeout = 60 * 1000
)
withr::defer(app$stop())

# Generous timeout: XGBoost training + CV on the subset is the slowest path.
TASK_TIMEOUT <- 180 * 1000

test_that("Train mode trains a Lasso model end-to-end", {
    app$set_inputs(navbar = "Train new Model")
    app$upload_file(ubTrainXlsx = make_subset_xlsx())
    app$set_inputs(rbMethod = "1") # 1 = Lasso
    app$click("btnTrain")
    val <- app$wait_for_value(output = "tblTrainResults", timeout = TASK_TIMEOUT)
    expect_false(is.null(val))
})

test_that("Train mode trains an XGBoost model end-to-end", {
    app$set_inputs(navbar = "Train new Model")
    app$upload_file(ubTrainXlsx = make_subset_xlsx())
    app$set_inputs(rbMethod = "2") # 2 = XGBoost
    app$click("btnTrain")
    val <- app$wait_for_value(output = "tblTrainResults", timeout = TASK_TIMEOUT)
    expect_false(is.null(val))
})

test_that("Selective Measuring computes medoids end-to-end", {
    app$set_inputs(navbar = "Selective Measuring")
    app$upload_file(ubSmXlsx = make_subset_xlsx())
    app$set_inputs(niK = 5)
    app$click("btnSM")
    val <- app$wait_for_value(output = "tblMedoids", timeout = TASK_TIMEOUT)
    expect_false(is.null(val))
})

test_that("Predict mode predicts retention times end-to-end", {
    app$set_inputs(navbar = "Predict Retention Times")
    app$upload_file(ubPredFRM = model_rds)
    app$upload_file(ubPredXlsx = adj_xlsx) # has NAME + SMILES columns
    app$click("btnPred")
    val <- app$wait_for_value(output = "tblPredResults", timeout = TASK_TIMEOUT)
    expect_false(is.null(val))
})

test_that("Adjust mode adjusts a model with every method end-to-end", {
    app$set_inputs(navbar = "Adjust existing Model")
    app$upload_file(ubAdjFRM = model_rds)
    app$upload_file(ubAdjXlsx = adj_xlsx)
    for (method in c("lasso", "lm", "gbtree")) {
        app$set_inputs(rbAdjMethod = method)
        app$click("btnAdj")
        val <- app$wait_for_value(output = "poAdjPerf", timeout = TASK_TIMEOUT)
        expect_false(is.null(val))
    }
})
