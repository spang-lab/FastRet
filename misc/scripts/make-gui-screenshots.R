#!/usr/bin/env Rscript
# Regenerate the GUI screenshots used in the documentation.
#
# Drives the real Shiny app with shinytest2 (headless Chrome via chromote),
# runs each of the four modes on a small subset of the bundled RP dataset, and
# captures a populated, whitespace-trimmed PNG of every page. Re-run this after
# any change to the GUI layout.
#
#   Rscript misc/scripts/make-gui-screenshots.R
#
# Requirements: shinytest2 + chromote + a Chrome/Chromium binary + Java (rcdk).
#
# Output: <pkg>/vignettes/GUI-Usage/{train,select,adjust,predict}.png plus a
# 2x2 overview montage gui-overview-2x2.png.

suppressMessages({
    devtools::load_all(quiet = TRUE)
    library(shinytest2)
    library(magick)
})

WIDTH <- 1400L            # capture width (high-res; downscale in the document)
TASK_TIMEOUT <- 180 * 1000
outdir <- "vignettes/GUI-Usage"
appdir <- "tests/testthat/apps/fastret"

subset_xlsx <- function(n = 20) {
    src <- system.file("extdata", "RP.xlsx", package = "FastRet")
    df <- openxlsx::read.xlsx(src, sheet = 1)[seq_len(n), c("RT", "NAME", "SMILES")]
    path <- tempfile(fileext = ".xlsx")
    openxlsx::write.xlsx(df, path)
    path
}
model <- system.file("extdata", "RP_lasso_model.rds", package = "FastRet")
adj   <- system.file("extdata", "RP_adj.xlsx", package = "FastRet")

app <- AppDriver$new(app_dir = appdir, name = "gui-screenshots",
                     width = WIDTH, height = 2200,
                     timeout = 30 * 1000, load_timeout = 60 * 1000)
on.exit(app$stop())

# NB: the console-log output calls invalidateLater(1000), so the app is never
# "idle" -- use wait_for_value() + a fixed paint delay, never wait_for_idle().
shot <- function(name) {
    Sys.sleep(2.5)
    tmp <- tempfile(fileext = ".png")
    app$get_screenshot(tmp, delay = 0.5)
    img <- image_trim(image_read(tmp), fuzz = 2)
    f <- file.path(outdir, paste0(name, ".png"))
    image_write(img, f)
    message(sprintf("wrote %s (%dx%d)", f, image_info(img)$width, image_info(img)$height))
}

app$set_inputs(navbar = "Train new Model")
app$upload_file(ubTrainXlsx = subset_xlsx())
app$set_inputs(rbMethod = "1") # Lasso (fast; still produces CV plot + table)
app$click("btnTrain")
app$wait_for_value(output = "tblTrainResults", timeout = TASK_TIMEOUT)
shot("train")

app$set_inputs(navbar = "Selective Measuring")
app$upload_file(ubSmXlsx = subset_xlsx())
app$set_inputs(niK = 5)
app$click("btnSM")
app$wait_for_value(output = "tblMedoids", timeout = TASK_TIMEOUT)
shot("select")

app$set_inputs(navbar = "Predict Retention Times")
app$upload_file(ubPredFRM = model)
app$upload_file(ubPredXlsx = adj)
app$click("btnPred")
app$wait_for_value(output = "tblPredResults", timeout = TASK_TIMEOUT)
shot("predict")

app$set_inputs(navbar = "Adjust existing Model")
app$upload_file(ubAdjFRM = model)
app$upload_file(ubAdjXlsx = adj)
app$set_inputs(rbAdjMethod = "lasso")
app$click("btnAdj")
app$wait_for_value(output = "poAdjPerf", timeout = TASK_TIMEOUT)
shot("adjust")

# 2x2 overview montage (train | select / adjust | predict)
files <- file.path(outdir, paste0(c("train", "select", "adjust", "predict"), ".png"))
H <- max(vapply(files, function(f) image_info(image_read(f))$height, integer(1)))
pad <- lapply(files, function(f) {
    image_extent(image_read(f), geometry = sprintf("%dx%d", WIDTH, H),
                 color = "white", gravity = "north")
})
row1 <- image_append(c(pad[[1]], pad[[2]]))
row2 <- image_append(c(pad[[3]], pad[[4]]))
grid <- image_border(image_append(c(row1, row2), stack = TRUE), "white", "6x6")
image_write(grid, file.path(outdir, "gui-overview-2x2.png"))
message("wrote ", file.path(outdir, "gui-overview-2x2.png"))
