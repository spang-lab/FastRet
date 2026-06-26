library(testthat)

# Helper: write a multi-sheet workbook to a temp file and return its path.
write_workbook <- function(sheets) {
    path <- tempfile(fileext = ".xlsx")
    openxlsx::write.xlsx(sheets, path)
    path
}

test_that("drop_unsupported_columns keeps only known columns", {
    df <- data.frame(
        NAME = "a", RT = 1, SMILES = "CCO", INCHIKEY = "x", RT_ADJ = 2,
        Comment = "note", BatchID = 3, mz = 1.5,
        stringsAsFactors = FALSE
    )
    df[[CDFeatures[1]]] <- 0.1
    out <- drop_unsupported_columns(df)
    expect_true(all(c("NAME", "RT", "SMILES", "INCHIKEY", "RT_ADJ", CDFeatures[1]) %in% colnames(out)))
    expect_false(any(c("Comment", "BatchID", "mz") %in% colnames(out)))
})

test_that("read_input_xlsx drops extra columns from the first sheet", {
    df <- data.frame(
        RT = c(1, 2, 3), NAME = c("a", "b", "c"), SMILES = c("CCO", "CCC", "CCCC"),
        Comment = "note", BatchID = 1:3, mz = rnorm(3),
        stringsAsFactors = FALSE
    )
    path <- write_workbook(list(Data = df))
    out <- read_input_xlsx(path, require = c("RT", "SMILES", "NAME"))
    expect_setequal(colnames(out), c("RT", "NAME", "SMILES"))
    expect_equal(nrow(out), 3)
})

test_that("read_input_xlsx finds the required columns on a non-first sheet", {
    cover <- data.frame(Title = "My experiment", Author = "FF", stringsAsFactors = FALSE)
    data <- data.frame(
        RT = c(1, 2), NAME = c("a", "b"), SMILES = c("CCO", "CCC"),
        Note = c("x", "y"), stringsAsFactors = FALSE
    )
    path <- write_workbook(list(Cover = cover, Measurements = data))
    out <- read_input_xlsx(path, require = c("RT", "SMILES", "NAME"))
    expect_setequal(colnames(out), c("RT", "NAME", "SMILES"))
    expect_equal(out$NAME, c("a", "b"))
})

test_that("read_input_xlsx errors when no sheet has the required columns", {
    s1 <- data.frame(Foo = 1:2, Bar = 3:4)
    s2 <- data.frame(RT = 1:2, NAME = c("a", "b")) # missing SMILES
    path <- write_workbook(list(One = s1, Two = s2))
    expect_error(
        read_input_xlsx(path, require = c("RT", "SMILES", "NAME")),
        "required columns"
    )
})

test_that("read_input_xlsx honours a reduced require set (prediction input)", {
    df <- data.frame(
        NAME = c("a", "b"), SMILES = c("CCO", "CCC"), Comment = "z",
        stringsAsFactors = FALSE
    )
    path <- write_workbook(list(Data = df))
    out <- read_input_xlsx(path, require = c("NAME", "SMILES"))
    expect_setequal(colnames(out), c("NAME", "SMILES"))
})

test_that("a workbook with extra columns trains end-to-end after sanitisation", {
    # Regression test for the GUI hang when uploading workbooks with extra
    # columns: the upload pipeline must reduce the data to the supported columns
    # so that train_frm() receives clean input.
    df <- FastRet::RP[1:30, c("RT", "NAME", "SMILES")]
    df$Comment <- "note"
    df$BatchID <- rep(1:3, length.out = nrow(df))
    df$mz <- rnorm(nrow(df))
    path <- write_workbook(list(Sheet1 = df))
    clean <- read_input_xlsx(path, require = c("RT", "SMILES", "NAME"))
    expect_false(any(c("Comment", "BatchID", "mz") %in% colnames(clean)))
    m <- train_frm(df = clean, method = "lasso", nfolds = 3, nw = 1, verbose = 0)
    expect_s3_class(m, "frm")
    expect_false(is.null(m$cv))
})

test_that("fastret_app raises the Shiny upload size limit to >= 500 MB", {
    withr::with_options(list(shiny.maxRequestSize = NULL), {
        app <- fastret_app()
        expect_gte(getOption("shiny.maxRequestSize"), 500 * 1024^2)
    })
})
