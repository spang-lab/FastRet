library(testthat)

opts <- options(testthat.show_trace = FALSE)

data_funcs <- c(
    # Full Datasets
    "read_rp_xlsx",
    "read_retip_hilic_data",
    ".read_hilic_xlsx",
    # ".read_rp_axmm_xlsx", # File missing. Fadi is looking for it.
    # Smaller Datasets for Model Adjustment
    ".read_rp_steep_xlsx",
    ".read_rp_flat_xlsx",
    ".read_rp_t25_xlsx",
    ".read_rp_fr25_xlsx",
    ".read_rp_t25_fr25_xlsx",
    ".read_rp_t25_fr25_steep_xlsx",
    # Validation Data Sets
    ".read_rp_steep_val_xlsx",
    ".read_rp_flat_val_xlsx",
    ".read_rp_t25_val_xlsx",
    ".read_rp_fr25_val_xlsx",
    ".read_rp_t25_fr25_val_xlsx",
    ".read_rp_t25_fr25_steep_val_xlsx",
    # For Package Testing
    "read_rpadj_xlsx",
    # "read_rp_lasso_model_rds", # Model, not a dataset.
    NULL
)
# (1) This excel file is missing. Fadi is looking for it.

if (pkg_file("misc/datasets") == "") {
    # 'misc/datasets' is not available, i.e., we are not in a
    # development environment and we can not use/test our private
    # functions depending on raw data files not shipped with the
    # package (e.g. the source files for lazy loaded datasets).
    # Therefore we filter out all functions starting with a dot.
    data_funcs <- data_funcs[!startsWith(data_funcs, ".read_")]
}

test_that("read_xxx_works", {
    # For every function we check wheter:
    # 1. It can be called without throwing an error
    # 2. It returns a data.frame object
    # 3. The data.frame has at least columns NAME and SMILES
    # 4. There should be now row with only NAs
    lapply(data_funcs, function(fn) {
        f <- get(fn)
        df <- try(f(), silent = TRUE)
        errmsg <- sprintf("Test for df <- %s() failed", fn)
        expect_false(inherits(df, "try-error"), errmsg)
        expect_true(is.data.frame(df), errmsg)
        expect_true(all(c("NAME", "SMILES") %in% toupper(colnames(df))), errmsg)
        expect_true(all(apply(df, 1, function(x) !all(is.na(x)))), errmsg)
    })
})
