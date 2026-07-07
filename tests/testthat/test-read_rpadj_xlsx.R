library(testthat)

test_that("read_rpadj_xlsx works", {
    x <- read_rpadj_xlsx()
    expect_equal(dim(x), c(25, 4))
    expect_setequal(colnames(x), c("RT", "NAME", "SMILES", "INCHIKEY"))
    # The adjustment demo requires every re-measured metabolite to also be in the
    # RP training set (matched by SMILES and INCHIKEY).
    expect_true(all(x$SMILES %in% RP$SMILES))
    expect_true(all(x$INCHIKEY %in% RP$INCHIKEY))
})
