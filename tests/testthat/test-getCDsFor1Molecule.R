library(testthat)

test_cached_cds <- function(smiles = "O=C(O)CCCCCCCCCO") {
    # Check whether the cached descriptors can be reproduced by non-cached
    # calculation
    x <- getCDsFor1Molecule(smiles, cache = FALSE, verbose = 0)
    y <- getCDsFor1Molecule(smiles, cache = TRUE, verbose = 0)
    rownames(y) <- smiles
    testthat::expect_equal(object = x, expected = y)
    testthat::expect_equal(colnames(x), CDFeatures)
    x
}

test_that("getCDs works correctly ", {

    # Test fixed SMILES String
    x <- test_cached_cds(smiles = "O=C(O)CCCCCCCCCO")
    testthat::expect_equal(which(is.na(x)), 29) # only descriptor 29 (`geomShape`) is NA

    # Randomly draw three more SMILES from the RP set and test them
    for (smi in sample(RP$SMILES, 3)) {
        test_cached_cds(smi)
    }
})
