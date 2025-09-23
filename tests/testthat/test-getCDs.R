library(testthat)

test_that("getCDs works correctly", {

    # Helpers
    measure <- function(expr) { system.time(expr)[[3]] }

    # Clear cached CDs
    options(FastRet.cachedCDs = NULL)
    expect_true(is.null(nrow(getOption("FastRet.cachedCDs"))))

    # Simple test case with known molecules from package data
    X <- RP[1:5, ]
    rt1 <- measure(Y1 <- getCDs(X, verbose = 0))
    nc <- ncol(Y1)
    expect_true(all.equal(Y1[, 1:3], X))
    expect_true(all(colnames(Y1)[4:nc] %in% CDFeatures))
    expect_equal(Y1$Fsp3, c(0.75, 1/3, 0.00, 0.50, 0.50))
    expect_true(rt1 < 0.05)
    expect_equal(nrow(getOption("FastRet.cachedCDs")), 1788)

    # Test with unknown SMILES (not in package data)
    smi <- c("CCC", "CCCC", "CCCCC", "CCC")
    df <- data.frame(NAME = paste0("test", 1:4), SMILES = smi, RT = 1:4)
    rt2 <- measure(Y2 <- getCDs(df, verbose = 0))
    expect_equal(nrow(Y2), 4)
    expect_equal(Y2$SMILES, Y2$SMILES)
    expect_equal(ncol(Y2), ncol(df) + length(CDFeatures))
    expect_true(rt2 > rt1)
    expect_equal(nrow(getOption("FastRet.cachedCDs")), 1791)

    # Test with multiple workers
    options(FastRet.cachedCDs = getOption("FastRet.cachedCDs")[1:1788, ]) # reset cache to original
    rt3 <- measure(Y3 <- getCDs(df, verbose = 0, nw = 2))
    expect_equal(Y3, Y2)
    expect_true(rt3 > rt2) # should take more time due to parallel overhead
})
