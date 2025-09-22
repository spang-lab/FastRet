library(testthat)

test_that("getCDs works correctly", {

    # Helpers
    rt <- function(expr) { system.time(expr)[[3]] }

    # Simple as possible case. One core, all SMILES cached.
    X <- RP[1:5, ] # Only check first 5 molecules for speedup
    rt1 <- rt(Y1 <- getCDs(X, verbose = 0))
    nc <- ncol(Y1)
    expect_true(all.equal(Y1[, 1:3], X))
    expect_true(all(colnames(Y1)[4:nc] %in% CDFeatures))
    expect_equal(Y1$Fsp3, c(0.75, 1/3, 0.00, 0.50, 0.50))
    expect_true(rt1 < 0.05) # should take less than 0.05s

    # Now with some uncached SMILES
    rc_clear(X$SMILES)
    dc_clear(X$SMILES)
    rt2 <- rt(Y2 <- getCDs(X, verbose = 0))
    expect_equal(Y2, Y1)
    expect_true(rt2 > 10 * rt1)

    # The following should NOT start a cluster, because force_nw is FALSE and
    # less than 20 SMILES are uncached.
    rt3 <- rt(Y3 <- getCDs(X, verbose = 0, nw = 2))
    expect_equal(Y3, Y1)
    expect_true(rt3 < 0.05)

    # Now we remove some SMILES from the cache again and set force_nw to TRUE.
    # I.e. we will actually use a cluster.
    rc_clear(X$SMILES[1:3])
    dc_clear(X$SMILES[1:3])
    rt4 <- rt(Y4 <- getCDs(X, verbose = 0, nw = 2, force_nw = TRUE))
    expect_equal(Y4, Y1)
    expect_true(rt4 > 10 * rt1)
})
