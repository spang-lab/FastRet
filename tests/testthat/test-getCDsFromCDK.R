library(testthat)

test_that("getCDsFromCDK works", {
    x <- sample(RP$SMILES, 2)
    y1 <- getCDsFromCDK(x[1]) # 0.16s
    y2 <- getCDsFromCDK(x) # 0.81 s
    y3 <- getCDsFromCDK(x, nw = 2) # 1.98 s (2)

    expect_equal(class(y1), "data.frame")
    expect_equal(class(y2), "data.frame")

    expect_equal(nrow(y1), 1)
    expect_equal(nrow(y2), 2)

    expect_equal(ncol(y1), length(CDFeatures))
    expect_equal(ncol(y2), length(CDFeatures))

    expect_equal(colnames(y1), CDFeatures)
    expect_equal(colnames(y2), CDFeatures)

    expect_equal(y1, y2[1, ])
    expect_equal(y2, y3)
})
