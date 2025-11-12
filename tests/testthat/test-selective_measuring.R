library(testthat)

test_that("selective_measuring works", {
    RP10 <- RP[1:10, ]
    obj <- selective_measuring(raw_data = RP10, k_cluster = 5, verbose = 0)
    onam <- names(obj)
    xnam <- c("clustering", "clobj", "coefs", "model", "df", "dfz", "dfzb")
    expect_identical(onam, xnam)
})
