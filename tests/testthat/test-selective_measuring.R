library(testthat)

test_that("selective_measuring works", {
    obj <- selective_measuring(raw_data = RP[1:10, ], k_cluster = 5, verbose = 0)
    expect_true(all(names(obj) == c("clustering", "clobj", "coefs", "model", "df", "dfz", "dfzb")))
})
