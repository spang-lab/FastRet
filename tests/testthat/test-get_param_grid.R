library(testthat)

test_that("get_param_grid tiny has expected shape", {
    g <- get_param_grid("tiny", nthread = 2)
    expect_true(all(colnames(g) %in% c("eta", "nthread")))
    expect_equal(nrow(g), 4)
    expect_true(all(g$nthread == 2))
})
