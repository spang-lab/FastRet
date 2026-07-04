library(testthat)

test_that("get_param_grid tiny has expected shape", {
    g <- get_param_grid("tiny", nthread = 2)
    expect_true(all(colnames(g) %in% c("eta", "nthread")))
    expect_equal(nrow(g), 4)
    expect_true(all(g$nthread == 2))
})

test_that("get_param_grid small is the 8-combination Retip-style grid", {
    g <- get_param_grid("small", nthread = 1)
    expect_equal(nrow(g), 8) # max_depth (4) x eta (2)
    expect_setequal(colnames(g), c("max_depth", "eta", "gamma", "colsample_bytree",
                                   "min_child_weight", "subsample", "nthread"))
    expect_setequal(unique(g$max_depth), c(2, 3, 4, 5))
    expect_setequal(unique(g$eta), c(0.01, 0.02))
    expect_true(all(g$min_child_weight == 10) && all(g$subsample == 0.5))
})

test_that("get_param_grid accepts a data.frame as an explicit grid", {
    d <- expand.grid(max_depth = c(3, 4), eta = 0.1)
    g <- get_param_grid(d, nthread = 2)
    expect_equal(nrow(g), 2)
    expect_true(all(g$nthread == 2))
})
