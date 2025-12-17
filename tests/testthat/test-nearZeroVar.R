library(testthat)

test_that("nearZeroVar works", {
    # Variables are filtered out if:
    # 1. freqRatio > freqCut (= 95/5 by default) AND
    # 2. percUniq <= uniqueCut (= 10 by default)
    # Where:
    # - freqRatio = frequency of most common value (MCV) / freq. of 2nd MCV
    # - percUniq = ((number of unique values) / (number of samples)) * 100
    df <- data.frame(
        both_false = 1:100,
        both_true = c(rep(1, 96), rep(2, 4)),
        only_unique = c(rep(1, 95), rep(2, 5)),
        only_freq = c(rep(0, 60), 1:40)
    )
    # Check our own implementation
    expect_equal(nearZeroVar(df), 2)
    # Check against original caret implementation
    skip_if_not(requireNamespace("caret", quietly = TRUE))
    expect_equal(caret::nearZeroVar(df), 2)
})

