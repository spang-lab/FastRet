library(testthat)

skip_on_ci()
skip_on_cran()

test_that("parLapply2 works correctly", {

    # Helpers
    square_fast <- function(x) { x^2 }
    square_slow <- function(x) { Sys.sleep(0.2) ; x^2 }
    rt <- function(expr) { system.time(expr)[[3]] }

    # Inputs
    x <- 1:4
    y <- as.list(x^2)

    # Call with 1 and 2 cores
    rt1 <- rt(y1 <- lapply(x, square_fast))
    rt2 <- rt(y2 <- parLapply2(NW = 1, ITERABLE = x, FUN = square_fast))
    rt3 <- rt(y3 <- parLapply2(NW = 2, ITERABLE = x, FUN = square_fast))
    rt4 <- rt(y4 <- parLapply2(NW = 2, ITERABLE = x, FUN = square_slow))

    # Expect correct results in all cases
    expect_equal(y1, y)
    expect_equal(y2, y)
    expect_equal(y3, y)
    expect_equal(y4, y)

    skip() # Runtime checks are always fragile, so skip them unless executed manually

    # Runtime checks
    cct <- rt(stopCluster(makeCluster(1))) # Time to create a cluster
    expect_true(rt1 <= 0.01) # Should be almost instant
    expect_true(rt2 <= 0.02) # Same but a little slower due to overhead
    expect_true(rt3 <= 2 * cct) # (1)
    # (1) Should be Cluster Creation time + Package Loading Time + Function
    # Execution Time + Function Overhead. We assume Cluster Creation is slower
    # than the other three combined, so we set a maximum of 2*cct.
    expect_true(rt4 <= rt3 + 0.6) # (2)
    # (2) We call Sys.sleep(0.2) four times, so if were sequential, we would
    # expect rt4 to be rt3 + 0.8s. However, since we use two cores, we expect
    # the rt4 to be roughly rt3 + 0.4s. To be on the safe side, we set the
    # maximum to rt3 + 0.6s.
})
