library(testthat)

test_fit_gbtree <- function() {
    set.seed(123)
    n <- 100
    Fsp3 <- rnorm(n)
    nSmallRings <- rnorm(n)
    nAromRings <- rnorm(n)
    noise <- rnorm(n)
    RT <- Fsp3^2 + sin(nSmallRings) + noise

    df <- data.frame(RT, Fsp3, nSmallRings, nAromRings)
    X <- df[, c("Fsp3", "nSmallRings", "nAromRings")]
    y <- df$RT

    # Perform test. Use fully qualified function names, so we can run this in a
    # separate R process.
    result <- FastRet:::fit_gbtree(X, y, verbose = 0)
    testthat::expect_true(inherits(result, "xgb.Booster"))
}

test_that("fit_gbtree works", {
    test_fit_gbtree()
})

skip_on_cran() # Test is very slow
skip_if(packageVersion("xgboost") < package_version("2.0.0")) # Older version already tested above

# Installs xgboost v1.7.9.1 (binary if possible, source otherwise) into new_lib
install_xgboost_1.7.9.1 <- function(new_lib) {
    binary_url <- "https://cran.r-project.org/bin/windows/contrib/4.3/xgboost_1.7.9.1.zip"
    source_url  <- paste0("https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.7.11.1.tar.gz")
    suppressMessages(tryCatch(
        install.packages(binary_url, type = "binary", lib = new_lib, repos = NULL, quiet = TRUE),
        error = function(e) install.packages(source_url, type = "source", lib = new_lib, repos = NULL, quiet = TRUE)
    ))
    testthat::expect_true(
        packageVersion("xgboost", lib.loc = new_lib) == package_version("1.7.9.1"),
        info = paste("Installed xgboost version:", as.character(packageVersion("xgboost", lib.loc = new_lib)))
    )
}

test_that("fit_gbtree works with xgboost 1.7", {
    new_lib <- file.path(tempdir(), "new_lib")
    unlink(new_lib, recursive = TRUE)
    dir.create(new_lib, showWarnings = FALSE)
    callr::r(
        func = function(new_lib, test_fit_gbtree, install_xgboost_1.7.9.1) {
            .libPaths(c(new_lib, .libPaths()))
            install_xgboost_1.7.9.1(new_lib)
            testthat::expect_true(
                packageVersion("xgboost") == package_version("1.7.9.1"),
                info = paste("Installed xgboost version:", as.character(packageVersion("xgboost")))
            )
            test_fit_gbtree()
        },
        args = list(new_lib, test_fit_gbtree, install_xgboost_1.7.9.1)
        # Rationale for using a separate process: loading xgboost modifies R's
        # search path, loaded DLLs, etc. and it turns out restoring the old
        # state is quite tricky. By using a separate R process we can avoid all
        # that hassle.
    )
    unlink(new_lib, recursive = TRUE)
})
