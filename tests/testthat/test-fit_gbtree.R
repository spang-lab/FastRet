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

test_that("fit_gbtree works with xgboost 1.7", {

    # Since installation of packages on CRAN servers is not allowed, we skip
    # this test if on CRAN. (In theory, we should be fine, because we install in
    # a temporary library and don't modify the user's search path, but still,
    # there are lots of things that could go wrong, like too much CPU time spent
    # in a subprocess. In the end, it's not worth going through all that
    # hassle.)
    skip_on_cran()

    # We also skip if the currently installed version of xgboost is already
    # below 2.0. In this case the test above ('fit_gbtree works') has already
    # verified that fit_gbtree works with such versions.
    skip_if(packageVersion("xgboost") < package_version("2.0.0"), "xgboost version is below 2.0")

    # And lastly, we also skip if we're not on Windows or macOS, because we only
    # have binaries for these platforms on CRAN and installing from source takes
    # too long (several minutes).
    skip_if(!(Sys.info()[["sysname"]] %in% c("Windows", "Darwin")), "Binary-install of xgboost not supported on this OS")

    # Function to install xgboost version 1.7.x into new_lib from binary.
    install_xgboost_v1 <- function(new_lib) {
        sys_info <- Sys.info()
        os <- sys_info[["sysname"]]
        architecture <- sys_info[["machine"]]
        if (os == "Windows") {
            url <- "https://cran.r-project.org/bin/windows/contrib/4.3/xgboost_1.7.9.1.zip"
            version <- "1.7.9.1"
        } else if (os == "Darwin" && architecture == "x86_64") {
            url <- "https://cran.r-project.org/bin/macosx/big-sur-x86_64/contrib/4.3/xgboost_1.7.11.1.tgz"
            version <- "1.7.11.1"
        } else if (os == "Darwin" && architecture == "arm64") {
            url <- "https://cran.r-project.org/bin/macosx/big-sur-arm64/contrib/4.3/xgboost_1.7.11.1.tgz"
            version <- "1.7.11.1"
        } else {
            fmt <- "Unsupported OS/architecture combination for installing xgboost binary: %s/%s"
            stop(sprintf(fmt, os, architecture))
        }
        message("Installing:", url)
        install.packages(url, type = "binary", lib = new_lib, repos = NULL, quiet = TRUE)
            testthat::expect_true(
            packageVersion("xgboost", lib.loc = new_lib) == package_version(version),
            info = paste("Installed xgboost version:", as.character(packageVersion("xgboost", lib.loc = new_lib)))
        )
    }

    new_lib <- file.path(tempdir(), "new_lib")
    unlink(new_lib, recursive = TRUE)
    dir.create(new_lib, showWarnings = FALSE)
    testthat::expect_no_error(callr::r(
        func = function(new_lib, test_fit_gbtree, install_xgboost_v1) {
            .libPaths(c(new_lib, .libPaths()))
            install_xgboost_v1(new_lib)
            testthat::expect_true(
                packageVersion("xgboost") < package_version("2.0.0"),
                info = paste("Installed xgboost version:", as.character(packageVersion("xgboost")))
            )
            test_fit_gbtree()
        },
        args = list(new_lib, test_fit_gbtree, install_xgboost_v1)
        # Rationale for using a separate process: loading xgboost modifies R's
        # search path, loaded DLLs, etc. and it turns out restoring the old
        # state is quite tricky. By using a separate R process we can avoid all
        # that hassle.
    ))
    unlink(new_lib, recursive = TRUE)
})
