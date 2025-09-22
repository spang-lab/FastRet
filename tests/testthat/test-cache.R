library("testthat")

#' @noRd
#' @title Check State of RAM Cache
#' @description Check that RAM cache has expected structure and content.
#' @param nsx Number of SMILES expected
check_state_of_ram_cache <- function(nsx = NULL) {
    expect_true(is.environment(ram_cache))
    expect_true(length(ram_cache) == 1)
    expect_true(names(ram_cache) == "CDs")
    expect_true(is.data.frame(ram_cache$CDs))
    if (!is.null(nsx)) {
        expect_true(nrow(ram_cache$CDs) == nsx)
    }
}

test_that("RAM Cache Functions Work", {
    rc_clear()
    check_state_of_ram_cache(nsx = 0)
    rc_set("CCC", getCDsFromCDK("CCC"))
    check_state_of_ram_cache(nsx = 1)
    rc_init()
    check_state_of_ram_cache(nsx = 1317)
    expect_equal(rc_get("CCC"), getCDsFromCDK("CCC"))
})

test_that("Disk Cache Functions Work", {
    # get_cache_dir
    # dc_path
    # dc_clear
    # dc_get
    # dc_set
})

skip()

test_that("update_cachedata_cds works", {
    CDs_calc <- update_cachedata_cds(overwrite = FALSE, nw = 1, nsmi = 4)
    CDs_read <- readRDS(pkg_file("cachedata/CDs.rds"))
    expect_equal(CDs_calc, CDs_read[1:4, ])
})
