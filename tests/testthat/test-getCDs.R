library(testthat)

# Count rows in the on-disk descriptor cache.
cache_rowcount <- function() {
    con <- cache_connect()
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM cds")$n
}

test_that("getCDs returns correct descriptors and caches new ones on disk", {

    # Known molecules from the shipped cache return correct descriptors.
    X <- RP[1:5, ]
    Y1 <- getCDs(X, verbose = 0)
    nc <- ncol(Y1)
    expect_true(all.equal(Y1[, 1:3], X))
    expect_true(all(colnames(Y1)[4:nc] %in% CDFeatures))
    expect_equal(Y1$Fsp3, c(0.75, 1 / 3, 0.00, 0.50, 0.50))

    n0 <- cache_rowcount()

    # Unknown SMILES (not in the shipped cache) are computed and written back.
    smi <- c("CCC", "CCCC", "CCCCC", "CCC")
    df <- data.frame(NAME = paste0("test", 1:4), SMILES = smi, RT = 1:4)
    Y2 <- getCDs(df, verbose = 0)
    expect_equal(nrow(Y2), 4)
    expect_equal(ncol(Y2), ncol(df) + length(CDFeatures))
    expect_true(all(CDFeatures %in% colnames(Y2)))
    expect_equal(cache_rowcount(), n0 + 3) # 3 distinct new molecules added

    # A second call returns the same result and writes nothing new.
    Y2b <- getCDs(df, verbose = 0)
    expect_equal(Y2b, Y2)
    expect_equal(cache_rowcount(), n0 + 3)

    # Parallel computation yields identical results.
    Y3 <- getCDs(df, verbose = 0, nw = 2)
    expect_equal(Y3, Y2)
})

test_that("getCDs(cache = FALSE) bypasses the cache", {

    n0 <- cache_rowcount()

    # A fresh, never-seen molecule computed with cache = FALSE must not be stored.
    df <- data.frame(NAME = "iso", SMILES = "CC(C)CC(C)C", RT = 1)
    Y <- getCDs(df, verbose = 0, cache = FALSE)
    expect_equal(nrow(Y), 1)
    expect_equal(ncol(Y), ncol(df) + length(CDFeatures))
    expect_true(all(CDFeatures %in% colnames(Y)))
    expect_equal(cache_rowcount(), n0) # nothing written to the cache
})
