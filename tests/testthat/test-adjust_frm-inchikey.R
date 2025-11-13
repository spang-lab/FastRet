library(testthat)

# This test ensures adjust_frm uses INCHIKEY+SMILES when available
# and succeeds even if NAMEs differ between new and old data.

test_that("adjust_frm merges by INCHIKEY when available", {
    # Load pretrained model and clone to avoid side-effects
    frm <- readRDS(pkg_file("extdata/RP_lasso_model.rds"))

    # Add an artificial INCHIKEY column to the training data
    # Ensure non-missing unique keys
    n <- nrow(frm$df)
    frm$df$INCHIKEY <- sprintf("IK%05d", seq_len(n))

    # Build new_data from a small subset and change some NAMES to avoid
    # NAME-matching
    idx <- seq(1, 401, 20)
    new <- frm$df[idx, c("NAME", "SMILES", "RT", "INCHIKEY")]
    new$RT <- new$RT + rnorm(nrow(new), mean = 3, sd = 0.5)  # Slightly perturb RTs

    # Run adjustment; this should use INCHIKEY+SMILES and not error
    afm <- adjust_frm(frm, new, nfolds = 2, verbose = 0, seed = 42, predictors = 1)

    # Basic checks
    expect_true("adj" %in% names(afm))
    expect_true(is.list(afm$adj))
    expect_true("INCHIKEY" %in% colnames(afm$adj$df))
    expect_equal(nrow(afm$adj$df), nrow(new))

    # Now remove some INCHIKEYs from the new data to force fallback to
    # SMILES+NAME mapping
    new$INCHIKEY[8:10] <- NA
    afm2 <- adjust_frm(frm, new, nfolds = 2, verbose = 0, seed = 42, predictors = 1)
    expect_identical(afm$model, afm2$model)

    # Now change a few NAMES as well, so that even the SMILES+NAME fallback
    # cannot find matches for all entries
    new$NAME[8:10] <- NA
    expect_error(
        afm3 <- adjust_frm(frm, new, nfolds = 2, verbose = 0, seed = 42, predictors = 1),
        "Could not map 3 new entries"
    )

    # Plot during interactive sessions
    if (identical(environment(), .GlobalEnv)) {
        par(mfrow = c(3, 1))
        plot(new$RT, frm$df$RT[idx], xlab = "New RT", ylab = "Old RT")
        plot_frm(afm, type = "scatter.train.adj")
        plot_frm(afm, type = "scatter.cv.adj")
        par(mfrow = c(1, 1))
    }
})
