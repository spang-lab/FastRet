library(testthat)

test_that("adjust_frm merges by INCHIKEY when available", {

    # Load pretrained model for speed-up
    frm <- readRDS(pkg_file("extdata/RP_lasso_model.rds"))

    # Add an artificial INCHIKEY column to the training data. Then introduce
    # some duplicate SMILES+INCHIKEY entries with slightly different RTs.
    n <- nrow(frm$df)
    frm$df$INCHIKEY <- sprintf("IK%05d", seq_len(n))
    frm$df$INCHIKEY[1:5] <- frm$df$INCHIKEY[1]
    frm$df$SMILES[1:5] <- frm$df$SMILES[1]
    frm$df$RT[1:5] <- frm$df$RT[1] + c(0, 0.1, 0.2, 0.3, 0.4)

    # Build new_data from a small subset and change some NAMES to avoid
    # NAME-matching. Duplicate the first entry, so we have the most complex
    # case. Multiple new entries (2) matching multiple old entries (5). Also
    # make sure that the duplicated entry is at different positions within
    # new_data.
    idx <- c(seq(1, 401, 40), 1)
    new <- frm$df[idx, c("NAME", "SMILES", "RT", "INCHIKEY")]
    new$RT <- new$RT + rnorm(nrow(new), mean = 3, sd = 0.5)  # Slightly perturb RTs

    # Run adjustment; this should use INCHIKEY+SMILES and not error
    afm <- adjust_frm(frm, new, nfolds = 2, verbose = 0, seed = 42, predictors = 1)
    afm42 <- adjust_frm(frm, new, nfolds = 2, verbose = 0, seed = 42, predictors = 1)
    afm99 <- adjust_frm(frm, new, nfolds = 2, verbose = 0, seed = 99, predictors = 1)

    # Basic checks
    expect_true("adj" %in% names(afm))
    expect_true(is.list(afm$adj))
    expect_true("INCHIKEY" %in% colnames(afm$adj$df))
    expect_equal(nrow(afm$adj$df), nrow(new))
    expect_true(isTRUE(all.equal(afm, afm42)))
    expect_false(isFALSE(all.equal(afm, afm99)))
    expect_equal(new$RT, afm$adj$df$RT_ADJ)

    # Ensure that afm$adj$df$RT is calculated correctly by mapping to frm$df via
    # INCHIKEY+SMILES and averaging RT over duplicates.
    n <- length(idx)
    old <- afm$df[idx, ] # (1)
    old$RT[c(1, n)] <- mean(afm$df$RT[1:5]) # (2)
    expect_equal(length(unique(c(new$SMILES[c(1, n)], afm$df$SMILES[1:5]))), 1) # (3)
    expect_equal(length(unique(c(new$INCHIKEY[c(1, n)], afm$df$INCHIKEY[1:5]))), 1) # (3)
    expect_equal(afm$adj$df$RT_ADJ, new$RT)
    expect_equal(afm$adj$df$RT, old$RT)
    # (1) We generated `new` from `afm$df[idx, ]`, so naivly, we would expect
    # these elements to be picked as corresponding "old" measurements, that will
    # be adjusted.
    # (2) However, because the first and last entry in `new` have the same
    # INCHIKEY+SMILES combination as the first five entries in `afm$df` (3), our
    # algorithm should calculate the average retention time over `afm$df[1:5,]`
    # and then map both `new[1,]` and `new[n,]` to that average.

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

test_that("adjust_frm works with do_cv = FALSE", {
    # Load pretrained model for speed-up
    frm <- readRDS(pkg_file("extdata/RP_lasso_model.rds"))
    
    # Build new_data from a small subset
    idx <- seq(1, 401, 40)
    new <- frm$df[idx, c("NAME", "SMILES", "RT")]
    new$RT <- new$RT + rnorm(nrow(new), mean = 3, sd = 0.5)
    
    # Test with do_cv=FALSE
    afm1 <- adjust_frm(frm, new, nfolds = 2, verbose = 0, seed = 42, predictors = 1, do_cv = FALSE)
    
    # Test with do_cv=TRUE
    afm2 <- adjust_frm(frm, new, nfolds = 2, verbose = 0, seed = 42, predictors = 1, do_cv = TRUE)
    
    # Model with do_cv=FALSE should have NULL cv element
    expect_null(afm1$adj$cv)
    expect_true(is.list(afm1$adj))
    expect_equal(names(afm1$adj), c("model", "df", "cv", "args"))
    
    # Model with do_cv=TRUE should have cv element
    expect_false(is.null(afm2$adj$cv))
    expect_true(is.list(afm2$adj$cv))
    expect_equal(names(afm2$adj$cv), c("folds", "models", "preds", "preds_adjonly"))
    
    # The actual adjustment models should be the same (only cv differs)
    expect_equal(coef(afm1$adj$model), coef(afm2$adj$model))
    expect_equal(afm1$adj$df, afm2$adj$df)
})
