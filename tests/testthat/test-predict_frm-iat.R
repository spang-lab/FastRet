library(testthat)

# This test ensures that prediction works with a model trained using
# degree_polynomial >= 2.

test_that("predict.frm works with interaction terms", {
    testthat::with_mocked_bindings(
        CDFeatures = c(
            "nRings9", "Zagreb", "nB", "VAdjMat", "topoShape", "nW", "MLogP"
            # Only use a very samll set of CDs for testing or fitting with all
            # interaction terms will become very slow (> 30s for the whole
            # test). Above CDs were chosen by first fitting a model on all CDs
            # interactively and then selecting only a few that were part of the
            # final model (so we will get an acceptable fit). Obviously, from a
            # statistical point of view this is not correct, but we're only
            # testing software functionality here, not prediction accuracy!
        ),
        code = {
            frm <- train_frm(
                df = RP[1:300, ],
                method = "lasso",
                nfolds = 2,
                nw = 1,
                verbose = 0,
                degree_polynomial = 1,
                interaction_terms = TRUE,
                rm_near_zero_var = FALSE,
                rm_na = TRUE
            )

            # Ensure at least one interaction term was part of the trained model
            preds <- get_predictors(frm)
            expect_true(any(grepl(":", preds, fixed = TRUE)))

            # Predict on holdout
            new_df <- RP[301:422, ]
            y <- new_df$RT
            yhat <- predict(frm, new_df, adjust = FALSE, verbose = 0)
            plot(y, yhat)
            expect_type(yhat, "double")
            expect_length(yhat, nrow(new_df))
            expect_true(all(is.finite(yhat)))
        }
    )
})
