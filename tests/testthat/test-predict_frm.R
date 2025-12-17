library(testthat)

test_that("predict.frm returns numeric vector (unadjusted)", {
    frm <- read_rp_lasso_model_rds()
    df <- frm$df[1:5, ]
    yhat <- predict(frm, df, adjust = FALSE, verbose = 0)
    expect_type(yhat, "double")
    expect_length(yhat, nrow(df))
})

test_that("predict.frm applies adjustment when available", {
    skip_on_cran()
    frm <- read_rp_lasso_model_rds()
    new_data <- read_rpadj_xlsx()
    frm2 <- adjust_frm(frm, new_data, predictors = 1:2, nfolds = 2, verbose = 0)
    df <- frm$df[1:5, ]
    yhat_adj <- predict(frm2, df, adjust = TRUE, verbose = 0)
    expect_type(yhat_adj, "double")
    expect_length(yhat_adj, nrow(df))
})

test_that("predict.frm imputes CDs if they are NA", {

    # Handcraft a frm model with only three descriptors, where one descriptor is
    # 'Kier3', as this one is sometimes NA (checked manually).
    frm <- read_rp_lasso_model_rds()
    frm$model$beta <- new(
        "dgCMatrix", i = c(2L), p = c(0L, 1L), Dim = c(3L, 1L),
        Dimnames = list(c("Fsp3", "nSmallRings", "Kier3"), "s0"),
        x = c(0.83), factors = list()
    )
    frm$model$df <- 1L
    frm$model$df <- 1L
    frm$model$dim <- c(3L, 1L)

    # Now handcraft new_data where Kier3 is NA for one entry (THIOUREA)
    new_data <- structure(
        list(
            NAME = c("THIOUREA", "L-SELENOMETHIONINE"),
            SMILES = c( "S=C(N)N", "O=C(O)C(N)CC[Se]C"),
            RT = c(1.3, 5.77),
            INCHIKEY = c( "UMGDCJDMYOKAJW-UHFFFAOYSA-N", "RJFAYQIBOAGBLC-BYPYZUCNSA-N")
        ),
        row.names = c( "HILIC_89", "HILIC_214"),
        class = "data.frame"
    )
    preproc_data <- preprocess_data(new_data, rm_near_zero_var = FALSE, rm_na = FALSE, verbose = FALSE)

    # Check Test Assumptions
    expect_equal(coef(frm$model)["Kier3", "s0"], 0.83)
    expect_true(any(is.na(preproc_data$Kier3)))

    # Now predict with and without imputation
    yhat <- predict.frm(frm, new_data)
    yhat_no_imputing <- predict.frm(frm, new_data, impute = FALSE)

    # Check Results
    expect_equal(round(yhat, 2), c(4.31, 4.54))
    expect_equal(round(yhat_no_imputing, 2), c(NaN, 4.54))
})
