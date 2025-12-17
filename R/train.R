# Fit Models #####

#' @export
#' @keywords public
#'
#' @title Train a new FastRet model (FRM) for retention time prediction
#'
#' @description
#' Trains a new model from molecule SMILES to predict retention times (RT) using
#' the specified method.
#'
#' @param df
#' A dataframe with columns "NAME", "RT", "SMILES" and optionally a set of
#' chemical descriptors. If no chemical descriptors are provided, they are
#' calculated using the function [preprocess_data()].
#'
#' @param method
#' A string representing the prediction algorithm. Either "lasso", "ridge",
#' "gbtree", "gbtreeDefault" or "gbtreeRP". Method "gbtree" is an alias for
#' "gbtreeDefault".
#'
#' @param verbose
#' A logical value indicating whether to print progress messages.
#'
#' @param nfolds
#' An integer representing the number of folds for cross validation.
#'
#' @param nw
#' An integer representing the number of workers for parallel processing.
#'
#' @param degree_polynomial
#' An integer representing the degree of the polynomial. Polynomials up to the
#' specified degree are included in the model.
#'
#' @param interaction_terms
#' A logical value indicating whether to include interaction terms in the model.
#'
#' @param rm_near_zero_var
#' A logical value indicating whether to remove near zero variance predictors.
#'
#' @param rm_na
#' A logical value indicating whether to remove NA values before training.
#' Highly recommended to avoid issues during model fitting. Setting this to
#' FALSE with `method = "lasso"` will most likely lead to errors.
#'
#' @param rm_ns
#' A logical value indicating whether to remove chemical descriptors that were
#' considered as not suitable for linear regression based on a previous analysis
#' of an independent dataset. Currently not used.
#'
#' @param seed
#' An integer value to set the seed for random number generation to allow for
#' reproducible results.
#'
#' @param do_cv
#' A logical value indicating whether to perform cross-validation. If FALSE,
#' the `cv` element in the returned object will be NULL.
#'
#' @return
#' A 'FastRet Model', i.e., an object of class `frm`. Components are:
#' + `model`: The fitted base model. This can be an object of class `glmnet`
#'   (for Lasso or Ridge regression) or `xgb.Booster` (for GBTree models).
#' + `df`: The data frame used for training the model. The data frame contains
#'   all user-provided columns (including mandatory columns RT, SMILES and NAME)
#'   as well the calculated chemical descriptors. (But no interaction terms or
#'   polynomial features, as these can be recreated within a few milliseconds).
#' + `cv`: A named list containing the cross validation results, or NULL if
#'   `do_cv = FALSE`. When not NULL, elements are:
#'   - `folds`: A list of integer vectors specifying the samples in each fold.
#'   - `models`: A list of models trained on each fold.
#'   - `stats`: A list of vectors with RMSE, Rsquared, MAE, pBelow1Min per fold.
#'   - `preds`: Retention time predictions obtained in CV as numeric vector.
#' + `seed`: The seed used for random number generation.
#' + `version`: The version of the FastRet package used to train the model.
#' + `args`: The value of function arguments besides `df` as named list.
#'
#' @examples
#' m <- train_frm(df = RP[1:40, ], method = "lasso", nfolds = 2, verbose = 0)
#' # For the sake of a short runtime, only the first 40 rows of the RP dataset
#' # are used in this example. In practice, you should always use the entire
#' # training dataset for model training.
train_frm <- function(df, method = "lasso", verbose = 1, nfolds = 5, nw = 1,
                      degree_polynomial = 1, interaction_terms = FALSE,
                      rm_near_zero_var = TRUE, rm_na = TRUE, rm_ns = FALSE,
                      seed = NULL, do_cv = TRUE) {

    # Check arguments
    stopifnot(
        is.data.frame(df),
        all(c("NAME", "RT", "SMILES") %in% colnames(df)),
        method %in% c("lasso", "ridge", "gbtree", "gbtreeDefault", "gbtreeRP"),
        is.numeric(nfolds), nfolds >= 2,
        is.numeric(nw), nw >= 1,
        is.numeric(degree_polynomial), degree_polynomial >= 1,
        is.logical(interaction_terms),
        is.logical(rm_near_zero_var),
        is.logical(rm_na),
        is.logical(rm_ns),
        is.null(seed) || (is.numeric(seed) && length(seed) == 1),
        is.logical(do_cv)
    )

    # Init variables
    if (is.numeric(seed)) set.seed(seed)
    if (method == "gbtree") method <- "gbtreeDefault"
    args <- named(
        method, verbose, nfolds, nw, degree_polynomial,
        interaction_terms, rm_near_zero_var,
        rm_na, seed, do_cv
    )
    logf <- if (verbose) catf else null
    dgp <- degree_polynomial
    iat <- interaction_terms
    rm_nzv <- rm_near_zero_var
    df <- preprocess_data(df, 1, FALSE, verbose, nw, FALSE, FALSE)
    # Only add CDs, but don't do any transformations yet, as this is part of
    # model training

    logf("Training a FastRet model with %s base", method)
    dfp <- preprocess_data(df, dgp, iat, verbose, nw, rm_nzv, rm_na)
    meta <- which(colnames(dfp) %in% c("NAME", "SMILES", "RT", "INCHIKEY"))
    M <- dfp[, meta]
    X <- as.matrix(dfp[, -meta])
    y <- M$RT
    model <- if (method %in% c("gbtreeDefault", "gbtreeRP")) {
        xgpar <- if (method == "gbtreeDefault") "default" else "rpopt"
        fit_gbtree(X, y, xgpar, seed, verbose, nw, nfolds, 2000, 1)
    } else {
        fit_glmnet(X, y, method, seed)
    }
    cv <- NULL
    version <- packageVersion("FastRet")
    frm <- named(model, df, cv, seed, version, args)
    frm <- structure(frm, class = "frm")
    logf("Finished training of the FastRet model")

    if (do_cv) {
        logf("Estimating model performance in CV using %d workers", nw)
        folds <- createFolds(seq_len(nrow(df)), k = nfolds)
        train_dfs <- lapply(folds, function(idx) df[-idx, ])
        test_dfs <- lapply(folds, function(idx) df[idx, ])
        verbose <- FALSE
        models <- parLapply2(
            nw, train_dfs, train_frm,
            method = method, verbose = verbose, nfolds = nfolds, nw = 1,
            degree_polynomial = dgp, interaction_terms = iat,
            rm_near_zero_var = rm_nzv, rm_na = rm_na, rm_ns = rm_ns,
            seed = seed, do_cv = FALSE
        )

        m <- models[[1]]
        toscutil::stub(
            predict.frm,
            object = models[[1]],
            df = test_dfs[[1]],
            adjust = NULL,
            verbose = 0,
            clip = TRUE,
            impute = TRUE,
            ... = list()
        )

        preds_per_fold <- mapply(predict.frm, models, test_dfs, SIMPLIFY = FALSE)


        preds <- unname(unlist(preds_per_fold)[order(unlist(folds))])
        stats <- mapply(get_stats, test_dfs, models, SIMPLIFY = FALSE)
        frm$cv <- named(folds, models, stats, preds)
        logf("Finished cross validation")
    }

    frm
}

#' @export
#' @keywords public
#'
#' @title Adjust an existing FastRet model for use with a new column
#'
#' @description
#' The goal of this function is to train a model that predicts RT_ADJ (retention
#' time measured on a new, adjusted column) from RT (retention time measured on
#' the original column) and to attach this adjustment model to an existing
#' FastRet model.
#'
#' @param frm An object of class `frm` as returned by [train_frm()].
#' @param new_data
#' Data frame with required columns "RT", "NAME", "SMILES"; optional "INCHIKEY".
#' "RT" must be the retention time measured on the adjusted column.
#' Each row must match at least one row in `frm$df`.
#' The exact matching behavior is described in 'Details'.
#' @param predictors
#' Numeric vector specifying which transformations to include in the model.
#' Available options are: 1=RT, 2=RT^2, 3=RT^3, 4=log(RT), 5=exp(RT),
#' 6=sqrt(RT). Note that predictor 1 (RT) is always included, even if not
#' specified explicitly.
#' @param nfolds The number of folds for cross validation.
#' @param verbose Show progress messages?
#' @param seed
#' An integer value to set the seed for random number generation to allow for
#' reproducible results.
#' @param do_cv
#' A logical value indicating whether to perform cross-validation. If FALSE,
#' the `cv` element in the returned adjustment object will be NULL.
#' @param adj_type
#' A string representing the adjustment model type. Either "lm", "lasso",
#' "ridge", or "gbtree".
#' @param add_cds
#' A logical value indicating whether to add chemical descriptors as predictors
#' to new data. Default is TRUE if `adj_type` is "lasso", "ridge" or "gbtree"
#' and FALSE if `adj_type` is "lm".
#'
#' @details
#' Matching is done via "SMILES"+"INCHIKEY" if both datasets have non-missing
#' INCHIKEYs for all rows; otherwise via "SMILES"+"NAME". If multiple rows in
#' `frm$df` match the same row in `new_data`, their RT values are averaged
#' first, and this average is used for training the adjustment model.
#'
#' Example: if `frm$df` equals data.frame OLD shown below and `new_data` equals
#' data.frame NEW, then the resulting, paired data.frame will look like PAIRED.
#'
#' ```R
#' OLD <- data.frame(
#'     NAME   = c("A", "B",  "B",  "C"  ),
#'     SMILES = c("C", "CC", "CC", "CCC"),
#'     RT     = c(5.0,  8.0,  8.2,  9.0 )
#' )
#' NEW <- data.frame(
#'     NAME   = c("A", "B",  "B",  "B"),
#'     SMILES = c("C", "CC", "CC", "CC"),
#'     RT     = c(2.5,  5.5,  5.7,  5.6)
#' )
#' PAIRED <- data.frame(
#'     NAME   = c("A", "B",  "B",  "B"),
#'     SMILES = c("C", "CC", "CC", "CC"),
#'     RT     = c(5.0,  8.1,  8.1,  8.1), # Average of OLD$RT[2:3]
#'     RT_ADJ = c(2.5,  5.5,  5.7,  5.6)  # Taken from NEW
#' )
#' ```
#'
#' @return
#' An object of class `frm`, as returned by [train_frm()], but with an
#' additional element `adj` containing the adjustment model. Components of `adj`
#' are:
#'
#' + `model`: The fitted adjustment model. Class depends on `adj_type` and is
#'   one of `lm`, `glmnet`, or `xgb.Booster`.
#'
#' + `df`: The data frame used for training the adjustment model. Including
#'   columns "NAME", "SMILES", "RT", "RT_ADJ" and optionally "INCHIKEY", as well
#'   as any additional predictors specified via the `predictors` argument.
#'
#' + `cv`: A named list containing the cross validation results (see 'Details'),
#'   or NULL if `do_cv = FALSE`. When not NULL, elements are:
#'
#'   - `folds`: A list of integer vectors specifying the samples in each fold.
#'   - `models`: A list of adjustment models trained on each fold.
#'   - `stats`: A list of vectors with RMSE, Rsquared, MAE, pBelow1Min per fold.
#'      Added with v1.3.0.
#'   - `preds`: Retention time predictions obtained during CV by applying the
#'      adjustment model to the hold-out data.
#'   - `preds_adjonly`: Removed (i.e. NULL) since v1.3.0.
#'
#'
#' + `args`: Function arguments used for adjustment (excluding `frm`, `new_data`
#'   and `verbose`). Added with v1.3.0.
#'
#' + `version`: The version of the FastRet package used to train the adjustment
#'   model. Added with v1.3.0.
#'
#' @details
#' If `do_cv` is TRUE, the adjustment procedure is evaluated in
#' cross-validation. However, care must be taken when interpreting the CV
#' results, as the model performance depends on both the adjustment layer and
#' the original model, which was trained on the full base dataset. Therefore,
#' the observed CV metrics should be read as "expected performance when
#' predicting RTs for molecules that were part of the base-model training but
#' not part of the adjustment set" instead of "expected performance when
#' predicting RTs for completely new molecules".
#'
#' @examples
#' frm <- read_rp_lasso_model_rds()
#' new_data <- read_rpadj_xlsx()
#' frm_adj <- adjust_frm(frm, new_data, verbose = 0)
#'
adjust_frm <- function(frm,
                       new_data,
                       predictors = 1:6,
                       nfolds = 5,
                       verbose = 1,
                       seed = NULL,
                       do_cv = TRUE,
                       adj_type = "lm",
                       add_cds = NULL) {

    # Stubs for interactive Development
    if (FALSE) {
        frm <- read_rp_lasso_model_rds()
        new_data <- read_rpadj_xlsx()
        stub(adjust_frm, frm = frm, new_data = new_data, verbose = 2)
    }

    # Check arguments
    stopifnot(
        inherits(frm, "frm"),
        is.data.frame(new_data), all(c("NAME", "SMILES", "RT") %in% colnames(new_data)),
        is.numeric(predictors), all(predictors %in% 1:7), length(predictors) >= 1,
        is.numeric(nfolds), nfolds >= 2,
        is.logical(verbose) || is.numeric(verbose),
        is.null(seed) || is.numeric(seed),
        is.logical(do_cv),
        is.character(adj_type), adj_type %in% c("lm", "lasso", "ridge", "gbtree"),
        is.null(add_cds) || is.logical(add_cds)
    )
    if (length(predictors) == 1 && predictors != 7 && adj_type != "lm") {
        fmt <- "Adjustment via '%s' requires 2+ predictors. Setting 'adj_type' to 'lm'."
        warning(sprintf(fmt, adj_type))
        adj_type <- "lm"
    }
    if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
    if (is.null(add_cds)) add_cds <- if (adj_type == "lm") FALSE else TRUE

    # Configure logging
    logf <- if (verbose >= 1) catf else null
    dbgf <- if (verbose >= 2) catf else null

    # Debug prints
    dbgf("Parameters for model adjustment:")
    dbgf("dim(original_data): %s", paste(dim(frm$df), collapse = " x "))
    dbgf("dim(new_data): %s", paste(dim(new_data), collapse = " x "))
    dbgf("predictors: %s", paste(predictors, collapse = ", "))
    dbgf("nfolds: %s", nfolds)
    dbgf("do_cv: %s", do_cv)
    dbgf("adj_type: %s", adj_type)

    # Map new and old data
    logf("Mapping new and old data")
    args <- named(predictors, nfolds, verbose, seed, do_cv, adj_type)
    df <- merge_dfs(old = frm$df, new = new_data)
    df <- preprocess_data(
        df, degree_polynomial = 1, interaction_terms = FALSE, verbose = verbose,
        nw = 1, rm_near_zero_var = FALSE, rm_na = FALSE, add_cds = add_cds,
        rm_ucs = TRUE, rt_terms = predictors
    )
    frm$adj <- named(model = NULL, df = df, cv = NULL, args = args)

    # Train adjustment model
    logf("Training adjustment model with %s base", adj_type)
    dfp <- preprocess_data(df, 1, FALSE, verbose, 1, TRUE, TRUE, FALSE, FALSE)
    meta <- match(c("NAME", "SMILES", "INCHIKEY", "RT_ADJ"), colnames(df))
    meta <- meta[!is.na(meta)]
    M <- dfp[, meta]
    X <- dfp[, -meta, drop = FALSE]
    y <- M$RT_ADJ
    frm$adj$model <- if (adj_type == "gbtree") {
        fit_gbtree(X, y, "default", seed, verbose, 1, nfolds, 2000, 1)
    } else if (adj_type %in% c("lasso", "ridge")) {
        fit_glmnet(X, y, adj_type, seed)
    } else if (adj_type == "lm") {
        fit_lm(X, y, "RT_ADJ", seed)
    }
    logf("Finished training of adjustment model")

    # Estimate performance in CV
    if (do_cv) {
        logf("Estimating performance of adjustment model in CV")
        folds <- createFolds(seq_len(nrow(df)), k = nfolds)
        train_dfs <- lapply(folds, function(idx) df[-idx, ])
        test_dfs <- lapply(folds, function(idx) df[idx, ])
        verbose <-  max(verbose - 1, 0)
        seeds <- sample.int(.Machine$integer.max, nfolds)
        dbgf("Training adjustment models for %d folds", nfolds)
        models <- mapply(
            adjust_frm, list(frm), train_dfs, list(predictors),
            nfolds, verbose, seeds, FALSE, adj_type,
            SIMPLIFY = FALSE, USE.NAMES = FALSE
        )
        dbgf("Predicting RT_ADJ for hold-out data in each fold")
        preds_per_fold <- mapply(predict.frm,
            models, test_dfs,
            MoreArgs = list(adjust = TRUE, verbose = verbose, clip = TRUE, impute = TRUE),
            SIMPLIFY = FALSE
        )
        preds <- unname(unlist(preds_per_fold)[order(unlist(folds))])
        stats <- mapply(get_stats, test_dfs, models, SIMPLIFY = FALSE)
        frm$adj$cv <- named(folds, models, stats, preds)
        logf("Finished cross validation")
    }

    frm
}

# Use Models #####

#' @export
print.frm <- function(x, ...) {
    model <- if (inherits(x$model, "glmnet")) "glmnet" else "xgboost"
    nr <- nrow(x$df)
    nc <- ncol(x$df)
    v <- as.character(x$version)
    seed <- if (is.null(x$seed)) "NULL" else as.character(x$seed)
    nfold <-  length(x$cv$folds)
    feat <- get_predictors(x, base = TRUE, adjust = FALSE)
    nfeat <- sum(feat != "(Intercept)")
    cv <- if (is.null(x$cv)) {
        "NULL (no cross validation performed)"
    } else {
        paste0(
            sprintf("results of %d-fold cross validation\n", nfold),
            sprintf("  $ folds: list of sample IDs per fold\n"),
            sprintf("  $ models: list of models trained per fold\n"),
            sprintf("  $ stats: list of RMSE, Rsquared, MAE and pBelow1Min per fold\n"),
            sprintf("  $ preds: numeric vector with CV predictions")
        )
    }
    msg <- paste0(
        sprintf("object of class 'frm'\n"),
        sprintf("$ model: %s (num. predictors: %d)\n", model, nfeat),
        sprintf("$ df: %d x %d\n", nr, nc),
        sprintf("$ cv: %s\n", cv),
        sprintf("$ version: %s\n", v),
        sprintf("$ seed: %s\n", seed)
    )
    if (!is.null(x$adj)) {
        cls <- as_str(class(x$adj$model))
        feats <- get_predictors(x, base = FALSE, adjust = TRUE)
        nfeat <- sum(feats != "(Intercept)")
        nr <- nrow(x$adj$df)
        nc <- ncol(x$adj$df)
        cv <- if (is.null(x$adj$cv)) {
            "NULL (no cross validation performed)"
        } else {
            sprintf("results of %d-fold cross validation", length(x$adj$cv$folds))
        }
        adj <- paste0(
            sprintf("$ adj: adjustment info\n"),
            sprintf("  $ model: %s (num. predictors: %s)\n", cls, nfeat),
            sprintf("  $ df: %d x %d\n", nr, nc),
            sprintf("  $ cv: %s\n", cv)
        )
        msg <- paste0(msg, adj)
    }
    cat(msg)
}

#' @export
#' @keywords public
#'
#' @title Predict retention times using a FastRet Model
#'
#' @description
#' Predict retention times for new data using a FastRet Model (FRM).
#'
#' @param object An object of class `frm` as returned by [train_frm()].
#' @param df A data.frame with the same columns as the training data.
#' @param adjust
#' If `object` was adjusted using [adjust_frm()], it will contain a property
#' `object$adj`. If `adjust` is TRUE, `object$adj` will be used to adjust
#' predictions obtained from `object$model`. If FALSE `object$adj` will be
#' ignored. If NULL, `object$model` will be used, if available.
#' @param verbose A logical value indicating whether to print progress messages.
#' @param clip Clip predictions to be within RT range of training data?
#' @param impute Impute missing predictor values using column means of training data?
#' @param ... Not used. Required to match the generic signature of `predict()`.
#'
#' @return
#' A numeric vector with the predicted retention times.
#'
#' @seealso [train_frm()], [adjust_frm()]
#'
#' @examples
#' object <- read_rp_lasso_model_rds()
#' df <- head(RP)
#' yhat <- predict(object, df)
predict.frm <- function(object = train_frm(),
                        df = object$df,
                        adjust = NULL,
                        verbose = 0,
                        clip = TRUE,
                        impute = TRUE,
                        ...) {

    # Load required packages
    pkgs <- get_req_pkgs(object)
    for (pkg in pkgs) withr::local_package(pkg)

    # Init locals
    logf <- if (verbose >= 1) catf else null
    cd_pds <- get_predictors(object, base = TRUE, adjust = FALSE)
    dgp <- get_dgp(cd_pds)
    iat <- any(grepl(":", cd_pds))
    rm_nzv <- FALSE
    rm_na <- FALSE
    add_cds <- TRUE

    # Check arguments
    if (isTRUE(adjust) && is.null(object$adj)) {
        errmsg <- "Model has not been adjusted yet. Please adjust the model first using `adjust_frm()`."
        stop(errmsg)
    }
    if (!all(cd_pds %in% colnames(df))) {
        logf("Chemical descriptors not found in newdata. Trying to calculate them from the provided SMILES.")
        df <- preprocess_data(
            df, dgp, iat, verbose, 1, rm_nzv, rm_na, add_cds,
            rm_ucs = TRUE, rt_terms = 1, mandatory = c("NAME", "SMILES")
        )

    }
    if (!all(cd_pds %in% colnames(df))) {
        missing <- paste(setdiff(cd_pds, colnames(df)), collapse = ", ")
        errmsg <- paste("The following cd_pds are missing in `df`: ", missing)
        stop(errmsg)
    }

    # Impute missing values
    if (impute && any(is.na(df[, cd_pds]))) {
        nap <- cd_pds[colSums(is.na(df[, cd_pds])) > 0]
        napstr <- paste(nap, collapse = ", ")
        logf("NA values found for following cd_pds: %s", napstr)
        logf("Replacing NA values by column means of training data")
        train_df <- object$df
        if (!all(cd_pds %in% colnames(train_df))) {
            # We don't store interaction terms and/or polynomial features in
            # the training data frame, so we need to preprocess it again to get
            # these features.
            train_df <- preprocess_data(
                train_df, dgp, iat, verbose,
                nw = 1, rm_near_zero_var = FALSE, rm_na = FALSE, add_cds = TRUE,
                rm_ucs = TRUE, rt_terms = 1, mandatory = c("NAME", "SMILES")
            )
        }
        for (p in nap) {
            col_mean <- mean(train_df[[p]], na.rm = TRUE)
            df[[p]][is.na(df[[p]])] <- col_mean
        }
    }

    # Predict retention times using base model
    adjust <- !is.null(object$adj) && (isTRUE(adjust) || is.null(adjust))
    logf("Predicting retention times")
    yhat <- as.numeric(predict(object$model, as.matrix(df[, cd_pds])))

    # Adjust predictions if requested
    if (adjust) {
        logf("Adjusting predictions using the adjustment model")
        adj_df <- df
        adj_df$RT <- clip_predictions(yhat, object$df$RT)
        adj_pds <- get_predictors(object, base = FALSE, adjust = TRUE)
        adj_cds <- intersect(adj_pds, CDFeatures)
        adj_rtts <- intersect(adj_pds, RTFeatures)
        missing_cds <- setdiff(adj_cds, colnames(adj_df))
        add_cds <- if (length(missing_cds) > 0) TRUE else FALSE
        if (any(adj_rtts %in% 2:6) || any(adj_pds %notin% colnames(adj_df))) {
            adj_df <- preprocess_data(
                data = adj_df, verbose = verbose, rm_near_zero_var = FALSE,
                rm_na = FALSE, add_cds = add_cds, rm_ucs = FALSE,
                rt_terms = adj_rtts
            )
        }
        adj_df <- adj_df[, adj_pds, drop = FALSE]
        X_adj <- if (inherits(object$adj$model, "lm")) adj_df else as.matrix(adj_df)
        x <- try(yhat <- as.numeric(predict(object$adj$model, X_adj)))
    }

    # Clip predictions if requested
    if (clip) {
        logf("Clipping predictions to be within RT range of training data")
        y <- if (adjust) object$adj$df$RT else object$df$RT
        yhat <- clip_predictions(yhat, y)
    }

    # Return predictions
    as.numeric(yhat)
}

get_req_pkgs <- function(frm) {
    model_types <- get_model_type(frm)
    pkgs <- c(
        if ("glmnet" %in% model_types) "glmnet" else NULL,
        if ("xgboost" %in% model_types) "xgboost" else NULL
    )
    pkgs <- unique(pkgs)
}

get_model_type <- function(frm) {
    bm <- if (inherits(frm$model, "glmnet")) "glmnet"
        else if (inherits(frm$model, "xgb.Booster")) "xgboost"
        else stop("Unknown model type")
    am <- if (is.null(frm$adj)) character(0)
        else if (inherits(frm$adj$model, "glmnet")) "glmnet"
        else if (inherits(frm$adj$model, "xgb.Booster")) "xgboost"
        else if (inherits(frm$adj$model, "lm")) "lm"
        else stop("Unknown adjustment model type")
    c(bm, am)
}

#' @export
#' @title Extract predictor names from an 'frm' object
#' @description Extracts the predictor names from an 'frm' object.
#' @param frm An object of class 'frm' from which to extract the predictor names.
#' @param base Logical indicating whether to include base model predictors.
#' @param adjust Logical indicating whether to include adjustment model predictors.
#' @return A character vector with the predictor names.
#' @keywords internal
#' @examples
#' frm <- read_rp_lasso_model_rds()
#' get_predictors(frm)
get_predictors <- function(frm, base = TRUE, adjust = FALSE) {
    bm <- frm$model
    bp <- if (isFALSE(base)) character(0)
        else if (inherits(bm, "glmnet")) rownames(bm$beta)
        else if (inherits(bm, "lm")) names(stats::coef(bm))[-1]
        else if (inherits(bm, "xgb.Booster")) {
            # For xgboost < 2.0, feature_names is present. For >= 2.0,
            # variable.names() works.
            bm$feature_names %||% variable.names(bm)
        } else {
            stop("Unknown model type")
        }
    am <- frm$adj$model
    ap <- if (isFALSE(adjust) || is.null(am)) character(0)
        else if (inherits(am, "glmnet")) rownames(am$beta)
        else if (inherits(am, "lm")) names(stats::coef(am))[-1]
        else if (inherits(am, "xgb.Booster")) {
            am$feature_names %||% variable.names(am)
        } else {
            stop("Unknown model type")
        }
    unique(c(bp, ap))
}

#' @export
#' @title Clip predictions to observed range
#'
#' @description
#' Clips predicted retention times by fitting a log-normal distribution to the
#' observed training RTs and bounding predictions to the central 99.99%
#' interval. All observed RTs must be positive to estimate the distribution.
#' If the estimated lower bound would be negative, it is replaced by 1% of the
#' observed minimum RT instead.
#'
#' @param yhat Numeric vector of predicted retention times.
#' @param y Numeric vector of observed retention times used to derive bounds.
#'
#' @return Numeric vector of clipped (bounded) predictions.
#'
#' @keywords public
#' @examples
#'
#' # Draw only a few samples (10) and clip based on these. The allowed range will
#' # be much bigger than the observed range.
#'
#' set.seed(42)
#' y <- rlnorm(n = 1000, meanlog = 2, sdlog = 0.1)
#' yhat <- y
#' yhat[1] <- -100 # way too low to be realistic
#' yhat[2] <- 1000 # way too high to be realistic
#' yhat <- clip_predictions(yhat, y)
#' range(y)  # [ 6.18,  8.93]
#' yhat[1:2] # [ 4.96, 10.61] # Limited by theoretical bounds
#'
#'
#' # Draw more samples (1000) and clip based on these. The allowed range will
#' # be almost identical to the observed range.
#'
#' set.seed(42)
#' y <- rnorm(n = 100, mean = 100, sd = 5)
#' yhat <- y
#' yhat[1] <- -100
#' yhat[2] <- 1000
#' yhat <- clip_predictions(yhat, y)
#' range(y)  # 83.14, 117.47
#' yhat[1:2] # 83.14, 117.72
#'
clip_predictions <- function(yhat, y) {
    if (any(y <= 0)) stop("Observed RTs must be strictly positive")
    log_y <- log(y)
    mu <- mean(log_y, na.rm = TRUE)
    sigma <- stats::sd(log_y, na.rm = TRUE)
    sigma <- ifelse(is.na(sigma) || sigma == 0, .Machine$double.eps, sigma)
    lower_bound <- stats::qlnorm(0.0005, meanlog = mu, sdlog = sigma)
    upper_bound <- stats::qlnorm(0.9995, meanlog = mu, sdlog = sigma)
    obsrvd_min <- min(y, na.rm = TRUE)
    obsrvd_max <- max(y, na.rm = TRUE)
    if (lower_bound > obsrvd_min) lower_bound <- obsrvd_min
    if (upper_bound < obsrvd_max) upper_bound <- obsrvd_max
    if (lower_bound < obsrvd_min * 0.01) lower_bound <- obsrvd_min * 0.01
    yhat <- pmin(yhat, upper_bound)
    yhat <- pmax(yhat, lower_bound)
    yhat
}

# Preprocessing #####

#' @noRd
#' @title Merge training and adjustment data.frames
#' @description
#' Merges the original data frame used to train a FastRet model with a new
#' data frame provided by the user for adjusting the original model.
#'
#' @param old The original data frame used to train `frm`. Must contain columns
#' "NAME", "SMILES", "RT" and optionally "INCHIKEY" as well the chemical
#' descriptors listed in `CDFeatures`
#' @param new The new data frame provided by the user. Must contain columns
#' "NAME", "SMILES", "RT" and optionally "INCHIKEY".
#'
#' @return A data frame with following columns:
#' - NAME: Molecule names from `new`
#' - SMILES: Molecule SMILES from `new`
#' - RT: Corresponding retention times taken from `old`
#' - RT_ADJ: Retention times from `new`
#'
#' @examples
#' new <- data.frame(
#'     NAME   = c("A", "B",  "B",  "B"),
#'     SMILES = c("C", "CC", "CC", "CC"),
#'     RT     = c(2.5,  5.5,  5.7,  5.6)
#' )
#' old <- data.frame(
#'      NAME   = c("A", "B",  "B",  "C"),
#'      SMILES = c("C", "CC", "CC", "CCC"),
#'      RT     = c(5.0,  8.0,  8.2,  9.0),
#'      nAtom  = c(1,    2,    2,    3)
#' )
#' merged <- merge_dfs(old, new)
#' print(merged)
#' ## data.frame(
#' ##      NAME = c("A", "B", "B", "B"),
#' ##      SMILES = c("C", "CC", "CC", "CC"),
#' ##      RT_ADJ = c(2.5, 5.5, 5.7, 5.6),
#' ##      RT = c(5, 8.1, 8.1, 8.1),
#' ##      nAtom = c(1, 2, 2, 2)
#' ## )
merge_dfs <- function(old, new) {
    use_inchi <- {
        !is.null(old$INCHIKEY) && all(!is.na(old$INCHIKEY)) &&
        !is.null(new$INCHIKEY) && all(!is.na(new$INCHIKEY))
    }
    if (use_inchi) {
        new <- data.frame(
            NAME = new$NAME,
            INCHIKEY = new$INCHIKEY,
            SMILES = new$SMILES,
            RT_ADJ = new$RT,
            KEY = paste0(new$SMILES, "_", new$INCHIKEY)
        )
        old <- data.frame(
            RT = old$RT,
            KEY = paste0(old$SMILES, "_", old$INCHIKEY)
        )
    } else {
        new <- data.frame(
            NAME = new$NAME,
            SMILES = new$SMILES,
            RT_ADJ = new$RT,
            KEY = paste0(new$SMILES, "_", new$NAME)
        )
        old <- data.frame(
            RT = old$RT,
            KEY = paste0(old$SMILES, "_", old$NAME)
        )
    }
    old <- old[old$KEY %in% unique(new$KEY), ]

    # At this point, we can still have multiple measurements per KEY in old. To
    # deal with this, we average their RTs first, then map every new entry to
    # the average.
    old_RT_avgs <- tapply(old$RT, old$KEY, mean)
    old <- old[!duplicated(old$KEY), ]
    old$RT <- as.numeric(old_RT_avgs[old$KEY]) # as.numeric drops dimnames attr

    # Now finally merge new and old data to get the paired data.frame for model
    # fitting. Make sure to restore the original order after merging, as
    # `merge()` returns rows in an 'unspecified' order (even with sort=FALSE).
    new$ID <- seq_len(nrow(new))
    df <- merge(new, old, by = "KEY", sort = FALSE, all = FALSE)
    df <- df[order(df$ID), ]
    df$ID <- NULL
    df$KEY <- NULL

    if (nrow(df) < nrow(new)) {
        key <- if (use_inchi) "SMILES+INCHIKEY" else "SMILES+NAME"
        fmt <- "Could not map %d new entries to original data based on %s."
        stop(sprintf(fmt, nrow(new) - nrow(df), key))
    } else if (nrow(df) > nrow(new)) {
        # This can't happen due to averaging above. We check anyway to be safe.
        msg <- sprintf("Ambiguous mapping: %d new -> %d old", nrow(new), nrow(df))
        stop(msg)
    }
    df
}

#' @noRd
#' @description Get RMSE, Rsquared, MAE and %below1min for a specific dataset and model.
#' @param data dataframe with retention time in the first column
#' @param model object useable as input for [predict()]
get_stats <- function(df, model) {
    X <- as.matrix(df[, colnames(df) %in% CDFeatures])
    y <- df$RT
    yhat <- predict(model, df, type = "response")
    measures <- c(
        RMSE = sqrt(mean((y - yhat)^2)),
        Rsquared = cor(y, yhat, on_zero_sd = 0)^2,
        MAE = mean(abs(y - yhat)),
        pBelow1Min = sum(abs(y - yhat) < 1.0) / length(y)
    )
    round(measures, 2)
}

#' @noRd
#' @title Get the max degree of polynomial features used in the model
#' @description
#' Extracts the maximum degree of polynomial features from the predictor names.
#' @param x A character vector with predictor names.
#' @return An integer representing the maximum degree of polynomial features.
#' @examples
#' x <- c("nK", "Fsp3", "XLogP", "nK^2", "Fsp3^3")
#' get_dgp(x) # 3
get_dgp <- function(x) {
    splt <- strsplit(x, "^", fixed = TRUE)
    lens <- lengths(splt)
    if (all(lens == 1)) return(1)
    degs <- sapply(splt[lens > 1], function(s) as.numeric(s[2]))
    max(degs, na.rm = TRUE)
}

# Fitting #####

fit_glmnet <- function(X, y = NULL, method = "lasso", seed = NULL) {
    if (is.numeric(seed)) set.seed(seed)
    alpha <- switch(method,
        "lasso" = 1,
        "ridge" = 0,
        stop("method must be 'lasso' or 'ridge', not '", method, "'.")
    )
    cvobj <- glmnet::cv.glmnet(
        x = as.matrix(X), y = y, alpha = alpha, standardize = TRUE,
        family = "gaussian", type.measure = "mse", grouped = FALSE
    )
    model <- glmnet::glmnet(
        x = as.matrix(X), y = y, alpha = alpha, standardize = TRUE,
        family = "gaussian", lambda = cvobj$lambda.min
    )
    model
}

#' @noRd
#' @title Fit a Gradient Boosted Tree model
#'
#' @description
#' Fits multiple GBTree (Gradiant Boosted Tree) models, evaluates their
#' performance in cross validation (CV) and then fits a final model using the
#' optimal set of parameters.
#'
#' @param X Matrix or dataframe with features
#' @param y Numeric vector with target variable
#' @param xgpar
#' Keyword defining how to set xgboost parameters. Available options are:
#' - default: Use default xgboost parameters.
#' - rpopt: Use parameters optimized for RT prediction on the [RP] dataset.
#' - tiny, small, large: Perform a grid search over a tiny, small or large grid
#'   to find the optimal set of parameters. Please note that searching the large
#'   grid can take multiple hours.
#' @param seed Seed for random number generation
#' @param verbose Verbosity level (0 = silent, 1 = progress, 2 = more details).
#' @param nw Number of workers to use for parallel processing.
#' @param nfolds Number of folds for cross validation.
#' @param nrounds Maximum number of boosting rounds.
#' @param nthread Number of threads used internally by xgboost.
#'
#' @return
#' An object of class `xgb.Booster`.
#'
#' @examples
#' df <- preprocess_data(RP, rm_near_zero_var = FALSE)
#' meta <- colnames(df) %in% c("NAME", "SMILES", "RT", "INCHIKEY")
#' M <- df[1:50, m]
#' X <- df[1:50, -m]
#' y <- M$RT
#' gbt <- fit_gbtree(X, y, xgpar = "rpopt", seed = 42)
fit_gbtree <- function(X, y, xgpar = "rpopt", seed = NULL, verbose = 1,
                       nw = 1, nfolds = 10, nrounds = 2000, nthread = 1) {
    if (is.numeric(seed)) set.seed(seed)
    params <- switch(
        xgpar,
        "default" = list(),
        "rpopt" = list(eta = 0.05, max_depth = 4, min_child_weight = 4, subsample = 0.5),
        find_params_best(X, y, xgpar, nfolds, nw, nthread, nrounds, verbose, seed)$best_params
    )
    params$nthread <- nthread
    params$objective <- "reg:squarederror"
    Xmat <- as.matrix(X)
    data <- xgboost::xgb.DMatrix(Xmat, label = y, nthread = nthread)
    xverb <- if (verbose == 2) TRUE else FALSE
    cvobj <- xgboost::xgb.cv(
        params, data, nrounds, nfold = nfolds, early_stopping_rounds = 20,
        verbose = xverb
    )
    nrounds <- cvobj$early_stop$best_iteration %||% cvobj$best_iteration
    xgboost::xgb.train(params, data, nrounds, verbose = xverb)
}

fit_lm <- function(X, y, yname, seed = NULL) {
    if (is.numeric(seed)) set.seed(seed)
    fm <- as.formula(paste(yname, "~", paste(colnames(X), collapse = " + ")))
    Xy <- as.data.frame(X)
    Xy[[yname]] <- y
    stats::lm(formula = fm, data = Xy)
}

RTFeatures <- c("RT", "I(RT^2)", "I(RT^3)", "log(RT)", "exp(RT)", "sqrt(RT)")

# Grid Search #####

#' @noRd
#' @title Find optimal xgboost parameters
#'
#' @description
#' Performs a grid search to find the optimal parameters for training a GBTree
#' model on X to predict y using cross validation.
#'
#' @param X Matrix or dataframe with features
#' @param y Numeric vector with target variable
#' @param searchspace Either "tiny", "small" or "large".
#' @param nfolds Number of folds for cross validation
#' @param nw Number of workers for parallel execution
#' @param nrounds Maximum number of boosting rounds
#' @param verbose Verbosity level
#' @param logf Logging function
#'
#' @return
#' A list with following elements:
#' + `best_params`: A named vector with the optimal parameters
#' + `cv_results`: A dataframe with one row per parame set and columns:
#'   - `nbest`: number of boosting rounds for which the best RMSE was achieved
#'   - `rmse_best`: the best RMSE achieved
#'   - `rmse_std`: standard deviation of the RMSE achieved
#' + `param_grid`: A dataframe with one row per param set that has been tested
#'
#' @examples
#' X <- getCDs(RP, keepdf = FALSE)
#' y <- RP$RT
#' find_params_best(X, y, "tiny", nfolds = 3, nw = 1, nthread = 3, seed = 42)
find_params_best <- function(
    X, y, searchspace = "tiny",
    nfolds = 10, nw = 1, nthread = 1, nrounds = 2000,
    verbose = 1, seed = NULL
) {

    logf <- if (verbose) catf else function(...) {}
    logf("Searching for optimal parameters using Grid Search")
    a <- Sys.time()
    if (is.numeric(seed)) set.seed(seed)
    X_mat <- as.matrix(X)
    param_grid <- get_param_grid(searchspace, nthread)
    foldids <- createFolds(seq_len(nrow(X_mat)), k = nfolds)
    nparams <- nrow(param_grid)
    nw <- min(nw, nparams)

    fmt <- "Evaluating %d param sets in %d-fold CV using %d workers"
    logf(sprintf(fmt, nparams, nfolds, nw))
    cv_result_list <- parLapply2(
        NW = min(nw, nparams),
        ITERABLE = seq_len(nparams),
        EXPORT = c("X_mat", "y", "logf", "param_grid", "nrounds", "foldids", "verbose"),
        BENCHMARK = TRUE,
        FUN = function(i) {
            logf(sprintf("Evaluating parameter set %d/%d", i, nparams))
            cv_obj <- xgboost::xgb.cv(
                params = as.list(param_grid[i, ]),
                data = xgboost::xgb.DMatrix(X_mat, label = y, nthread = nthread),
                nrounds = nrounds,
                folds = foldids,
                early_stopping_rounds = 20,
                objective = "reg:squarederror",
                verbose = if (verbose == 2) TRUE else FALSE
            )
            niter <- cv_obj$niter
            nbest <- cv_obj$best_iteration
            rmse <- cv_obj$evaluation_log$test_rmse_mean[nbest]
            sd <- cv_obj$evaluation_log$test_rmse_std[nbest]
            c(niter = niter, nbest = nbest, rmse = rmse, sd = sd)
        }
    )
    cv_results <- data.frame(do.call(rbind, cv_result_list))
    ntrees <- sum(cv_results$niter)
    secs <- as.numeric(Sys.time() - a, units = "secs")
    logf("Finished evaluation of all parameter sets")
    logf("Total time:  %.2fs", secs)
    logf("Total trees: %d",    ntrees)
    logf("Time / tree: %.2fs", secs / ntrees)
    logf("Time / set:  %.2fs", secs / nparams)

    logf("Finding best parameter set")
    best_index <- which.min(cv_results$rmse)
    best_rmse <- cv_results$rmse[best_index]
    best_params <- unlist(param_grid[best_index, ])
    logf("Best RMSE:   %s", best_rmse)
    logf("Best params: %s", as_str(best_params))

    logf("Finished parameter search in %s", format(Sys.time() - a))
    list(
        best_params = best_params,
        cv_results = cv_results,
        param_grid = param_grid,
        total_secs = secs,
        total_trees = ntrees
    )
}

#' @noRd
#' @description Plot cross validation results of a GBTree model.
#' @param x Object as returned by [fit_gbtree()]
#' @param print Print the plots to the console?
#' @param pdfpath Path to save the plots as PDF
#' @examples
#' X <- head(getCDs(RP, keepdf = FALSE), 20)
#' y <- head(RP$RT, 20)
#' x <- find_params_best(X, y, "tiny", nfolds = 2, seed = 42)
#' plot_params_perf(x)
#' \dontrun{plot_params_perf(x, pdfpath = "misc/cvgbtree.pdf")}
plot_params_perf <- function(x = fit_gbtree(),
                             print = TRUE,
                             pdfpath = NULL,
                             type = "box") {
    .data <- ggplot2::.data
    df <- cbind(x$param_grid, x$cv_results)
    params <- colnames(x$param_grid)
    df[, params] <- lapply(df[, params, drop = FALSE], as.factor)
    df$run <- seq_len(nrow(df))
    # Convert to long format
    dfs <- lapply(params, function(p) {
        data.frame(run = df$run, param = p, value = df[[p]], rmse = df$rmse)
    })
    ggdf <- do.call(rbind, dfs)
    mapping <- ggplot2::aes(x = .data[["value"]], y = .data[["rmse"]])
    p <- ggplot2::ggplot(ggdf, mapping)
    p <- p + ggplot2::facet_wrap(~ .data[["param"]], scales = "free_x")
    p <- p + ggplot2::labs(x = "value", y = "rmse")
    if (type == "violin") p <- p + ggplot2::geom_violin(trim = FALSE)
    if (type == "box") p <- p + ggplot2::geom_boxplot()
    if (!is.null(pdfpath)) ggplot2::ggsave(pdfpath, p, "pdf")
    if (print) print(p)
    invisible(p)
}

#' @noRd
#' @title Get parameter grid
#' @description
#' Get the parameter grid for [find_params_best()].
#' @param size Size of the param grid. Either "tiny", "small" or "large".
#' @return A dataframe with a unique combination of parameters per row.
get_param_grid <- function(size = "large", nthread = NULL) {
    space <- if (size == "large") list(
        max_depth = 1:6, # max depth of a tree
        eta = c(0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40), # learning rate
        gamma = 0:3, # min loss reduction to allow further partition
        colsample_bytree = (8:10) / 10, # subsample ratio of columns when constructing each tree
        subsample = (4:10) / 10, # subsample ratio of the training instance
        min_child_weight = 1:5 # minimum number of instances needed to be in each node
    ) else if (size == "small") list(
        max_depth = 3:6,
        eta = c(0.10, 0.20, 0.30, 0.40)
    ) else if (size == "tiny") list(
        eta = c(0.10, 0.20, 0.30, 0.40)
    ) else {
        stop("Invalid grid size: ", size)
    }
    space$nthread = nthread # threads used internally by xgboost
    expand.grid(space, stringsAsFactors = FALSE)
}

#' @noRd
#' @description
#' Calls [find_params_best()] multiple times with different numbers of workers
#' and threads to benchmark runtime of the grid search with respect to these two
#' parameters.
#' @examples
#'
#' \dontrun{benchmark_find_params}
#'
#' # Param Sets: 16 - Total trees: 1490
#' #
#' # |   w/t   secs  |  w/t   secs |  w/t   secs | w/t    secs |
#' # | ------------- | ----------- | ----------- | ----------- |
#' # |   1/1   22.5  |   -       - |   -      -  |  -       -  |
#' # |   2/1   16.5  |  1/2   17.4 |   -      -  |  -       -  |
#' # |   4/1    9.8  |  2/2   12.4 |  1/4   14.5 | 1/8    16.3 |
#' # |   8/1    8.2  |  4/2    8.9 |  2/4   12.8 | 2/8    14.0 |
#' # |  16/1    8.4  |  8/2    8.2 |  4/4    9.0 | 1/16   24.0 |
#' #
#' # Summary: using lots of threads (t) doesn't help if the number of trees per
#' # model is low. Using lots of workers (w) is more efficient. That's nice,
#' # because it means implementing multi-processing was not a waste of time.
#'
benchmark_find_params <- function() {
    .data <- ggplot2::.data
    RP <- FastRet::RP
    steps <- c(1, 2, 4, 8, 16)
    pgrid <- expand.grid(nw = steps, nthread = steps, KEEP.OUT.ATTRS = FALSE)
    pgrid <- pgrid[(pgrid$nw * pgrid$nthread) %in% steps, ]
    ngrid <- nrow(pgrid)
    X <- getCDs(RP, keepdf = FALSE)
    y <- RP$RT
    results <- lapply(seq_len(nrow(pgrid)), function(i) {
        nw <- pgrid$nw[i]
        nthread <- pgrid$nthread[i]
        catf("Benchmark %d/%d (nw=%d, nthread=%d)", i, ngrid, nw, nthread)
        find_params_best(X, y, "small", 4, nw, nthread, verbose = FALSE, seed = 42)
    })
    secs <- sapply(results, `[[`, "total_secs")
    n_params <- nrow(results[[1]]$param_grid)
    n_trees <- results[[1]]$total_trees
    main <- sprintf("Param sets: %d - Total trees: %s", n_params, n_trees)
    df <- data.frame(nw = pgrid$nw, nthread = pgrid$nthread, secs = secs)
    mapping <- ggplot2::aes(x = factor(.data$nw), y = factor(.data$nthread), fill = secs)
    p <- ggplot2::ggplot(df, mapping)
    p <- p + ggplot2::geom_tile()
    p <- p + ggplot2::scale_fill_viridis_c(name = "secs")
    p <- p + ggplot2::labs(x = "nw", y = "nthread", title = main)
    p <- p + ggplot2::theme_minimal()
    print(p)
    df
}
