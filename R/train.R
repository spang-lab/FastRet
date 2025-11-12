# Fitting #####

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
#' "gbtreeRP".
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
#' @return
#' A 'FastRet Model', i.e., an object of class `frm`. Components are:
#' + `model`: The fitted base model. This can be an object of class `glmnet`
#'   (for Lasso or Ridge regression) or `xgb.Booster` (for GBTree models).
#' + `df`: The data frame used for training the model. The data frame contains
#'   all user-provided columns (including mandatory columns RT, SMILES and NAME)
#'   as well the calculated chemical descriptors. (But no interaction terms or
#'   polynomial features, as these can be recreated within a few milliseconds).
#' + `cv`: A named list containing the cross validation results. Elements are:
#'   - `folds`: A list of integer vectors specifying the samples in each fold.
#'   - `models`: A list of models trained on each fold.
#'   - `stats`: A list of vectors with RMSE, Rsquared, MAE, pBelow1Min per fold.
#'   - `preds`: Retention time predictions obtained in CV as numeric vector.
#' + `seed`: The seed used for random number generation.
#' + `version`: The version of the FastRet package used to train the model.
#'
#' @examples
#' m <- train_frm(df = RP[1:40, ], method = "lasso", nfolds = 2, verbose = 0)
#' # For the sake of a short runtime, only the first 40 rows of the RP dataset
#' # are used in this example. In practice, you should always use the entire
#' # training dataset for model training.
train_frm <- function(df, method = "lasso", verbose = 1, nfolds = 5, nw = 1,
                      degree_polynomial = 1, interaction_terms = FALSE,
                      rm_near_zero_var = TRUE, rm_na = TRUE, rm_ns = FALSE,
                      seed = NULL) {

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
        is.null(seed) || (is.numeric(seed) && length(seed) == 1)
    )

    # Init variables
    if (is.numeric(seed)) set.seed(seed)
    if (method == "gbtree") method <- "gbtreeRP"
    logf <- if (verbose) catf else null
    dgp <- degree_polynomial
    iat <- interaction_terms
    rm_nzv <- rm_near_zero_var
    df <- preprocess_data(df, 1, FALSE, verbose, nw, FALSE, FALSE)
    # Only add CDs, but don't do any transformations yet, as this is part of
    # model training

    # Train FastRet model on full data
    logf("Training a FastRet model with %s base", method)
    frm <- train_frm_internal(
        df,
        method, verbose, nfolds, nw, dgp, iat, rm_nzv, rm_na, seed
    )
    logf("Finished training of the FastRet model")

    # Estimate performance in cross validation
    logf("Estimating model performance in CV using %d workers", nw)
    folds <- createFolds(seq_len(nrow(df)), k = nfolds)
    train_dfs <- lapply(folds, function(idx) df[-idx, ])
    test_dfs <- lapply(folds, function(idx) df[idx, ])
    verbose <- FALSE
    models <- parLapply2(
        nw, train_dfs, train_frm_internal,
        method, verbose, nfolds, nw, dgp, iat, rm_nzv, rm_na, seed
    )
    preds_per_fold <- mapply(predict, models, test_dfs, SIMPLIFY = FALSE)
    preds <- unname(unlist(preds_per_fold)[order(unlist(folds))])
    stats <- mapply(get_stats, test_dfs, models, SIMPLIFY = FALSE)
    frm$cv <- named(folds, models, stats, preds)
    logf("Finished cross validation")
    frm
}

#' @noRd
#' @title Trains a FastRet model without cross validation
#' @description
#' Trains a `frm` model, but does NOT estimate performance in CV.
train_frm_internal <- function(df, method, verbose, nfolds, nw, dgp, iat,
                               rm_nzv, rm_na, seed) {

    ## Prepare data for model fitting
    logf <- if (verbose) catf else null
    dfp <- preprocess_data(df, dgp, iat, verbose, nw, rm_nzv, rm_na)
    meta <- which(colnames(dfp) %in% c("NAME", "SMILES", "RT", "INCHIKEY"))
    M <- dfp[, meta]
    X <- as.matrix(dfp[, -meta])
    y <- M$RT

    # Start model fitting
    logf("Training a %s model", method)
    pkg <- if (grepl("gbtree", method)) "xgboost" else "glmnet"
    withr::local_package(pkg)
    model <- if (method %in% c("gbtreeDefault", "gbtreeRP")) {
        xpar <- if (method == "gbtreeDefault") "default" else "rpopt"
        fit_gbtree(X, y, xpar, seed, verbose, nw, nfolds, 2000, 1)
    } else {
        alpha <- switch(method, "lasso" = 1, "ridge" = 0, NA)
        fit_glmnet(X, y, alpha)
    }
    model$feature_names <- colnames(X) # Required prediction of new data
    logf("Finished training of the %s model", method)

    cv <- NULL
    version <- packageVersion("FastRet")
    frm <- named(model, df, cv, seed, version)
    frm <- structure(frm, class = "frm")
}

#' @export
#' @keywords public
#'
#' @title Adjust an existing FastRet model for use with a new column
#'
#' @description
#' The goal of this function is to train a model that predicts RT_ADJ (retention
#' time measured on a new, adjusted column) from RT (retention time measured on
#' the original column) and to attach this "adjustmodel" to an existing FastRet
#' model.
#'
#' @param frm
#' An object of class `frm` as returned by [train_frm()].
#'
#' @param new_data
#' Dataframe with columns "RT", "NAME", "SMILES" and optionally a set of
#' chemical descriptors.
#'
#' @param predictors
#' Numeric vector specifying which predictors to include in the model in
#' addition to RT. Available options are: 1=RT, 2=RT^2, 3=RT^3, 4=log(RT),
#' 5=exp(RT), 6=sqrt(RT).
#'
#' @param nfolds
#' An integer representing the number of folds for cross validation.
#'
#' @param verbose
#' A logical value indicating whether to print progress messages.
#'
#' @return
#' An object of class `frm`, which is a list with the following elements:
#'
#' * `model`: A list containing details about the original model.
#' * `df`: The data frame used for training the model.
#' * `cv`: A list containing the cross validation results.
#' * `seed`: The seed used for random number generation.
#' * `version`: The version of the FastRet package used to train the model.
#' * `adj`: A list containing details about the adjusted model.
#'
#' @examples
#' frm <- read_rp_lasso_model_rds()
#' new_data <- read_rpadj_xlsx()
#' frmAdjusted <- adjust_frm(frm, new_data, verbose = 0)
adjust_frm <- function(frm = train_frm(),
                       new_data = read_rpadj_xlsx(),
                       predictors = 1:6,
                       nfolds = 5,
                       verbose = 1) {
    if (!is.numeric(predictors) || length(predictors) < 1 || !all(predictors %in% 1:6)) {
        stop("Invalid predictors. Please provide a vector of integers between 1 and 6.")
    }
    if (isFALSE(verbose) || verbose == 0) catf <- function(...) {}

    catf("Starting model Adjustment")
    catf("dim(original_data): %s", paste(dim(frm$df), collapse = " x "))
    catf("dim(new_data): %s", paste(dim(new_data), collapse = " x "))
    catf("predictors: %s", paste(predictors, collapse = ", "))
    catf("nfolds: %s", nfolds)

    catf("Preprocessing data")
    new <- data.frame(NAME = new_data$NAME, SMILES = new_data$SMILES, RT_ADJ = new_data$RT)
    old <- frm$df
    df <- merge(new, old) # Only keep rows where we have NAME and SMILES in both datasets
    pvec <- c("RT", "I(RT^2)", "I(RT^3)", "log(RT)", "exp(RT)", "sqrt(RT)")[predictors]
    pstr <- paste(pvec, collapse = " + ")
    fmstr <- paste("RT_ADJ ~", pstr)
    fm <- as.formula(fmstr)
    catf("Formula: %s", fmstr)
    cv <- list(
        folds = createFolds(y = df$RT, k = nfolds),
        models = vector("list", nfolds),
        preds = rep(NA, nrow(df)),
        preds_adjonly = rep(NA, nrow(df))
    )

    catf("Fitting adjustment model on full new data set")
    model <- lm(formula = fm, data = df)

    catf("Estimating performance of adjusted model in CV")
    for (i in seq_along(cv$folds)) {
        train <- unname(unlist(cv$folds[-i]))
        cv$models[[i]] <- lm(formula = fm, data = df[train, ])
        # First create adjusted predictions using the true RT from the original
        # column. This allows us to see how the adjustment model would perform
        # if we could predict RT's for the original column with 100% accuracy.
        test <- cv$folds[[i]]
        testdf <- df[test, ]
        cv$preds_adjonly[test] <- predict(cv$models[[i]], testdf)
        # Now create adjusted predictions using the RT values predicted from the
        # orginal model. This is the scenario we will encounter for future data
        # as well.
        testdf$RT <- predict(frm, testdf, adjust = FALSE, verbose = 0)
        cv$preds[test] <- predict(cv$models[[i]], testdf)
    }

    catf("Returning adjusted frm object")
    frm$adj <- list(model = model, df = df, cv = cv)
    frm
}

#' @export
#' @keywords public
#'
#' @title Predict retention times using a FastRet Model
#'
#' @description
#' Predict retention times for new data using a FastRet Model (FRM).
#'
#' @param object
#' An object of class `frm` as returned by [train_frm()].
#'
#' @param df
#' A data.frame with the same columns as the training data.
#'
#' @param adjust
#' If `object` was adjusted using [adjust_frm()], it will contain a property
#' `object$adj`. If `adjust` is TRUE, `object$adj` will be used to adjust
#' predictions obtained from `object$model`. If FALSE `object$adj` will be
#' ignored. If NULL, `object$model` will be used, if available.
#'
#' @param verbose
#' A logical value indicating whether to print progress messages.
#'
#' @param ...
#' Not used. Required to match the generic signature of `predict()`.
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
predict.frm <- function(object = train_frm(), df = object$df, adjust = NULL, verbose = 0, ...) {

    pkg <- if (inherits(object$model, "glmnet")) "glmnet" else "xgboost"
    withr::local_package(pkg)
    logf <- if (verbose == 1) catf else null
    predictors <- get_predictors(object)
    dgp <- get_dgp(predictors)
    iat <- any(grepl(":", predictors))

    if (isTRUE(adjust) && is.null(object$adj)) {
        errmsg <- "Model has not been adjusted yet. Please adjust the model first using `adjust_frm()`."
        stop(errmsg)
    }

    if (sum(predictors %in% colnames(df)) == 0 && "SMILES" %in% colnames(df)) {
        logf("Chemical descriptors not found in newdata. Trying to calculate them from the provided SMILES.")
        cds <- preprocess_data(df, dgp, iat, verbose, 1, FALSE, FALSE, TRUE)
        df <- cbind(df, cds)
    }

    if (!all(predictors %in% colnames(df))) {
        missing <- paste(setdiff(predictors, colnames(df)), collapse = ", ")
        errmsg <- paste("The following predictors are missing in `df`: ", missing)
        stop(errmsg)
    }

    logf("Predicting retention times")
    yhat <- c(predict(object$model, as.matrix(df[, predictors])))
    logf("Predictions: %s", paste(round(yhat, 2), collapse = ", "))

    if (!is.null(object$adj) && (isTRUE(adjust) || is.null(adjust))) {
        logf("Adjusting predictions using the adjustment model")
        yhatadj <- predict(object$adj$model, data.frame(RT = yhat))
        logf("Adjusted predictions: %s", paste(round(yhatadj, 2), collapse = ", "))
        logf("Returning adjusted predictions")
        return(yhatadj)
    } else {
        logf("Returning (unadjusted) predictions")
        return(yhat)
    }
}

#' @export
#' @title Extract predictor names from an 'frm' object
#' @description Extracts the predictor names from an 'frm' object.
#' @param frm An object of class 'frm' from which to extract the predictor names.
#' @return A character vector with the predictor names.
#' @keywords internal
#' @examples
#' frm <- read_rp_lasso_model_rds()
#' get_predictors(frm)
get_predictors <- function(frm = train_frm()) {
    m <- frm$model
    if (inherits(m, "glmnet")) rownames(m$beta) else m$feature_names
}


# Preprocessing #####

#' @noRd
#' @description Get RMSE, Rsquared, MAE and %below1min for a specific dataset and model.
#' @param data dataframe with retention time in the first column
#' @param model object useable as input for [stats::predict()]
get_stats <- function(df, model) {
    X <- as.matrix(df[, colnames(df) %in% CDFeatures])
    y <- df$RT
    yhat <- stats::predict(model, X, type = "response")
    measures <- c(
        RMSE = sqrt(mean((y - yhat)^2)),
        Rsquared = cor(y, yhat)^2,
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

validate_inputdata <- function(df,
                               require = c("RT", "SMILES", "NAME"),
                               min_cds = 1) {
    missing_cols <- setdiff(require, colnames(df))
    if (length(missing_cols) > 0) stop(sprintf("missing columns: %s", paste(missing_cols, collapse = ", ")))
    n_cds <- sum(colnames(df) %in% CDFeatures)
    if (n_cds < min_cds) {
        msg <- sprintf("At least %d chemical descriptors are required, but only %d are present", min_cds, n_cds)
        stop(msg)
    }
    unnown_cols <- setdiff(colnames(df), c("RT", "SMILES", "NAME", CDFeatures))
    if (length(unnown_cols) > 0) {
        msg <- sprintf("Unknown columns present: %s", paste(unnown_cols, collapse = ", "))
        stop(msg)
    }
    invisible(df)
}

validate_inputmodel <- function(model) {
    model_nams <- names(model)
    expected_names <- c("model", "df", "cv")
    n_missing <- sum(!expected_names %in% model_nams)
    if (n_missing > 0) {
        if (n_missing < length(expected_names)) {
            missing <- paste(setdiff(expected_names, model_nams), collapse = ", ")
            errmsg1 <- sprintf("Model object is missing required elements: %s.", missing)
        } else {
            errmsg1 <- sprintf("Model object is invalid.")
        }
        errmsg2 <- sprintf("Please upload a model trained with FastRet version %s or greater.", packageVersion("FastRet"))
        errmsg <- paste(errmsg1, errmsg2)
        stop(errmsg)
    }
    invisible(model)
}

# Fitting #####

fit_glmnet <- function(X, y = NULL, method = "lasso") {
    alpha <- switch(method, "lasso" = 1, "ridge" = 0)
    cvobj <- glmnet::cv.glmnet(
        X, y, alpha = alpha, standardize = TRUE, family = "gaussian",
        type.measure = "mse", grouped = FALSE
    )
    model <- glmnet::glmnet(
        X, y, alpha = alpha, standardize = TRUE, family = "gaussian",
        lambda = cvobj$lambda.min
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
#' @param xpar
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
#' gbt <- fit_gbtree(X, y, xpar = "rpopt", seed = 42)
fit_gbtree <- function(X, y, xpar = "rpopt", seed = NULL, verbose = 1,
                       nw = 1, nfolds = 10, nrounds = 2000, nthread = 1) {
    if (is.numeric(seed)) set.seed(seed)
    params <- switch(
        xpar,
        "default" = list(),
        "rpopt" = list(eta = 0.05, max_depth = 4, min_child_weight = 4, subsample = 0.5),
        find_params_best(X, y, xpar, nfolds, nw, nthread, nrounds, verbose, seed)$best_params
    )
    params$nthread <- nthread
    Xmat <- as.matrix(X)
    data <- xgboost::xgb.DMatrix(Xmat, label = y, nthread = nthread)
    xverb <- if (verbose == 2) TRUE else FALSE
    cvobj <- xgboost::xgb.cv(
        params, data, nrounds, nfold = nfolds, early_stopping_rounds = 20,
        objective = "reg:squarederror", verbose = xverb
    )
    nrounds <- cvobj$best_iteration
    xgboost::xgb.train(params, data, nrounds, verbose = xverb)
}

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
