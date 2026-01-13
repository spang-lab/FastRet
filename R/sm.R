# Public #####

#' @export
#' @keywords public
#'
#' @title Selective Measuring
#'
#' @description
#' The function [adjust_frm()] is used to modify existing FastRet models based
#' on changes in chromatographic conditions. It requires a set of molecules with
#' measured retention times on both the original and new column. This function
#' selects a sensible subset of molecules from the original dataset for
#' re-measurement. The selection process includes:
#'
#' 1. Generating chemical descriptors from the SMILES strings of the molecules.
#'    These are the features used by [train_frm()] and [adjust_frm()].˝
#' 2. Standardizing chemical descriptors to have zero mean and unit variance.
#' 3. Training a Ridge Regression model with the standardized chemical
#'    descriptors as features and the retention times as the target variable.
#' 4. Scaling the chemical descriptors by coefficients of the Ridge Regression
#'    model.
#' 5. Clustering the entire dataset, which includes the scaled chemical
#'    descriptors and the retention times.
#' 6. Returning the clustering results, which include the cluster assignments,
#'    the medoid indicators, and the raw data.
#'
#' @param raw_data
#' The raw data to be processed.
#' Must be a dataframe with columns NAME, RT and SMILES.
#'
#' @param k_cluster
#' The number of clusters for PAM clustering.
#'
#' @param verbose
#' The level of verbosity.
#'
#' @param seed
#' An optional random seed for reproducibility, set at the beginning of the
#' function.
#'
#' @param rt_coef
#' Which coefficient to use for scaling RT before clustering. Options are:
#' - `0`: exclude RT from the clustering.
#' - 'max' / 'max_ridge_coef': scale with the maximum absolute coefficient
#'   obtained in ridge regression. I.e., RT will have approximately the same
#'   weight as the most important chemical descriptor.
#' - `1`: do not scale RT any further, i.e., use standardized RT. The effect of
#'   leaving RT unscaled is kind of unpredictable, as the ridge coefficients
#'   depend on the dataset. If the maximum absolute coefficient is much smaller
#'   than 1, RT will dominate the clustering. If it is much larger than 1, RT
#'   will have little influence on the clustering.
#' - 'inf': set all chemical descriptor values to zero, i.e., RT is "infinitely"
#'   more important than any chemical descriptor.
#'
#' @return
#' A list containing the following elements:
#'
#' - `clustering`: A data frame with columns RT, SMILES, NAME, CLUSTER and
#'   IS_MEDOID.
#' - `clobj`: The clustering object. The object returned by the clustering
#'   function. Depends on the `method` parameter.
#' - `coefs`: The coefficients from the Ridge Regression model.
#' - `model`: The Ridge Regression model.
#' - `df`: The preprocessed data.
#' - `dfz`: The standardized features.
#' - `dfzb`: The features scaled by the coefficients (betas) of the Ridge
#'   Regression model.
#'
#' @examples
#' x <- selective_measuring(RP[1:50, ], k = 5, verbose = 0)
#' # For the sake of a short runtime, only the first 50 rows of the RP dataset
#' # were used in this example. In practice, you should always use the entire
#' # dataset to find the optimal subset for re-measurement.
#'
selective_measuring <- function(raw_data,
                                k_cluster = 25,
                                verbose = 1,
                                seed = NULL,
                                rt_coef = "max_ridge_coef"
                                ) {

    stopifnot(
        is.data.frame(raw_data),
        all(c("NAME", "SMILES", "RT") %in% colnames(raw_data)),
        is.numeric(k_cluster), length(k_cluster) == 1, !is.na(k_cluster), k_cluster >= 2,
        is.logical(verbose) || is.numeric(verbose),
        is.null(seed) || is.numeric(seed),
        is.numeric(rt_coef) || is.character(rt_coef), length(rt_coef) == 1, !is.na(rt_coef)
    )
    rt_coef <- try_as_numeric(rt_coef, fallback = rt_coef)
    logf <- if (verbose >= 1) catf else null

    logf("Starting Selective Measuring")
    if (is.numeric(seed)) set.seed(seed)

    logf("Preprocessing input data")
    df <- preprocess_data(raw_data, verbose = verbose)
    # Now df contains only NAME, SMILES, RT [,INCHIKEY] and CDs

    logf("Standardizing features")
    nonmeta <- !colnames(df) %in% c("NAME", "SMILES", "INCHIKEY")
    df <- df[, nonmeta, drop = FALSE]
    dfz <- as.data.frame(scale(df)) # z-score standardized (mean 0, sd 1)

    if (rt_coef %in% c(Inf, "inf")) {
        logf("Setting CDs to zero because rt_coef is Inf")
        coefs <- rep(0, ncol(dfz) - 1)
        names(coefs) <- colnames(dfz)[colnames(dfz) != "RT"]
        model <- NULL
        dfzb <- data.frame(RT = dfz$RT)
    } else {
        logf("Training Ridge Regression model")
        Xz <- as.matrix(dfz[, colnames(dfz) != "RT", drop = FALSE])
        y <- as.numeric(dfz[, "RT"])
        model <- fit_glmnet(Xz, y, method = "ridge", seed = seed)

        logf("Scaling features by coefficients of Ridge Regression model")
        coef_mat <- glmnet::coef.glmnet(model) # (m+1) x 1 matrix
        coefs <- as.numeric(coef_mat)[-1] # remove intercept
        coefs <- setNames(coefs, rownames(coef_mat)[-1])
        coefs <- coefs[colnames(Xz)] # ensure correct order
        Xzb <- sweep(Xz, 2, coefs, `*`) # z-score standardized and beta-scaled

        logf("Scaling RT by %s before clustering", rt_coef)
        rtc <- if (grepl("max", rt_coef)) max(abs(coefs), na.rm = TRUE)
            else if (is.numeric(rt_coef)) as.numeric(rt_coef)
            else stop("Invalid value for rt_coef: ", rt_coef)
        dfzb <- data.frame(RT = dfz$RT * rtc, Xzb)
    }

    logf("Applying PAM clustering")
    clobj <- cluster::pam(dfzb, k = as.numeric(k_cluster))

    logf("Returning clustering results")
    CLUSTER = clobj$clustering
    IS_MEDOID = seq_len(nrow(raw_data)) %in% clobj$id.med
    ret <- list(
        clustering = data.frame(raw_data, CLUSTER, IS_MEDOID),
        clobj = clobj, coefs = coefs, model = model,
        df = df, dfz = dfz, dfzb = dfzb
    )
    ret
}

try_as_numeric <- function(x, fallback = x) {
    if (length(x) != 1) stop("scalar expected")
    y <- suppressWarnings(as.numeric(x))
    if (!is.na(y)) y else fallback
}