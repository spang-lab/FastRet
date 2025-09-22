# Public #####

#' @export
#' @title Selective Measuring
#' @description
#' The function [adjust_frm()] is used to modify existing FastRet models based
#' on changes in chromatographic conditions. It requires a set of molecules with
#' measured retention times on both the original and new column. This function
#' selects a sensible subset of molecules from the original dataset for
#' re-measurement. The selection process includes:
#'
#' 1. Generating chemical descriptors from the SMILES strings of the molecules.
#'    These are the features used by [train_frm()] and [adjust_frm()].
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
#' @param k_cluster
#' The number of clusters for PAM clustering.
#' @param verbose
#' The level of verbosity.
#' @param seed
#' An optional random seed for reproducibility, set at the beginning of the
#' function.
#' @return
#' A list containing the following elements:
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
#' @keywords public
#' @examples
#' x <- selective_measuring(RP[1:50, ], k = 5, verbose = 0)
#' # For the sake of a short runtime, only the first 50 rows of the RP dataset
#' # were used in this example. In practice, you should always use the entire
#' # dataset to find the optimal subset for re-measurement.
selective_measuring <- function(raw_data,
                                k_cluster = 25,
                                verbose = 1,
                                seed = NULL) {


    if (!is.null(seed)) set.seed(seed)

    # Configure logging behaviours
    catf <- if (verbose >= 1) catf else function(...) {}
    catf("Starting Selective Measuring")

    # Return pregenerated results if mocking is enabled for this function
    if ("selective_measuring" %in% getOption("FastRet.mocks", c())) {
        catf("Mocking is enabled. Returning 'mockdata/clustering.rds'")
        return(readRDS(pkg_file("mockdata/clustering.rds")))
    }

    catf("Preprocessing input data")
    validate_inputdata(raw_data, min_cds = 0)
    df <- preprocess_data(raw_data, verbose = verbose)

    catf("Standardizing features")
    dfz <- scale(df[, -which(colnames(df) %in% c("NAME", "SMILES"))])
    dfz_noRT <- dfz[, -which(colnames(df) == "RT")]

    catf("Training Ridge Regression model")
    model <- fit_ridge(dfz, verbose = verbose)
    coefs <- glmnet::coef.glmnet(model)@x[-1] # remove intercept

    catf("Scaling features by coefficients of Ridge Regression model")
    dfzb <- data.frame(t(apply(dfz_noRT, 1, function(z) z * coefs)))
    dfzb <- `colnames<-`(dfzb, colnames(dfz_noRT))

    catf("Applying PAM clustering")
    RTcol <- which(colnames(df) == "RT")
    clobj <- cluster::pam(dfzb[, -RTcol], k = as.numeric(k_cluster))

    catf("Returning clustering results")
    CLUSTER = clobj$clustering
    IS_MEDOID = seq_len(nrow(raw_data)) %in% clobj$id.med
    ret <- list(
        clustering = data.frame(raw_data, CLUSTER, IS_MEDOID),
        clobj = clobj,
        coefs = coefs,
        model = model,
        df = df,
        dfz = dfz,
        dfzb = dfzb
    )
    ret
}
