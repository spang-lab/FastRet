#' @title Selective Measuring Function
#' @description This function performs selective measuring on the provided raw data. It includes preprocessing, standardizing features, training a Ridge Regression model, scaling features by coefficients of the Ridge Regression model, applying PAM clustering, and returning clustering results.
#' @param raw_data The raw data to be processed. Default is the result of `read_rp_xlsx()`.
#' @param k_cluster The number of clusters for PAM clustering. Default is 25.
#' @param verbose The level of verbosity. Default is 1.
#' @return A list containing the following elements:
#' * clustering: a data frame with raw data, cluster assignments, and medoid indicators
#' * clobj: the PAM clustering object
#' * coefs: the coefficients from the Ridge Regression model
#' * model: the Ridge Regression model
#' * df: the preprocessed data
#' * dfz: the standardized features
#' * dfzb: the features scaled by coefficients of the Ridge Regression model
#' @keywords public
#' @export
selective_measuring <- function(raw_data = read_rp_xlsx(), k_cluster = 25, verbose = 1) {

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
    df <- preprocess_data(raw_data)

    catf("Standardizing features")
    dfz <- scale(df[, -which(colnames(df) %in% c("NAME", "SMILES"))])
    dfz_noRT <- dfz[, -which(colnames(df) == "RT")]

    catf("Training Ridge Regression model")
    model <- fit_ridge(dfz)
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
