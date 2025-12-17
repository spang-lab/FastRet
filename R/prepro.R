# Public #####

#' @export
#' @keywords public
#'
#' @title Preprocess data
#'
#' @description
#' Preprocess data so they can be used as input for [train_frm()].
#'
#' @param data
#' Dataframe with following columns:
#' - Mandatory: NAME, RT and SMILES.
#' - Recommmended: INCHIKEY.
#' - Optional: Any of the chemical descriptors listed in [CDFeatures].
#' All other columns will be removed.
#' See 'Details'.
#' @param degree_polynomial
#' Add predictors with polynomial terms up to the specified degree, e.g. 2 means
#' "add squares", 3 means "add squares and cubes". Set to 1 to leave descriptors
#' unchanged.
#' @param interaction_terms
#' Add interaction terms? Polynomial terms are not included in the generation of
#' interaction terms.
#' @param verbose 0: no output, 1: show progress, 2: progress and warnings.
#' @param nw Number of workers to use for parallel processing.
#' @param rm_near_zero_var Remove near zero variance predictors?
#' @param rm_na Remove NA values?
#' @param add_cds Add chemical descriptors using [getCDs()]? See 'Details'.
#' @param rm_ucs Remove unsupported columns?
#' @param rt_terms
#' Which retention-time transformations to append as extra predictors. Supply a
#' numeric vector referencing predefined rt_terms (1=RT, 2=I(RT^2),
#' 3=I(RT^3), 4=log(RT), 5=exp(RT), 6=sqrt(RT)) or a character vector with the
#' explicit transformation terms. Character values are passed to [model.frame()],
#' so they must use valid formula syntax (e.g. "I(RT^2)" rather than "RT^2").
#' @param mandatory
#' Character vector of mandatory columns that must be present in `data`. If any
#' of these columns are missing, an error is raised.
#'
#' @details
#' If `add_cds = TRUE`, chemical descriptors are added using [getCDs()]. If
#' **all** chemical descriptors listed in [CDFeatures] are already present in
#' the input `data` object, [getCDs()] will leave them unchanged. If one or more
#' chemical descriptors are missing, **all** chemical descriptors will be
#' recalculated and existing ones will be overwritten.
#'
#' @return
#' A dataframe with the preprocessed data.
#'
#' @examples
#' data <- head(RP, 3)
#' pre <- preprocess_data(data, verbose = 0)
preprocess_data <- function(data,
                            degree_polynomial = 1,
                            interaction_terms = FALSE,
                            verbose = 1,
                            nw = 1,
                            rm_near_zero_var = TRUE,
                            rm_na = TRUE,
                            add_cds = TRUE,
                            rm_ucs = TRUE,
                            rt_terms = 1,
                            mandatory = c("NAME", "RT", "SMILES")
                            ) {

    if (FALSE) stub(preprocess_data, data=head(RP, 3), rt_terms=c(1:2,6))

    logf <- if (verbose >= 1) catf else null
    dbgf <- if (verbose >= 2) catf else null

    logf("Preprocessing dataframe with dimension %d x %d", nrow(data), ncol(data))

    # Convert rt_term numbers to names
    all_terms <- c("RT", "I(RT^2)", "I(RT^3)", "log(RT)", "exp(RT)", "sqrt(RT)")
    if (is.numeric(rt_terms)) {
        if (!all(rt_terms %in% 1:6)) {
            stop("Invalid numbers in 'rt_terms'. Valid: 1 to 6.")
        }
        rt_terms <- all_terms[rt_terms]
    } else if (is.character(rt_terms)) {
        if (!all(rt_terms %in% all_terms)) {
            stop("Invalid names in 'rt_terms'. Valid: ", as_str(all_terms))
        }
    } else {
        stop("'rt_terms' must be either numeric or character vector.")
    }
    rt_terms <- rt_terms[rt_terms != "RT"] # RT is 'mandatory' anyways

    # Check for missing columns
    optional <- c("INCHIKEY", "RT_ADJ")
    supported <- c(mandatory, optional, rt_terms, CDFeatures)
    is_supported <- colnames(data) %in% supported
    if (!all(mandatory %in% colnames(data))) {
        missing <- mandatory[!mandatory %in% colnames(data)]
        stop("Missing mandatory columns: ", as_str(missing))
    }

    # Remove unsupported columns if requested and necessary
    if (!all(is_supported) && rm_ucs) {
        logf("Removing %d unsupported columns", sum(!is_supported))
        data <- data[, is_supported, drop = FALSE]
    }

    # Add chemical descriptors if requested
    if (add_cds) {
        logf("Obtaining chemical descriptors using %d workers", nw)
        df <- getCDs(data, verbose, nw)
    } else {
        df <- data
    }

    # Split data into metadata and chemical descriptors
    iscd <- colnames(df) %in% CDFeatures
    M <- df[, !iscd, drop = FALSE] # Metadata (NAME, RT, SMILES, INCHIKEY, )
    X <- df[, iscd, drop = FALSE] # Chemical Descriptors

    # Backup original X for use in upcoming transformations
    Xorig <- X

    # Add polynomial terms if requested
    if (degree_polynomial >= 2 && ncol(X) > 0) {
        logf("Adding polynomial predictors up to degree %d", degree_polynomial)
        for (p in seq(2, degree_polynomial)) {
            P <- as.data.frame(lapply(Xorig, function(col) col^p))
            colnames(P) <- paste0(colnames(Xorig), "^", p)
            X <- cbind(X, P)
        }
    }

    # Add interaction terms if requested
    np <- ncol(Xorig)
    if (interaction_terms && ncol(X) > 0 && np >= 2) {
        logf("Adding interaction terms")
        withr::local_options(list(na.action='na.pass'))
        I <- stats::model.matrix(~ .^2 - . - 1, data = Xorig, na.action = stats::na.pass)
        X <- cbind(X, I)
    }

    # Add transformed RT predictors if requested
    if (length(rt_terms) >= 1 && !all(rt_terms %in% colnames(M))) {
        fm <- paste0("~", paste(rt_terms, collapse = " + "))
        T <- model.frame(fm, data = M)
        M[, colnames(T)] <- T
    }

    # Remove columnbs with NAs if requested
    if (rm_na && ncol(X) > 0) {
        logf("Removing CDs with NAs")
        keep <- !apply(X, 2, function(col) any(is.na(col)))
        X <- X[, keep, drop = FALSE]
    }

    # Remove columns with near zero variance if requested
    if (rm_near_zero_var && ncol(X) > 0) {
        logf("Removing CDs with variance close to zero")
        idx <- nearZeroVar(X)
        if (length(idx)) X <- X[, -idx, drop = FALSE]
    }

    # Combine metadata and processed chemical descriptors
    df <- cbind(M, X)
    logf("Preprocessing finished")
    df
}

# Private #####

#' @noRd
#' @keywords internal
#'
#' @title Checks which chemical descriptors are suitable for linear models
#'
#' @description
#' Checks which chemical descriptors are suitable for use in linear model.
#' Chemical descriptors with missing values, near-zero variance or strong
#' outlier values are considered as not suitable.
#' @param df
#' Input data for performing the analysis. Must be a data frame with columns
#' NAME, RT and SMILES.
#'
#' @param verbose
#' A logical value indicating whether to print verbose output.
#'
#' @param nw
#' The number of workers to use for parallel processing.
#'
#' @return
#' A data frame with the predictors and their suitability status.
#'
#' @seealso [plot_lm_suitability()]
#'
#' @examples
#' x <- check_lm_suitability(head(RP, 3), verbose = FALSE, nw = 1)
#'
check_lm_suitability <- function(df = read_retip_hilic_data(),
                                 verbose = FALSE,
                                 nw = 2) {
    y <- df$RT
    cds <- getCDs(df, verbose = verbose, nw = nw)
    X <- cds[5:ncol(cds)]
    predictors <- colnames(X)
    n <- ncol(X)
    hasNAs <- apply(X, 2, function(x) any(is.na(x)))
    isAlmostConstant <- seq_len(ncol(X)) %in% nearZeroVar(X)
    hasOutliers <- apply(X, 2, function(x) any(abs(x - median(x)) > 50 * mad(x)))
    isSuitable <- !(hasNAs | isAlmostConstant | hasOutliers)
    V <- data.frame(predictors, hasNAs, isAlmostConstant, hasOutliers, isSuitable)
    list(df = df, X = X, V = V)
}

#' @noRd
#' @keywords internal
#'
#' @title Plot the suitability of predictors for linear models
#'
#' @description
#' Creates one pdf page for every predictor inside `slist$X`. The pdf page
#' consists of the following three plots shown next to each other:
#'
#' 1. Histogram
#' 2. Density plot
#' 3. Scatterplot against `slist$df$RT`
#'
#' The name of the predictor, its suitability, and the status of the checks for
#' missing values, near-zero variance, and outliers are shown in the title of
#' each plot.
#'
#' @param slist
#' A list containing the data frame `df`, the matrix `X`, and the data frame `V`
#' from [check_lm_suitability()].
#'
#' @param pdfpath
#' The path to the pdf file to save the plots.
#'
#' @param descs
#' Index of chemical descriptors to plot. Leave at NULL to plot all chemical
#' descriptors.
#'
#' @return
#' No return value. The function is used for its side effect of creating a pdf
#' file with the plots.
#'
#' @seealso [check_lm_suitability()]
#'
#' @examples
#' df <- head(RP, 3)
#' slist <- check_lm_suitability(df, verbose = FALSE, nw = 1)
#' plot_lm_suitability(slist, descs = 1:5)
#'
plot_lm_suitability <- function(slist = check_lm_suitability(),
                                pdfpath = tempfile("lm_suitability", fileext = ".pdf"),
                                descs = NULL) {
    catf("Plotting suitability of predictors for linear models to file '%s'", pdfpath)
    pdf(pdfpath, width = 9, height = 3) # A4 size in inches
    on.exit(dev.off(), add = TRUE)
    opar <- par(mfrow = c(1, 3), oma = c(2, 0, 2, 0), mar = c(1, 4, 1, 2))
    on.exit(par(opar), add = TRUE, after = FALSE)
    if (is.null(descs)) descs <- seq_len(ncol(slist$X))
    V <- slist$V[descs, ]
    X <- slist$X[, descs]
    RT <- slist$df$RT
    for (i in seq_len(ncol(X))) {
        x <- X[, i]
        name <- colnames(X)[i]
        state <- ifelse(V$isSuitable[i], "Ok", "Bad")
        reasons <- c("Has NAs", "Almost Constant", "Has Outliers")
        reasons <- reasons[unlist(V[i, c("hasNAs", "isAlmostConstant", "hasOutliers")])]
        reasons <- paste(reasons, collapse = ", ")
        if (reasons != "") reasons <- sprintf(" (%s)", reasons)
        title <- sprintf("%s: %s%s", name, state, reasons)
        tryCatch(
            {
                hist(x, main = "", xlab = "", xaxt = "n") # removed x-axis
                plot(density(x), main = "", xlab = "", xaxt = "n") # removed x-axis
                plot(x, RT, main = "", xlab = "", ylab = "RT", xaxt = "n") # removed x-axis
                axis(
                    side = 1,
                    at = seq(min(x), max(x), length.out = 5),
                    labels = seq(min(x), max(x), length.out = 5),
                    outer = TRUE
                ) # added single x-axis at the bottom
                title(main = title, outer = TRUE)
                catf("[%d/%d] %s (PLOTTED SUCCESSFULLY)", i, ncol(X), name)
            },
            error = function(e) {
                catf("[%d/%d] %s (FAILED)", i, ncol(X), name)
            }
        )
    }
}
