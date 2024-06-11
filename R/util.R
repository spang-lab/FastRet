# Standard Lib Imports
#' @import graphics
#' @import grDevices
#' @import parallel
#' @import stats
#' @import utils

# External Imports
#' @import ggplot2
#' @import rcdk
#' @import shiny
#' @import shinybusy
#' @import shinyhelper
#' @import xgboost
#' @import glmnet
# Note: xgboost and glmnet must be imported because we use their S3 methods for predict. If we remove the imports, we need to use :: to access their predict explicitly.

globalVariables(".data") # to avoid warnings about NSE in ggplot2 calls (https://ggplot2.tidyverse.org/articles/ggplot2-in-packages.html)

# Colors #####

GREY <- "\033[1;30m"
RED <- "\033[1;31m"
GREEN <- "\033[1;32m"
YELLOW <- "\033[1;33m"
BLUE <- "\033[1;34m"
RESET <- "\033[0m"

# Print #####

#' @title now function
#' @description This function returns the current system time formatted according to the provided format string.
#' @param format A string representing the desired time format. Default is "%Y-%m-%d %H:%M:%OS2".
#' @return A string representing the current system time in the specified format.
#' @keywords internal
#' @export
now <- function(format = "%Y-%m-%d %H:%M:%OS2") {
    format(Sys.time(), format)
}

#' @title catf function
#' @description This function prints a formatted string with optional prefix and end strings.
#' @param ... Arguments to be passed to sprintf for string formatting.
#' @param prefix A function returning a string to be used as the prefix. Default is a timestamp.
#' @param end A string to be used as the end of the message. Default is a newline character.
#' @return No return value. This function is called for its side effect of printing a message.
#' @keywords internal
#' @export
catf <- function(..., prefix = .Options$FastRet.catf.prefix, end = .Options$FastRet.catf.end) {
    prefixstr <- if (is.null(prefix)) sprintf("%s%s%s ", GREY, now(), RESET) else prefix()
    endstr <- if (is.null(end)) "\n" else end
    middlestr <- sprintf(...)
    msg <- sprintf("%s%s%s", prefixstr, middlestr, endstr)
    cat(msg)
}

# Multicore #####

#' @description Calculate the number of workers for parallel processing
#' @param mult A multiplier for the number of cores. Default is 0.5.
#' @param nmax The maximum number of workers allowed. Default is 16.
#' @return The number of workers to be used for parallel processing.
#' @examples
#' get_n_workers()
#' get_n_workers(2)
#' get_n_workers(0.5, 10)
#' @noRd
get_n_workers <- function(mult = 0.5, nmax = 16) {
    n <- parallel::detectCores()
    nmul <- ceiling(n * mult)
    min(nmul, nmax)
}

# Datasets #####

#' @title Retention Times (RT) Measured on a Reverse Phase (RP) Column
#' @description Retention time data from a reverse phase liquid chromatography measured with a temperature of 35 degree and a flowrate of 0.3ml/min. The same data is available as an xlsx file in the package. To read it into R use `read_rp_xlsx()`.
#' @format A dataframe of 442 metabolites with the following columns:
#' \describe{
#'   \item{RT}{Retention time}
#'   \item{SMILES}{SMILES notation of the metabolite}
#'   \item{NAME}{Name of the metabolite}
#' }
#' @source Measured by functional genomics lab at the University of Regensburg.
#' @keywords dataset
#' @seealso read_rp_xlsx
"RP"

update_RP <- function() {
    RP <- read_rp_xlsx()
    usethis::use_data(RP, overwrite = TRUE)
}

#' @title Read retention times (RT) measured on a reverse phase (RP) column
#' @description Read retention time data from a reverse phase liquid chromatography measured with a temperature of 35 degree and a flowrate of 0.3ml/min. The data also exists as dataframe in the package. To use it directly in R just enter `RP`.
#' @return A dataframe of 442 metabolites with columns `RT`, `SMILES` and `NAME`.
#' @keywords dataset
#' @source Measured by functional genomics lab at the University of Regensburg.
#' @seealso RP
#' @export
read_rp_xlsx <- function() {
    xlsx::read.xlsx(pkg_file("extdata/RP.xlsx"), 1)
}

#' @title Hypothetical retention times (RT) measured on a reverse phase (RP) column
#' @description Subset of the data from [read_rp_xlsx()] with some slight modifications to simulate changes in temperature and/or flowrate.
#' @format A dataframe of 25 metabolites and columns `RT`, `SMILES` and `NAME`.
#' @keywords dataset
#' @export
read_rpadj_xlsx <- function() {
    xlsx::read.xlsx(pkg_file("extdata/RP_adj.xlsx"), 1)
}

# Caching #####

#' @title Get cache directory
#' @description Creates and returns the cache directory for the FastRet package.
#' @param subdir Optional subdirectory within the cache directory.
#' @return The path to the cache directory or subdirectory.
#' @keywords internal
#' @export
get_cache_dir <- function(subdir = NULL) {
    cache_dir <- tools::R_user_dir("FastRet", which = "cache")
    if (!is.null(subdir)) cache_dir <- file.path(cache_dir, subdir)
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    normalizePath(cache_dir, winslash = "/", mustWork = FALSE)
}

# CONTINUE HERE:
# Error: object 'CDNames' not found

#' @title Get package file
#' @description Returns the path to a file within the FastRet package.
#' @param path The path to the file within the package.
#' @param mustWork If TRUE, an error is thrown if the file does not exist.
#' @return The path to the file.
#' @keywords internal
#' @export
pkg_file <- function(path, mustWork = FALSE) {
    system.file(path, package = "FastRet", mustWork = mustWork)
}

# Misc #####

#' @title Collect elements from a list of lists
#' @description This function takes a list of lists where each inner list has the same names. It returns a list where each element corresponds to a name of the inner list that is extracted from each inner list. Especially useful for collecting results from lapply.
#' @param xx A list of lists where each inner list has the same names.
#' @return A list where each element corresponds to a name of the inner list that is extracted from each inner list.
#' @examples
#' \dontrun{
#' xx <- lapply(1:3, function(i) list(a = i, b = i^2, c = i^3))
#' ret <- collect(xx)
#' }
#' @keywords internal
#' @export
collect <- function(xx) {
    ns <- names(xx[[1]])
    ret <- lapply(ns, function(n) sapply(xx, function(x) x[[n]]))
    `names<-`(ret, ns)
}

# Test Helpers #####

serve_docs <- function() {
    servr::httd("docs")
}

make_docs <- function(reload = TRUE) {
    if (reload) {
        pkgdown::build_site(lazy = TRUE)
        vignette_files <- dir(pkg_file("vignettes"), recursive = TRUE, full.names = TRUE)
        timestamps_old <- file.mtime(vignette_files)
        while (TRUE) {
            Sys.sleep(0.1)
            timestamps_new <- file.mtime(vignette_files)
            if (!identical(timestamps_old, timestamps_new)) {
                pkgdown::build_site(lazy = TRUE)
                timestamps_old <- timestamps_new
            }
        }
    } else {
        pkgdown::build_site(lazy = TRUE)
    }
}
update_mockdata <- function(getCD = FALSE,
                            preprocess_data = FALSE,
                            train_frm = FALSE,
                            selective_measuring = FALSE,
                            all = FALSE) {
    options(FastRet.mocks = c()) # reset all mocks
    if (getCD || all) {
        getCDs <- TRUE
        preprocess_data <- TRUE
        train_frm <- TRUE
        cds <- getCDs(verbose = 1, nw = 4)
        save_mockdata(cds, "RPCD")
    }
    if (preprocess_data || all) {
        options(FastRet.mocks = c("getCDs")) # mock getCDs to speed up below calls
        df <- preprocess_data()
        save_mockdata(df, "RPCD_prepro")
    }
    if (train_frm || all) {
        options(FastRet.mocks = c("preprocess_data")) # mock preprocess_data to speed up below calls
        df <- preprocess_data()
        ridge_model <- train_frm(df, method = "ridge", nw = 5)
        lasso_model <- train_frm(df, method = "lasso", nw = 5)
        gbtree_model <- train_frm(df, method = "gbtree", nw = 5)
        save_mockdata(ridge_model, "ridge_model")
        save_mockdata(lasso_model, "lasso_model")
        save_mockdata(gbtree_model, "gbtree_model")
    }
    if (selective_measuring || all) {
        options(FastRet.mocks = c("preprocess_data")) # mock preprocess_data to speed up below calls
        obj <- selective_measuring(raw_data = read_rp_xlsx(), k_cluster = 25)
        save_mockdata(obj, "clustering")
    }
}

save_mockdata <- function(obj, name, xlsx = FALSE) {
    mockdatapath <- system.file("mockdata", package = "FastRet")
    rdspath <- file.path(mockdatapath, paste0(name, ".rds"))
    xlsxpath <- file.path(mockdatapath, paste0(name, ".xlsx"))
    callstr <- deparse(substitute(obj))
    catf("Evaluating %s", callstr)
    obj <- force(obj)
    catf("Saving %s", rdspath)
    saveRDS(obj, rdspath)
    if (xlsx) {
        catf("Saving %s", xlsxpath)
        xlsx::write.xlsx(obj, xlsxpath, row.names = FALSE)
    }
}

#' @description To enable either call add "eval(.Options$FastRet.onFuncEntry)" at the beginning of each function or call `trace(f, quote(eval(.Options$FastRet.onFuncEntry)))` for each function `f` you want to trace.
#' @noRd
enable_function_tracing <- function() {
    onFuncEntry <- quote({
        funcname <- as.character(sys.call(-2))
        cat(sprintf("Start: %s\n", funcname))
        on.exit(cat(sprintf("Exit: %s\n", funcname)), add = TRUE)
    })
    options(FastRet.onFuncEntry = onFuncEntry)
}
