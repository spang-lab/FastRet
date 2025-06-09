# Lazy Loaded Objects #####

## RP #####

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

# Create Lazy Loaded Objects #####

.update_RP <- function() {
    RP <- read_rp_xlsx()
    usethis::use_data(RP, overwrite = TRUE)
}

# Read Full Training Datasets #####

#' @export
#' @title Read retention times (RT) measured on a reverse phase (RP) column
#' @description Read retention time data from a reverse phase liquid chromatography measured with a temperature of 35 degree and a flowrate of 0.3ml/min. The data also exists as dataframe in the package. To use it directly in R just enter `RP`.
#' @return A dataframe of 442 metabolites with columns `RT`, `SMILES` and `NAME`.
#' @keywords dataset
#' @source Measured by functional genomics lab at the University of Regensburg.
#' @seealso RP
#' @examples
#' x <- read_rp_xlsx()
#' all.equal(x, RP)
read_rp_xlsx <- function() {
    xlsx::read.xlsx(pkg_file("extdata/RP.xlsx"), 1)
}

#' @export
#' @title Read the HILIC dataset from the Retip package
#' @description If the [Retip](https://www.retip.app/) package is installed, its HILIC dataset is read from its installation directory. If not, the dataset is downloaded from `https://github.com/oloBion/Retip/raw/master/data/HILIC.RData` to a temporary file and then read from there.
#' @param verbose Verbosity level. 1 == print progress messages, 0 == no progress messages.
#' @return df A data frame containing the HILIC dataset.
#' @examples
#' df <- read_retip_hilic_data(verbose = 0)
#' @references
#' Retip: Retention Time Prediction for Compound Annotation in Untargeted Metabolomics Paolo Bonini, Tobias Kind, Hiroshi Tsugawa, Dinesh Kumar Barupal, and Oliver Fiehn Analytical Chemistry 2020 92 (11), 7515-7522 DOI: 10.1021/acs.analchem.9b05765
#' @keywords dataset
read_retip_hilic_data <- function(verbose = 1) {
    if (eval(parse(text = "requireNamespace('Retip', quietly = TRUE)"))) {
        return(eval(parse(text = "Retip::HILIC")))
        # We do this absurd eval trick to avoid R CMD check
        # complaining about Retip not being a suggested package.
        # (We would suggest it, but it is not on CRAN or any other
        # CRAN-like repository, so we cannot suggest it without
        # triggering further warnings.)
    }
    url <- "https://github.com/oloBion/Retip/raw/master/data/HILIC.RData"
    destfile <- tempfile("HILIC", fileext = ".RData")
    download.file(url, destfile, mode = "wb", quiet = !verbose)
    HILIC <- NULL # will be loaded in the next line
    load(destfile)
    df <- HILIC
}

.read_hilic_xlsx <- function() {
    file <- "misc/datasets/20210702_RT_Prediction_Hilic_Library_FF.xlsx"
    df <- xlsx::read.xlsx(pkg_file(file), 1)
}

.read_rp_axmm_xlsx <- function() {
    stop("Excel file with measurements from RP-AXMM column is missing")
}

# Read Smaller Datasets for Model Adjustment #####

.read_rp_steep_xlsx <- function() {
    file <- "misc/datasets/R8P_Medoids_input_Steeper gradient.xlsx"
    xlsx::read.xlsx(pkg_file(file), 1)
}

.read_rp_flat_xlsx <- function() {
    file <- "misc/datasets/R8P_Medoids_input_Flatter gradient.xlsx"
    xlsx::read.xlsx(pkg_file(file), 1)
}

.read_rp_t25_xlsx <- function() {
    file <- "misc/datasets/R8P_Medoids_input_Lower Temperature 25.xlsx"
    xlsx::read.xlsx(pkg_file(file), 1)
}

.read_rp_fr25_xlsx <- function() {
    file <- "misc/datasets/R8P_Medoids_input_Lower flow rate 025.xlsx"
    xlsx::read.xlsx(pkg_file(file), 1)
}

.read_rp_t25_fr25_xlsx <- function() {
    file <- "misc/datasets/R8P_Medoids_input_Lower Temperature 25 + Lower Flow rate 025.xlsx"
    xlsx::read.xlsx(pkg_file(file), 1)
}

.read_rp_t25_fr25_steep_xlsx <- function() {
    file <- "misc/datasets/R8P_Medoids_input_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx"
    xlsx::read.xlsx(pkg_file(file), 1)
}

# Read Validation Data Sets #####

.read_rp_steep_val_xlsx <- function() {
    file <- "misc/datasets/ValidSet_Steeper gradient.xlsx"
    xlsx::read.xlsx(pkg_file(file), 1)
}

.read_rp_flat_val_xlsx <- function() {
    file <- "misc/datasets/ValidSet_Flatter gradient.xlsx"
    xlsx::read.xlsx(pkg_file(file), 1)
}

.read_rp_t25_val_xlsx <- function() {
    file <- "misc/datasets/ValidSet_Lower Temperature 25.xlsx"
    xlsx::read.xlsx(pkg_file(file), 1)
}

.read_rp_fr25_val_xlsx <- function() {
    file <- "misc/datasets/ValidSet_Lower flow rate 025.xlsx"
    xlsx::read.xlsx(pkg_file(file), 1)
}

.read_rp_t25_fr25_val_xlsx <- function() {
    file <- "misc/datasets/ValidSet_Lower Temperature 25 + Lower Flow rate 025.xlsx"
    xlsx::read.xlsx(pkg_file(file), 1)
}

.read_rp_t25_fr25_steep_val_xlsx <- function() {
    file <- "misc/datasets/ValidSet_Steeper gradient+lower Temperature 25+Lower Flow rate 025.xlsx"
    xlsx::read.xlsx(pkg_file(file), 1)
}

# Read Models and Datasets for Package Testing #####

#' @export
#' @title Hypothetical retention times
#' @description Subset of the data from [read_rp_xlsx()] with some slight modifications to simulate changes in temperature and/or flowrate.
#' @return A dataframe with 25 rows (metabolites) and 3 columns `RT`, `SMILES` and `NAME`.
#' @examples \donttest{
#' x <- read_rpadj_xlsx()
#' }
#' @keywords dataset
read_rpadj_xlsx <- function() {
    xlsx::read.xlsx(pkg_file("extdata/RP_adj.xlsx"), 1)
}

#' @export
#' @title LASSO Model trained on RP dataset
#' @description Read a LASSO model trained on the [RP] dataset using [train_frm()].
#' @return A `frm` object.
#' @keywords dataset
#' @examples
#' frm <- read_rp_lasso_model_rds()
read_rp_lasso_model_rds <- function() {
    readRDS(pkg_file("extdata/RP_lasso_model.rds"))
}
