# Read From Disk (Public) #####

#' @export
#' @title Read retention times (RT) measured on a reverse phase (RP) column
#' @description
#' Reads retention times from a reverse phase liquid chromatography experiment,
#' performed at 35\eqn{^\circ}C and a flow rate of 0.3 mL/min. The data is also available
#' as a dataframe in the package; to access it directly, use [RP].
#' @return
#' A dataframe of 442 metabolites with columns `RT`, `SMILES` and `NAME`.
#' @keywords dataset
#' @source
#' Measured by the Institute of Functional Genomics at the University of
#' Regensburg.
#' @seealso RP
#' @examples
#' x <- read_rp_xlsx()
#' all.equal(x, RP)
read_rp_xlsx <- function() {
    openxlsx::read.xlsx(pkg_file("extdata/RP.xlsx"), 1)
}

#' @export
#' @title Read the HILIC dataset from the Retip package
#' @description
#' If the [Retip](https://www.retip.app/) package is installed, its HILIC
#' dataset is read from its installation directory. If not, the dataset is
#' downloaded from
#' `https://github.com/oloBion/Retip/raw/master/data/HILIC.RData` to a temporary
#' file and then read from there.
#' @param verbose Verbosity level. Set to 1 to enable progress messages,
#' 0 to suppress them.
#' @return A data frame containing the Retip HILIC dataset.
#' @examples
#' df <- read_retip_hilic_data(verbose = 0)
#' @references
#' Retip: Retention Time Prediction for Compound Annotation in Untargeted
#' Metabolomics Paolo Bonini, Tobias Kind, Hiroshi Tsugawa, Dinesh Kumar
#' Barupal, and Oliver Fiehn Analytical Chemistry 2020 92 (11), 7515-7522 DOI:
#' 10.1021/acs.analchem.9b05765
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

#' @export
#' @title Hypothetical retention times
#' @description
#' Subset of the data from [read_rp_xlsx()] with some slight modifications to
#' simulate changes in temperature and/or flowrate.
#' @return
#' A dataframe with 25 rows (metabolites) and 3 columns: RT, SMILES and NAME.
#' @examples \donttest{
#' x <- read_rpadj_xlsx()
#' }
#' @keywords dataset
read_rpadj_xlsx <- function() {
    openxlsx::read.xlsx(pkg_file("extdata/RP_adj.xlsx"), 1)
}

#' @export
#' @title LASSO Model trained on RP dataset
#' @description
#' Read a LASSO model trained on the [RP] dataset using [train_frm()].
#' @return A `frm` object.
#' @keywords dataset
#' @examples
#' frm <- read_rp_lasso_model_rds()
read_rp_lasso_model_rds <- function() {
    readRDS(pkg_file("extdata/RP_lasso_model.rds"))
}


# Lazyload Objects (Public) #####

## RP #####

#' @title Retention Times (RT) measured on a Reverse Phase (RP) Column
#' @description
#' Retention time data from a reverse phase liquid chromatography measured with
#' a temperature of 35\eqn{^\circ}C and a flowrate of 0.3ml/min. The same data
#' is available as an xlsx file in the package. To read it into R use
#' `read_rp_xlsx()`. @format A dataframe of 442 metabolites with the following
#' columns:
#' \describe{
#'   \item{RT}{Retention time}
#'   \item{SMILES}{SMILES notation of the metabolite}
#'   \item{NAME}{Name of the metabolite}
#' }
#' @source
#' Measured by the Institute of Functional Genomics at the University of
#' Regensburg.
#' @keywords dataset
#' @seealso read_rp_xlsx
"RP"

# Make Lazyload Objects (Private) #####

make_lazyload_objs <- function(HILIC, RP, RP_Mod, RP_Val) {
    # Read in data from Excel files
    RP <- read_rp_xlsx()
    # Store objects in the data folder
    usethis::use_data(RP, overwrite = TRUE)
}
