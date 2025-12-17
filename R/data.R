# Public #####

#' @export
#' @keywords dataset
#'
#' @title Read retention times (RT) measured on a reverse phase (RP) column
#'
#' @description
#' Reads retention times from a reverse phase liquid chromatography experiment,
#' performed at 35\eqn{^\circ}C and a flow rate of 0.3 mL/min. The data is also available
#' as a dataframe in the package; to access it directly, use [RP].
#'
#' @return
#' A dataframe of 442 metabolites with columns `RT`, `SMILES` and `NAME`.
#'
#' @source
#' Measured by the Institute of Functional Genomics at the University of
#' Regensburg.
#'
#' @seealso RP
#'
#' @examples
#' x <- read_rp_xlsx()
#' all.equal(x, RP)
#'
read_rp_xlsx <- function() {
    openxlsx::read.xlsx(pkg_file("extdata/RP.xlsx"), 1)
}

#' @export
#' @keywords dataset
#'
#' @title Read the HILIC dataset from the Retip package
#'
#' @description
#' Reads the `Retip::HILIC` dataset (CC BY 4.0) from the Retip package or, if
#' Retip is not installed, downloads the dataset directly from the [Retip GitHub
#' repository](https://github.com/oloBion/Retip). Before returning the dataset,
#' SMILES strings are canonicalized and the original `tibble` object is
#' converted to a base R `data.frame`.
#'
#' @param verbose Verbosity. 1 for messages, 0 to suppress them.
#'
#' @return A data frame with 970 rows and the following columns:
#' - `NAME`: Molecule name
#' - `INCHIKEY`: InChIKey
#' - `SMILES`: Canonical SMILES string
#' - `RT`: Retention time in Minutes
#'
#' @details
#' Attribution as required by CC BY 4.0:\cr
#' + Original dataset by:
#'   Paolo Bonini, Tobias Kind, Hiroshi Tsugawa, Dinesh Kumar Barupal, and
#'   Oliver Fiehn as part of the Retip project.\cr
#' + Source repository: <https://github.com/oloBion/Retip>\cr
#' + Original file: <https://github.com/oloBion/Retip/raw/master/data/HILIC.RData>\cr
#' + License: CC BY 4.0 (<https://creativecommons.org/licenses/by/4.0/>)\cr
#' + Modifications in FastRet:
#'   - converted tibble to data.frame
#'   - canonicalized SMILES using [as_canonical()]
#'   - renamed column 'INCHKEY' to 'INCHIKEY'
#'
#' @source
#' <https://github.com/oloBion/Retip/raw/master/data/HILIC.RData>
#'
#' @references
#' Retip: Retention Time Prediction for Compound Annotation in Untargeted Metabolomics\cr
#' Paolo Bonini, Tobias Kind, Hiroshi Tsugawa, Dinesh Kumar Barupal, and Oliver Fiehn\cr
#' Analytical Chemistry 2020 92 (11), 7515-7522 DOI: 10.1021/acs.analchem.9b05765
#'
#' @examples
#' df <- read_retip_hilic_data(verbose = 0)
read_retip_hilic_data <- function(verbose = 1) {
    logf <- if (verbose) catf else null
    if (is.data.frame(getOption("FastRet.HILIC_Retip"))) {
        logf("Reading Retip::HILIC from RAM cache")
        return(getOption("FastRet.HILIC_Retip"))
    } else if (eval(parse(text = "requireNamespace('Retip', quietly = TRUE)"))) {
        logf("Reading Retip::HILIC from Retip installation.")
        tbl <- eval(parse(text = "Retip::HILIC"))
        # We do this absurd eval trick to avoid R CMD check complaining about
        # Retip not being a suggested package. (We would suggest it, but it is
        # not on CRAN or any other CRAN-like repository, so we cannot suggest it
        # without triggering further warnings.)
    } else {
        logf("Retip package not found. Downloading from GitHub instead.")
        logf("Tip: install Retip from GitHub to avoid this download in the future.")
        logf("See https://github.com/oloBion/Retip for installation instructions.")
        # Encourage installation of Retip, to support the original project.
        url <- "https://github.com/oloBion/Retip/raw/master/data/HILIC.RData"
        destfile <- tempfile("HILIC", fileext = ".RData")
        download.file(url, destfile, mode = "wb", quiet = !verbose)
        HILIC <- NULL # will be loaded in the next line
        load(destfile)
        tbl <- HILIC
    }
    logf("Converting from tibble to data.frame")
    df <- as.data.frame(tbl)
    logf("Canonicalizing SMILES strings")
    df$SMILES <- as_canonical(df$SMILES)
    logf("Renaming column 'INCHKEY' to 'INCHIKEY'")
    names(df)[names(df) == "INCHKEY"] <- "INCHIKEY"
    logf("Caching Retip::HILIC in RAM for future calls")
    options(FastRet.HILIC_Retip = df) # cache for future calls
    df
}

#' @export
#' @keywords dataset
#'
#' @title Hypothetical retention times
#'
#' @description
#' Subset of the data from [read_rp_xlsx()] with some slight modifications to
#' simulate changes in temperature and/or flowrate.
#'
#' @return
#' A dataframe with 25 rows (metabolites) and 3 columns: RT, SMILES and NAME.
#'
#' @examples \donttest{
#' x <- read_rpadj_xlsx()
#' }
read_rpadj_xlsx <- function() {
    openxlsx::read.xlsx(pkg_file("extdata/RP_adj.xlsx"), 1)
}

#' @export
#' @keywords dataset
#' @title LASSO Model trained on RP dataset
#'
#' @description
#' Read a LASSO model trained on the [RP] dataset using [train_frm()].
#'
#' @return A `frm` object.
#'
#' @examples
#' frm <- read_rp_lasso_model_rds()
read_rp_lasso_model_rds <- function() {
    readRDS(pkg_file("extdata/RP_lasso_model.rds"))
}

## RP #####

#' @keywords dataset
#'
#' @title Retention Times (RT) measured on a Reverse Phase (RP) Column
#'
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
#'
#' @source
#' Measured by the Institute of Functional Genomics at the University of
#' Regensburg.
#'
#' @seealso read_rp_xlsx
"RP"

# Private #####

make_lazyload_objs <- function(HILIC, RP, RP_Mod, RP_Val) {
    # Read in data from Excel files
    RP <- read_rp_xlsx()
    # Store objects in the data folder
    usethis::use_data(RP, overwrite = TRUE)
}
