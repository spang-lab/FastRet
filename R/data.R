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
    xlsx::read.xlsx(pkg_file("extdata/RP.xlsx"), 1)
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
    xlsx::read.xlsx(pkg_file("extdata/RP_adj.xlsx"), 1)
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

# Read From Disk (Private) #####

read_rp_sub <- function() {
    path <- "misc/datasets/20211022_R8_dif_conditions_Medoids_validSet.xlsx"
    df <- xlsx::read.xlsx(
        file = pkg_file(path),
        sheetName = "R8_RT_Medoids",
        startRow = 2,
        colIndex = 2:18
    )
    # Original column names: PLATE, NROW, NCOL, CNAME, SMILES, FORMULA, M.H (M+H
    # in Excel), M.H.1 (M-H in Excel), HMDB, RT, Code, Steeper, Flatter, T25,
    # FR025, T25_FR025, T25_FR025_Steeper
    RP_Mod <- list(
        Steep = df[, c("CNAME", "SMILES", "Steeper")],
        Flat = df[, c("CNAME", "SMILES", "Flatter")],
        T25 = df[, c("CNAME", "SMILES", "T25")],
        FR25 = df[, c("CNAME", "SMILES", "FR025")],
        T25_FR25 = df[, c("CNAME", "SMILES", "T25_FR025")],
        T25_FR25_Steep = df[, c("CNAME", "SMILES", "T25_FR025_Steeper")]
    )
    for (x in names(RP_Mod)) {
        colnames(RP_Mod[[x]]) <- c("NAME", "SMILES", "RT")
    }
    RP_Mod
}

read_rp_val <- function() {
    path <- "misc/datasets/20211022_R8_dif_conditions_Medoids_validSet.xlsx"
    df <- xlsx::read.xlsx(
        file = pkg_file(path),
        sheetName = "R8_RT_Validation set",
    )
    # Original column names: SMILES, Name, RT.Normal, RT.Steep, RT.Flatter, T25,
    # FR025, T25_FR025, T25_FR025_Steeper
    RP_Val <- list(
        Normal = df[, c("Name", "SMILES", "RT.Normal")],
        Steep = df[, c("Name", "SMILES", "RT.Steep")],
        Flat = df[, c("Name", "SMILES", "RT.Flatter")],
        T25_Flat = df[, c("Name", "SMILES", "T25")],
        FR25_Flat = df[, c("Name", "SMILES", "FR025")],
        T25_FR25_Flat = df[, c("Name", "SMILES", "T25_FR025")],
        T25_FR25_Steep = df[, c("Name", "SMILES", "T25_FR025_Steeper")]
    )
    for (x in names(RP_Val)) {
        colnames(RP_Val[[x]]) <- c("NAME", "SMILES", "RT")
    }
    RP_Val
}

read_hilic_xlsx <- function() {
    file <- "misc/datasets/20210702_RT_Prediction_Hilic_Library_FF.xlsx"
    df <- xlsx::read.xlsx(pkg_file(file), 1)
    df[, c("NAME", "SMILES", "RT")]
}

read_rp_axmm_xlsx <- function() {
    stop("Excel file with measurements from RP-AXMM column is missing")
}

# Lazyload Objects (Private) #####

make_lazyload_objs <- function(HILIC, RP, RP_Mod, RP_Val) {

    # Read in data from Excel files
    HILIC <- read_hilic_xlsx()
    RP <- read_rp_xlsx()
    RP_Mod <- read_rp_sub()
    RP_Val <- read_rp_val()

    # Prepend the normal conditions to the RP_Mod datasets
    mod_idx <- match(RP_Mod[[1]]$SMILES, RP$SMILES)
    RP_Sub <- RP[mod_idx, c("NAME", "SMILES", "RT")]
    rownames(RP_Sub) <- NULL
    RP_Mod <- c(list(Normal = RP_Sub), RP_Mod)

    # Store objects in the data folder
    usethis::use_data(HILIC, overwrite = TRUE)
    usethis::use_data(RP, overwrite = TRUE)
    usethis::use_data(RP_Mod, overwrite = TRUE)
    usethis::use_data(RP_Val, overwrite = TRUE)
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

## HILIC #####

#' @title Retention Times (RT) measured on a Hydrophilic Interaction Chromatography (HILIC) Column
#' @description Retention time data from a hydrophilic interaction chromatography experiment.
#' @format A dataframe of 393 metabolites with the following columns:
#' \describe{
#'   \item{NAME}{Name of the metabolite}
#'   \item{SMILES}{SMILES notation of the metabolite}
#'   \item{RT}{Retention time}
#' }
#' @source
#' Measured by the Institute of Functional Genomics at the University of
#' Regensburg.
#' @keywords dataset
"HILIC"

## RP_Mod #####

#' @title Subset of the 'RP' dataset measured under different conditions
#' @description
#' A subset of 25 metabolites from the original [RP] dataset was re-measured
#' under modified chromatographic conditions (gradient, temperature, and/or flow
#' rate). This object contains the results as data frames. See 'Details' for an
#' overview of the modifications.
#' @format
#' A named list with elements: 'Normal', 'Steep', 'Flat', 'T25_Flat',
#' 'FR25_Flat', 'T25_FR25_Flat', and 'T25_FR25_Steep'. Each element is a data
#' frame with 25 metabolites (rows) and the following columns:
#' \describe{
#'   \item{RT}{Retention time}
#'   \item{SMILES}{SMILES notation of the metabolite}
#'   \item{NAME}{Name of the metabolite}
#' }
#' The 'Normal' condition is a true subset of the original [RP] dataset: it
#' contains the same metabolites with identical retention times as in [RP]. All
#' other elements ('Steep', 'Flat', etc.) contain the same metabolites, but
#' their retention times were measured under different chromatographic
#' conditions.
#' @details
#' The following conditions were used during the experiment:
#'
#' | Name            | Temperature | Flow Rate | Gradient |
#' | --------------- | ----------- | --------- | -------- |
#' | Normal          | 35          | 0.3       | normal   |
#' | Steep           | 35          | 0.3       | steep    |
#' | Flat            | 35          | 0.3       | flat     |
#' | T25_Flat        | 25          | 0.3       | flat     |
#' | FR025_Flat      | 35          | 0.25      | flat     |
#' | T25_FR025_Flat  | 25          | 0.25      | flat     |
#' | T25_FR025_Steep | 25          | 0.25      | steep    |
#'
#' For all conditions, formic acid 0.1% (v/v) in Water was used as mobile phase
#' A and formic acid 0.1% (v/v) in ACN as mobile phase B. Temperatures are given
#' in \eqn{^\circ}C and flow rates in mL/min. For the exact definitions of a
#' 'normal', 'steep' and 'flat' gradient please refer to the FastRet
#' publication.
#' @source
#' Measured by the Institute of Functional Genomics at the University of
#' Regensburg.
#' @keywords dataset
"RP_Mod"

## RP_Val #####

#' @title Validation Dataset
#' @description
#' A set of 22 metabolites measured on the same Reverse Phase (RP)
#' chromatography column that was used for the [RP] and [RP_Mod] datasets. None
#' of the 22 metabolites is included in the [RP] dataset. The measurements were
#' performed under varying conditions, such as different temperatures, flow
#' rates and/or gradients. For a detailed description of the applied conditions,
#' see section 'Details'.
#' @format
#' A named list with elements: 'Normal', 'Steep', 'Flat', 'T25_Flat',
#' 'FR25_Flat', 'T25_FR25_Flat', and 'T25_FR25_Steep'. Each element is a data
#' frame with 25 metabolites (rows) and the following columns:
#' \describe{
#'   \item{RT}{Retention time}
#'   \item{SMILES}{SMILES notation of the metabolite}
#'   \item{NAME}{Name of the metabolite}
#' }
#' @details
#' The following conditions were used during the experiment:
#'
#' | Name            | Temperature | Flow Rate | Gradient |
#' | --------------- | ----------- | --------- | -------- |
#' | Normal          | 35          | 0.3       | normal   |
#' | Steep           | 35          | 0.3       | steep    |
#' | Flat            | 35          | 0.3       | flat     |
#' | T25_Flat        | 25          | 0.3       | flat     |
#' | FR025_Flat      | 35          | 0.25      | flat     |
#' | T25_FR025_Flat  | 25          | 0.25      | flat     |
#' | T25_FR025_Steep | 25          | 0.25      | steep    |
#'
#' For all conditions, formic acid 0.1% (v/v) in Water was used as mobile phase
#' A and formic acid 0.1% (v/v) in ACN as mobile phase B. Temperatures are given
#' in \eqn{^\circ}C and flow rates in mL/min. For the exact definitions of a
#' 'normal', 'steep' and 'flat' gradient please refer to the FastRet
#' publication.
#' @source
#' Measured by the Institute of Functional Genomics at the University of
#' Regensburg.
#' @keywords dataset
"RP_Val"
