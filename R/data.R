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

read_rp_mod <- function() {
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
        T25_Flat = df[, c("CNAME", "SMILES", "T25")],
        FR25_Flat = df[, c("CNAME", "SMILES", "FR025")],
        T25_FR25_Flat = df[, c("CNAME", "SMILES", "T25_FR025")],
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
    path <- "misc/datasets/RP-AXMM_FastRet_Input.xlsx"
    df <- xlsx::read.xlsx(pkg_file(path), 1)
    df[, c("NAME", "SMILES", "RT")]
}

# Hardcoded (Private) #####

make_columns_df <- function() {
    Columns = data.frame(ID = character(4), Description = character(4))
    Columns[1, ] = c("RP",          "ACQUITY Premier HSS T3 (1.8 µm 150 x 2.1 mm id)")
    Columns[2, ] = c("RP_AXMM",     "Atlantis PREMIER BEH C18 AX (1.7 µm 150 x 2.1 mm id)")
    Columns[3, ] = c("HILIC",       "ACQUITY PREMIER BEH Amide (1.7 µm 150 x 2.1 mm id)")
    Columns[4, ] = c("HILIC_Retip", "TODO")
    Columns
}

make_eluents_df <- function() {
    Eluents = data.frame(ID = character(6), Description = character(6))
    Eluents[1, ] = c("FA100_W",         "Formic Acid 0.1%   (v/v) in Water")
    Eluents[2, ] = c("FA100_A",         "Formic Acid 0.1%   (v/v) in ACN")
    Eluents[3, ] = c("FA100_AF10_W",    "Formic Acid 0.1%   (v/v) and 10 mM Ammonium Formate in Water")
    Eluents[4, ] = c("FA125_AF10_W",    "Formic Acid 0.125% (v/v) and 10 mM Ammonium Formate in Water")
    Eluents[5, ] = c("FA100_AF10_WA10", "Formic Acid 0.1%   (v/v) and 10 mM Ammonium Formate in Water:ACN (10:90, v/v)")
    Eluents[6, ] = c("FA125_AF10_WA5",  "Formic Acid 0.125% (v/v) and 10 mM Ammonium Formate in Water:ACN (5:95, v/v)")
    Eluents
}

make_gradients_df <- function() {
    # All Gradients are given as B% (i.e. the percentage of mobile phase B)
    # Pretty-Printable via following command:
    # pander::pandoc.table(Gradients, justify="right", keep.line.breaks=TRUE)
    Gradients = data.frame(ID = character(6), Description = character(6))
    Gradients[1, "ID"] <- "RP_Normal"
    Gradients[1, "Description"] <- paste(sep = "\n",
        "  3% to  50%  in  8.0 min",
        " 50% to 100%  in  1.0 min",
        "        100% for  3.5 min",
        "100% to   0%  in  1.0 min",
        "          0% for  8.5 min"
    )
    Gradients[2, "ID"] <- "RP_Steep"
    Gradients[2, "Description"] <- paste(sep = "\n",
        "  3% to  50%  in  5.0 min",
        " 50% to 100%  in  4.0 min",
        "        100% for  3.5 min",
        "100% to   0%  in  1.0 min",
        "          0% for  8.5 min"
    )
    Gradients[3, "ID"] <- "RP_Flat"
    Gradients[3, "Description"] <- paste(sep = "\n",
        "  3% to  30%  in  8.0 min",
        " 30% to 100%  in  1.0 min",
        "        100% for  3.5 min",
        "100% to   0%  in  1.0 min",
        "          0% for  8.5 min"
    )
    Gradients[4, "ID"] <- "RPAXMM_Normal"
    Gradients[4, "Description"] <- paste(sep = "\n",
        "  0% to  30%  in 10.0 min",
        " 30% to 100%  in  5.0 min",
        "        100% for  5.0 min",
        "100% to   0%  in  2.0 min",
        "          0% for  8.0 min"
    )
    Gradients[5, "ID"] <- "HILIC_Normal"
    Gradients[5, "Description"] <- paste(sep = "\n",
        " 97% to  70%  in  8.0 min",
        " 70% to  40%  in  2.0 min",
        "         40% for  1.5 min",
        " 40% to  97%  in  1.5 min",
        "         97% for  6.0 min"
    )
    Gradients[6, "ID"] <- "HILIC_Retip_Normal"
    Gradients[6, "Description"] <- "TODO"
    Gradients
}

# The temperature is given in degree Celsius.
# The flow rate is given in ml/min.
# The gradient is given as percentage of mobile phase B (%B) in the eluent (which is a mixture of mobile phase A and B).
make_conditions_df <- function() {
    Conditions <- data.frame(
        ID = character(10),
        Temp = numeric(10),
        FlowRate = numeric(10),
        Gradient = character(10),
        PhaseA = character(10),
        PhaseB = character(10)
    )
    Conditions[1, ]  <- c("RP_Normal",          35,  0.30,     "RP_Normal",     "FA100_W",      "FA100_A")
    Conditions[2, ]  <- c("RP_Steep",           35,  0.30,     "RP_Steep",      "FA100_W",      "FA100_A")
    Conditions[3, ]  <- c("RP_Flat",            35,  0.30,     "RP_Flat",       "FA100_W",      "FA100_A")
    Conditions[4, ]  <- c("RP_T25_Flat",        25,  0.30,     "RP_Flat",       "FA100_W",      "FA100_A")
    Conditions[5, ]  <- c("RP_FR25_Flat",       35,  0.25,     "RP_Flat",       "FA100_W",      "FA100_A")
    Conditions[6, ]  <- c("RP_T25_FR25_Flat",   25,  0.25,     "RP_Flat",       "FA100_W",      "FA100_A")
    Conditions[7, ]  <- c("RP_T25_FR25_Steep",  25,  0.25,     "RP_Steep",      "FA100_W",      "FA100_A")
    Conditions[8, ]  <- c("HILIC",              45,  0.40,     "HILIC_Normal",  "FA125_AF10_W", "FA125_AF10_WA5")
    Conditions[9, ]  <- c("RPAXMM",             35,  0.30,     "RPAXMM_Normal", "FA100_AF10_W", "FA100_AF10_WA10")
    Conditions[10, ] <- c("HILIC_Retip_Normal", NA,  NA,       NA,              NA,             NA)
    Conditions
}

#' @noRd
#' @title Get Normalized Datasets
#' @details Runtime: 1.26 seconds
#' @examples
#' make_datasets_list(update_cache = TRUE) # update cache
#' dsl <- make_datasets_list()
#' str(dsl, 1)
make_datasets_list <- function(read_cache = TRUE, update_cache = FALSE) {
    cachefile <- file.path(get_cache_dir("make_datasets_list"), "nds.rds")
    if (!update_cache && read_cache && file.exists(cachefile)) {
        return(readRdsVerbose(cachefile))
    }
    dsl <- list(
        RP = FastRet::RP,
        RP_Steep = FastRet::RP_Mod$Steep,
        RP_Flat = FastRet::RP_Mod$Flat,
        RP_T25_Flat = FastRet::RP_Mod$T25_Flat,
        RP_FR25_Flat = FastRet::RP_Mod$FR25_Flat,
        RP_T25_FR25_Flat = FastRet::RP_Mod$T25_FR25_Flat,
        RP_T25_FR25_Steep = FastRet::RP_Mod$T25_FR25_Steep,
        RP_Val = FastRet::RP_Val$Normal,
        RP_Val_Steep = FastRet::RP_Val$Steep,
        RP_Val_Flat = FastRet::RP_Val$Flat,
        RP_Val_T25_Flat = FastRet::RP_Val$T25_Flat,
        RP_Val_FR25_Flat = FastRet::RP_Val$FR25_Flat,
        RP_Val_T25_FR25_Flat = FastRet::RP_Val$T25_FR25_Flat,
        RP_Val_T25_FR25_Steep = FastRet::RP_Val$T25_FR25_Steep,
        RP_AXMM = FastRet::RP_AXMM,
        HILIC = FastRet::HILIC,
        HILIC_Retip = FastRet::read_retip_hilic_data()
    )
    ids <- names(dsl)
    cols <- c("NUM", "RT", "NAME", "DUPOF", "SMILES", "CANONICAL")
    for (i in seq_along(ids)) {
        catf("Processing dataset %s (%d/%d)", ids[i], i, length(ids))
        dsl[[i]]$NUM <- seq_len(nrow(dsl[[i]]))
        dsl[[i]]$CANONICAL <- get_canonical_smiles(dsl[[i]]$SMILES) # slow part
        dsl[[i]]$DUPOF <- get_dupof(dsl[[i]])
        dsl[[i]] <- as.data.frame(dsl[[i]][, cols])
    }
    if (update_cache || (read_cache && !file.exists(cachefile))) {
        saveRdsVerbose(dsl, cachefile)
    }
    invisible(dsl)
}

# Relationship between datasets, columns, conditions, eluents and gradients
# - 1 Dataset
#   - 1 ID
#   - 1 Purpose
#   - 1 Column
#     - 1 ID
#     - 1 Description
#   - 1 Condition
#     - 1 ID
#     - 1 Temperature
#     - 1 Flow Rate
#     - 1 Gradient
#     - 1 Mobile Phase A (Eluent)
#     - 1 Mobile Phase B (Eluent)
#   - n Measurments
#     - 1 NUM
#     - 1 RT
#     - 1 NAME
#     - n DUPOF (Duplicate of)
#     - 1 SMILES (Input SMILES)
#     - 1 CANON (Canonical SMILES)
make_datasets_df <- function(dsl = make_datasets_list()) {
    ds <- data.frame(
        ID = character(0),
        Purpose = character(0),
        ColumnID = character(0),
        ConditionID = character(0),
        nMeas = integer(0),
        nNames = integer(0),
        nSMILES = integer(0),
        nMets = integer(0)
    )
    #              ID                       Purpose,      ColumnID       ConditionID
    ds[1, ]  <- c("RP",                    "Training",    "RP",          "RP_Normal"          )
    ds[2, ]  <- c("RP_Steep",              "Adjustment",  "RP",          "RP_Steep"           )
    ds[3, ]  <- c("RP_Flat",               "Adjustment",  "RP",          "RP_Flat"            )
    ds[4, ]  <- c("RP_T25_Flat",           "Adjustment",  "RP",          "RP_T25_Flat"        )
    ds[5, ]  <- c("RP_FR25_Flat",          "Adjustment",  "RP",          "RP_FR25_Flat"       )
    ds[6, ]  <- c("RP_T25_FR25_Flat",      "Adjustment",  "RP",          "RP_T25_FR25_Flat"   )
    ds[7, ]  <- c("RP_T25_FR25_Steep",     "Adjustment",  "RP",          "RP_T25_FR25_Steep"  )
    ds[8, ]  <- c("RP_Val",                "Validation",  "RP",          "RP_Normal"          )
    ds[9, ]  <- c("RP_Val_Steep",          "Validation",  "RP",          "RP_Steep"           )
    ds[10, ] <- c("RP_Val_Flat",           "Validation",  "RP",          "RP_Flat"            )
    ds[11, ] <- c("RP_Val_T25_Flat",       "Validation",  "RP",          "RP_T25_Flat"        )
    ds[12, ] <- c("RP_Val_FR25_Flat",      "Validation",  "RP",          "RP_FR25_Flat"       )
    ds[13, ] <- c("RP_Val_T25_FR25_Flat",  "Validation",  "RP",          "RP_T25_FR25_Flat"   )
    ds[14, ] <- c("RP_Val_T25_FR25_Steep", "Validation",  "RP",          "RP_T25_FR25_Steep"  )
    ds[15, ] <- c("RP_AXMM",               "Training",    "RP_AXMM",     "HILIC"              )
    ds[16, ] <- c("HILIC",                 "Training",    "HILIC",       "RPAXMM"             )
    ds[17, ] <- c("HILIC_Retip",           "Training",    "HILIC_Retip", "HILIC_Retip_Normal" )
    ds[, "nMeas"] <- sapply(dsl, nrow)
    ds[, "nNames"] <- sapply(dsl, function(x) length(unique(x$NAME)))
    ds[, "nSMILES"] <- sapply(dsl, function(x) length(unique(x$SMILES)))
    ds[, "nMets"] <- sapply(dsl, function(x) length(unique(x$CANONICAL)))
    ds
}

make_measurements_df <- function(dsl = make_datasets_list()) {
    meas <- do.call(rbind, dsl)
    meas$GLOBNUM <- seq_len(nrow(meas))
    meas$DSID <- sapply(strsplit(rownames(meas), ".", fixed = TRUE), "[", 1)
    meas$DSINUM <- meas$NUM
    meas[, c("GLOBNUM", "DSID", "DSINUM", "RT", "NAME", "DUPOF", "SMILES", "CANONICAL")]
}

make_metabolites_df <- function(dsl = make_datasets_list()) {
    n <- length(dsl)
    canonicals <- unique(unlist(sapply(dsl, `[[`, "CANONICAL")))
    counts <- data.frame(CANONICAL = canonicals)
    counts[2:(n+1)] <- as.data.frame(sapply(dsl, function(x)
        table(factor(x, levels = canonicals)
    )))
    counts
}

# Lazyload Objects (Private) #####

make_lazyload_objs <- function(HILIC, RP, RP_Mod, RP_Val) {
    # Read in data from Excel files
    HILIC <- read_hilic_xlsx()
    RP <- read_rp_xlsx()
    RP_Mod <- read_rp_mod()
    RP_Val <- read_rp_val()
    RP_AXMM <- read_rp_axmm_xlsx()
    # Store objects in the data folder
    usethis::use_data(HILIC, overwrite = TRUE)
    usethis::use_data(RP, overwrite = TRUE)
    usethis::use_data(RP_Mod, overwrite = TRUE)
    usethis::use_data(RP_Val, overwrite = TRUE)
    usethis::use_data(RP_AXMM, overwrite = TRUE)
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

## RP_AXMM #####

#' @title
#' Retention Times (RT) measured on a Reverse Phase Anion Exchange Mixed Mode
#' (RP-AXMM) Column
#' @description
#' Retention time data from a Reverse Phase Anion Exchange Mixed Mode (RP-AXMM)
#' Column.
#' @format A dataframe of 438 metabolites with the following
#' columns:
#' \describe{
#'   \item{SMILES}{SMILES notation of the metabolite}
#'   \item{NAME}{Name of the metabolite}
#'   \item{RT}{Retention time}
#' }
#' @source
#' Measured by the Institute of Functional Genomics at the University of
#' Regensburg.
#' @keywords dataset
#' @seealso read_rp_xlsx
"RP_AXMM"

## HILIC #####

#' @title
#' Retention Times (RT) measured on a Hydrophilic Interaction Chromatography
#' (HILIC) Column
#' @description
#' Retention time data from a hydrophilic interaction chromatography experiment.
#' @format A dataframe of 392 metabolites with the following columns:
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
#' A named list with elements: 'Steep', 'Flat', 'T25_Flat',
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
#' in \eqn{^\circ}C and flow rates in mL/min. The 'Normal' condition was used to
#' measured the original `RP` dataset and is not included in the `RP_Mod`
#' object. For the exact definitions of a 'normal', 'steep' and 'flat' gradient
#' please refer to the FastRet publication.
#'
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
