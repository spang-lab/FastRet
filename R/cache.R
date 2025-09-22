# RAM Cache (RC) #####

#' @export
#' @keywords internal
#'
#' @title RAM Cache Environment
#'
#' @description
#' An environment used for caching data in RAM.
#'
#' @details
#' This environment is used by [getCDsFor1Molecule()] to store the results of
#' previous calculations to speed up subsequent calls. It gets initalized upon
#' the first call of [getCDsFor1Molecule()] with the chemical descriptors for
#' all molecules available in the [RP] dataset and the `HILIC` dataset of the
#' [Retip](https://www.retip.app/) package.
#'
#' @format
#' An environment with the following elements:
#' - `CDs`: A named list where each element contains the chemical descriptors
#'   for one SMILES string. The names of the list are SMILES strings, and each
#'   element is a data.frame with one row containing the chemical descriptors
#'   listed in [CDFeatures].
#'
#' @references
#' Retip: Retention Time Prediction for Compound Annotation in Untargeted
#' Metabolomics Paolo Bonini, Tobias Kind, Hiroshi Tsugawa, Dinesh Kumar
#' Barupal, and Oliver Fiehn Analytical Chemistry 2020 92 (11), 7515-7522 DOI:
#' 10.1021/acs.analchem.9b05765
#'
#' @examples
#' dim(ram_cache$CDs) # 1316 241
ram_cache <- as.environment(list(CDs = data.frame())) # Initialized in .onLoad()

#' @noRd
#' @title Initialize RAM Cache
#' @description
#' Initialize the RAM cache with pre-computed chemical descriptors from the
#' package data.
#' @param logf Function for logging messages
rc_init <- function(logf = null) {
    if (nrow(ram_cache$CDs) > 0) {
        pkgCDs <- readRDS(pkg_file("cachedata/CDs.rds"))
        newCDs <- pkgCDs[rownames(pkgCDs) %notin% rownames(ram_cache$CDs), ]
        CDs <- rbind(ram_cache$CDs, newCDs)
    } else {
        CDs = readRDS(pkg_file("cachedata/CDs.rds"))
    }
    assign("CDs", CDs, envir = ram_cache)
}

#' @noRd
#' @title Clear RAM Cache
#' @description
#' Clear the RAM cache either completely or for specific SMILES strings.
#' @param smiles Vector of SMILES strings to remove, or NULL to clear all
rc_clear <- function(smiles = NULL) {
    if (is.null(smiles)) {
        assign("CDs", ram_cache$CDs[0, ], envir = ram_cache)
    } else {
        keep <- rownames(ram_cache$CDs) %notin% smiles
        assign("CDs", ram_cache$CDs[keep, ], envir = ram_cache)
    }
}

#' @noRd
#' @title Get Chemical Descriptors from RAM Cache
#' @description
#' Retrieve chemical descriptors for a SMILES string from the RAM cache.
#' @param smi SMILES strings
#' @param logf Function for logging messages
#' @return Data frame with chemical descriptors
rc_get <- function(smi, logf = null) {
    n_smi <- length(smi)
    n_miss <- sum(smi %notin% rownames(ram_cache$CDs))
    if (n_miss > 0) {
        logf("CDs for %d/%d SMILES not in RAM cache. Returning NULL.", n_miss, n_smi)
        return(NULL)
    } else {
        logf("Loading chemical descriptors for %d SMILES from RAM cache", smi)
    }
    X <- ram_cache$CDs[smi, ]
}

#' @noRd
#' @title Write Chemical Descriptors to RAM Cache
#' @description
#' Store chemical descriptors for a SMILES string in the RAM cache.
#' @param smi SMILES string
#' @param cds Chemical descriptors data frame (1 SMILES x 241 CDs)
#' @param logf Function for logging messages
rc_set <- function(smi, cds, logf = null) {
    logf("Storing chemical descriptors for '%s' in RAM cache", smi)
    CDs <- ram_cache$CDs # (1)
    CDs[smi, ] <- cds
    assign("CDs", CDs, envir = ram_cache)
    # (1) We need to make a local copy first, as we cannot modify objects inside
    # ram_cache directly. That's very unfortunate, as copying the whole
    # data.frame every time is inefficient. But since we only update the RAM
    # cache in blocks in [getCDs()] it's acceptable for now.
}

#' @noRd
#' @title Make cachedata/CDs.rds
#' @description
#' The FastRet package comes with a set of pre-calculated chemical descriptors
#' for metabolites stored in cachedata/CDS.rds, allowing all examples, tests
#' etc. to run much faster. This function creates cachedata/CDs.rds.
#' @param nw Number of workers for parallel processing
#' @param nsmi Number of SMILES to process. NULL means "all". Useful to speed up
#' testing.
#' @param overwrite Whether to overwrite existing cachedata/CDs.rds file
#' @details
#' Since this function is used to create the cachedata/CDs.rds file, that is
#' shipped with the package, it MUST be called during package development, i.e.
#' after cloning the package sources and invoking [devtools::load_all()]. It
#' makes no sense to call if after installation and loading of the package via
#' [library()].
update_cachedata_cds <- function(nw = 16, nsmi = NULL, overwrite = FALSE) {
    rc_clear()
    hilic <- read_retip_hilic_data()
    df <- rbind(FastRet::RP, hilic[colnames(FastRet::RP)])
    smis <- unname(unique(df$SMILES))
    if (!is.null(nsmi)) { # Mostly for Testing
        catf("ATTENTION: Only the first %d SMILES will be processed.", nsmi)
        smis <- smis[seq_len(nsmi)]
    }
    cds <- getCDsFromCDK(smis, nw = nw, null)
    cachedata_path <- pkg_file("cachedata")
    rdspath <- file.path(cachedata_path, "CDs.rds")
    if (file.exists(rdspath)) {
        if (overwrite) {
            catf("Overwriting '%s'.", rdspath)
            saveRDS(cds, rdspath)
        } else {
            catf("File '%s' already exists. Not overwriting.", rdspath)
        }
    } else {
        catf("Writing '%s'.", rdspath)
        saveRDS(cds, rdspath)
    }
    invisible(cds)
}

# Disk Cache #####

#' @export
#' @title Get cache directory
#' @description Creates and returns the cache directory for the FastRet package.
#' @param subdir Optional subdirectory within the cache directory.
#' @return The path to the cache directory or subdirectory.
#' @keywords internal
#' @examples
#' path <- get_cache_dir()
get_cache_dir <- function(subdir = NULL, obj = NULL) {
    # This function was already exported with before the disk_cache_* function
    # were introduced in 1.2, so the name doesn't fit the disk_cache_* pattern,
    # but renaming is not worth the breaking change.
    cache_dir <- tools::R_user_dir("FastRet", which = "cache")
    if (!is.null(subdir)) cache_dir <- file.path(cache_dir, subdir)
    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    normalizePath(cache_dir, winslash = "/", mustWork = FALSE)
}

#' @noRd
#' @title Get Cache File Path
#' @description
#' Generate the file path for caching a specific object.
#' @param obj Object to be cached (used for generating hash)
#' @param func Function name for cache directory organization
#' @return File path for cache file
dc_path <- function(obj, func = "getCDsFor1Molecule") {
    dir <- get_cache_dir(func)
    hash <- digest::digest(obj, algo = "md5")
    name <- paste0(hash, ".rds")
    file.path(dir, name)
}

#' @noRd
#' @title Clear Disk Cache
#' @description
#' Clear the disk cache either completely or for specific SMILES strings.
#' @param smiles Vector of SMILES strings to remove, or NULL to clear all
dc_clear <- function(smiles = NULL) {
    path <- get_cache_dir("getCDsFor1Molecule")
    smiles <- unique(smiles)
    if (is.null(smiles)) {
        unlink(path, recursive = TRUE)
        invisible(path)
    } else {
        paths <- vapply(smiles, dc_path, character(1), func = "getCDsFor1Molecule")
        file.remove(paths[file.exists(paths)])
        invisible(paths)
    }
}

#' @noRd
#' @title Get Chemical Descriptors from Disk Cache
#' @description
#' Retrieve chemical descriptors for a SMILES string from the disk cache.
#' @param smi SMILES string
#' @param logf Function for logging messages
#' @return Data frame with chemical descriptors or NULL if not found
dc_get <- function(smi, logf = null) {
    f <- dc_path(smi, "getCDsFor1Molecule")
    if (file.exists(f)) {
        logf("Loading chemical descriptors for '%s' from disk cache", smi)
        return(readRDS(f))
    } else {
        logf("Chemical descriptors for '%s' not found in disk cache", smi)
        return(NULL)
    }
}

#' @noRd
#' @title Write Chemical Descriptors to Disk Cache
#' @description
#' Store chemical descriptors for a SMILES string in the disk cache.
#' @param smi SMILES string
#' @param cds Chemical descriptors data frame
#' @param catf Function for logging messages
#' @param overwrite Whether to overwrite existing cache files
dc_set <- function(smi, cds, logf = null, overwrite = FALSE) {
    path <- dc_path(smi)
    logf("Storing chemical descriptors for '%s' at '%s'", smi, path)
    if (!file.exists(path) || overwrite) {
        dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
        saveRDS(cds, path)
    } else {
        logf("File '%s' already exists. Not overwriting.", path)
    }
    invisible(NULL)
}
