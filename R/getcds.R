# Public #####

#' @export
#' @keywords public
#'
#' @title Get Chemical Descriptors for a list of molecules
#'
#' @description
#' Calculate Chemical Descriptors (CDs) for a list of molecules. Molecules can
#' appear multiple times in the list.
#'
#' @param df
#' dataframe with two mandatory columns: "NAME" and "SMILES"
#'
#' @param verbose
#' 0: no output, 1: progress, 2: more progress and warnings
#'
#' @param nw
#' number of workers for parallel processing
#'
#' @param keepdf
#' If TRUE, `cbind(df, CDs)` is returned. Else `CDs`.
#'
#' @return
#' A dataframe with all input columns (if `keepdf` is TRUE) and chemical
#' descriptors as remaining columns.
#'
#' @examples
#' cds <- getCDs(head(RP, 3), verbose = 1, nw = 1)
#'
getCDs <- function(df,
                   verbose = 1,
                   nw = 1,
                   keepdf = TRUE) {
    catf <- if (verbose >= 1) catf else null
    a <- Sys.time()
    cachedCDs <- getOption("FastRet.cachedCDs")
    if (is.null(cachedCDs) || nrow(cachedCDs) < 1788) {
        cachedCDs <- readRDS(pkg_file("cachedata/CDs.rds"))
        options(FastRet.cachedCDs = cachedCDs)
    }
    smi <- unique(df$SMILES[df$SMILES %notin% rownames(cachedCDs)])
    if (length(smi) > 0) {
        catf("Computing CDs for %d new SMILES using rCDK", length(smi))
        if (nw > 1) {
            chunkids <- cut(seq_along(smi), nw, labels = FALSE)
            DF <- split(data.frame(SMILES = smi), chunkids)
            CDS <- parLapply2(nw, DF, getCDs, verbose = 0, nw = 1, keepdf = FALSE)
            newCDs <- as.data.frame(data.table::rbindlist(CDS, use.names = TRUE))
            rownames(newCDs) <- smi
        } else {
            objs <- withr::with_options(list(warn = 2), rcdk::parse.smiles(smi))
            lapply(objs, rcdk::convert.implicit.to.explicit) # Add H atoms
            lapply(objs, rcdk::generate.2d.coordinates) # Gen 2D coords
            newCDs <- suppressWarnings(rcdk::eval.desc(objs, CDNames, verbose = FALSE))
        }
        cachedCDs <- rbind(cachedCDs, newCDs) # Add new CDs to already existing ones
        options(FastRet.cachedCDs = cachedCDs) # Update cached CDs
    }
    cds <- cachedCDs[df$SMILES, ]
    retdf <- if (keepdf) cbind(df, cds) else cds
    b <- Sys.time()
    catf("Finished calculating chemical descriptors in %s", format(b - a))
    invisible(retdf)
}

# Constants #####

#' @noRd
#' @title Make cachedata/CDs.rds
#' @description
#' The FastRet package comes with a set of pre-calculated chemical descriptors
#' for metabolites stored in cachedata/CDS.rds, allowing all examples, tests
#' etc. to run much faster. This function creates cachedata/CDs.rds.
#' @details
#' Since this function is used to create the cachedata/CDs.rds file, that is
#' shipped with the package, it MUST be called during package development, i.e.
#' after cloning the package sources and invoking [devtools::load_all()]. It
#' makes no sense to call if after installation and loading of the package via
#' [library()].
updateCachedCDs <- function() {
    cols <- c("NAME", "SMILES", "RT")
    hilic <- read_retip_hilic_data()[, cols]
    meas8 <- openxlsx::read.xlsx(pkg_file("extdata/Measurements_v8.xlsx"))[, cols]
    RPold <- openxlsx::read.xlsx(pkg_file("extdata/RP.xlsx"))[, cols]
    df <- rbind(hilic, meas8, RPold)
    df <- df[!duplicated(df$SMILES), ]
    CDs <- getCDs(df, verbose = 1, nw = 8, keepdf = FALSE)
    cachedata_path <- pkg_file("cachedata")
    rdspath <- file.path(cachedata_path, "CDs.rds")
    saveRDS(CDs, rdspath)
}


#' @export
#' @keywords internal
#'
#' @title Chemical Descriptors Names
#'
#' @description
#' This object contains the names of various chemical descriptors.
#'
#' @details
#' One descriptor can be associated with multiple features, e.g. the BCUT
#' descriptor corresponds to the following features: BCUTw.1l, BCUTw.1h,
#' BCUTc.1l, BCUTc.1h, BCUTp.1l, BCUTp.1h. Some descriptors produce warnings for
#' certain molecules., e.g. "The AtomType null could not be found" or "Molecule
#' must have 3D coordinates" and return NA in such cases. Descriptors that
#' produce only NAs in our test datasets will be excluded. To see which
#' descriptors produce only NAs, run `analyzeCDNames`. The
#' "LongestAliphaticChain" descriptors sometimes even produces `Error: segfault
#' from C stack overflow` error, e.g. for SMILES `c1ccccc1C(Cl)(Cl)Cl` (==
#' `rcdk::bpdata$SMILES[200]`) when using `OpenJDK Runtime Environment (build
#' 11.0.23+9-post-Ubuntu-1ubuntu122.04.1)`. Therefore, this descriptor is also
#' excluded.
#'
#' @seealso [analyzeCDNames()], [CDFeatures]
#'
#' @examples
#' str(CDNames)
#'
CDNames <- c(
    "org.openscience.cdk.qsar.descriptors.molecular.FractionalCSP3Descriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.SmallRingDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.FractionalPSADescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.MDEDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.HybridizationRatioDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.FMFDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.ChiClusterDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.BPolDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability",
    "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass",
    "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge",
    "org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.AminoAcidCountDescriptor"
)

updateCDNames <- function() {
    x <- rcdk::get.desc.names(type = "all")
    skipPattern <- paste0(
        "(WHIM",
        "|VABC",
        "|MomentOfInertia",
        "|LengthOverBreadth",
        "|GravitationalIndex",
        "|CPSA",
        "|TaeAminoAcid",
        "|LongestAliphaticChain)"
    )
    x[!grepl(skipPattern, x)]
}

#' @export
#' @keywords internal
#'
#' @title Analyze Chemical Descriptors Names
#'
#' @description
#' Analyze the chemical descriptor names and return a dataframe with their names
#' and a boolean column indicating if all values are NA.
#'
#' @details
#' This function is used to analyze the chemical descriptor names and to
#' identify which descriptors produce only NAs in the test datasets. The
#' function is used to generate the CDNames object.
#'
#' @param df
#' dataframe with two mandatory columns: "NAME" and "SMILES"
#'
#' @param descriptors
#' Vector of chemical descriptor names
#'
#' @return
#' A dataframe with two columns `descriptor` and `all_na`. Column `descriptor`
#' contains the names of the chemical descriptors. Column `all_na` contains a
#' boolean value indicating if all values obtained for the corresponding
#' descriptor are NA.
#'
#' @examples
#' X <- analyzeCDNames(df = head(RP, 2), descriptors = CDNames[1:2])
#'
analyzeCDNames <- function(df, descriptors = rcdk::get.desc.names(type = "all")) {
    n <- nrow(df)
    k <- length(descriptors)
    dfs <- list()
    for (j in seq_along(descriptors)) {
        catf("Descriptor %d/%d", j, k)
        desc <- descriptors[j]
        dfs[[j]] <- df[, c()]
        for (i in seq_len(n)) {
            smi <- df$SMILES[i]
            obj <- rcdk::parse.smiles(smi)
            cds <- rcdk::eval.desc(obj, desc, verbose = FALSE)
            dfs[[j]][i, colnames(cds)] <- cds
        }
    }
    data.frame(descriptor = descriptors, all_na = sapply(dfs, function(x) all(is.na(x))))
}

#' @export
#' @keywords internal
#'
#' @title Chemical Descriptor Features
#'
#' @description
#' Vector containing the feature names of the chemical descriptors listed in
#' [CDNames].
#'
#' @seealso [CDNames]
#'
CDFeatures <- c(
    "Fsp3", "nSmallRings", "nAromRings", "nRingBlocks",
    "nAromBlocks", "nRings3", "nRings4", "nRings5", "nRings6", "nRings7",
    "nRings8", "nRings9", "tpsaEfficiency", "Zagreb", "XLogP", "WPATH", "WPOL",
    "WTPT.1", "WTPT.2", "WTPT.3", "WTPT.4", "WTPT.5", "MW", "VAdjMat",
    "TopoPSA", "LipinskiFailures", "nRotB", "topoShape", "geomShape",
    "PetitjeanNumber", "MDEC.11", "MDEC.12", "MDEC.13", "MDEC.14", "MDEC.22",
    "MDEC.23", "MDEC.24", "MDEC.33", "MDEC.34", "MDEC.44", "MDEO.11", "MDEO.12",
    "MDEO.22", "MDEN.11", "MDEN.12", "MDEN.13", "MDEN.22", "MDEN.23", "MDEN.33",
    "MLogP", "nAtomP", "nAtomLC", "khs.sLi", "khs.ssBe", "khs.ssssBe",
    "khs.ssBH", "khs.sssB", "khs.ssssB", "khs.sCH3", "khs.dCH2", "khs.ssCH2",
    "khs.tCH", "khs.dsCH", "khs.aaCH", "khs.sssCH", "khs.ddC", "khs.tsC",
    "khs.dssC", "khs.aasC", "khs.aaaC", "khs.ssssC", "khs.sNH3", "khs.sNH2",
    "khs.ssNH2", "khs.dNH", "khs.ssNH", "khs.aaNH", "khs.tN", "khs.sssNH",
    "khs.dsN", "khs.aaN", "khs.sssN", "khs.ddsN", "khs.aasN", "khs.ssssN",
    "khs.sOH", "khs.dO", "khs.ssO", "khs.aaO", "khs.sF", "khs.sSiH3",
    "khs.ssSiH2", "khs.sssSiH", "khs.ssssSi", "khs.sPH2", "khs.ssPH",
    "khs.sssP", "khs.dsssP", "khs.sssssP", "khs.sSH", "khs.dS", "khs.ssS",
    "khs.aaS", "khs.dssS", "khs.ddssS", "khs.sCl", "khs.sGeH3", "khs.ssGeH2",
    "khs.sssGeH", "khs.ssssGe", "khs.sAsH2", "khs.ssAsH", "khs.sssAs",
    "khs.sssdAs", "khs.sssssAs", "khs.sSeH", "khs.dSe", "khs.ssSe", "khs.aaSe",
    "khs.dssSe", "khs.ddssSe", "khs.sBr", "khs.sSnH3", "khs.ssSnH2",
    "khs.sssSnH", "khs.ssssSn", "khs.sI", "khs.sPbH3", "khs.ssPbH2",
    "khs.sssPbH", "khs.ssssPb", "Kier1", "Kier2", "Kier3", "HybRatio", "nHBDon",
    "nHBAcc", "fragC", "FMF", "ECCEN", "SP.0", "SP.1", "SP.2", "SP.3", "SP.4",
    "SP.5", "SP.6", "SP.7", "VP.0", "VP.1", "VP.2", "VP.3", "VP.4", "VP.5",
    "VP.6", "VP.7", "SPC.4", "SPC.5", "SPC.6", "VPC.4", "VPC.5", "VPC.6",
    "SC.3", "SC.4", "SC.5", "SC.6", "VC.3", "VC.4", "VC.5", "VC.6", "SCH.3",
    "SCH.4", "SCH.5", "SCH.6", "SCH.7", "VCH.3", "VCH.4", "VCH.5", "VCH.6",
    "VCH.7", "C1SP1", "C2SP1", "C1SP2", "C2SP2", "C3SP2", "C1SP3", "C2SP3",
    "C3SP3", "C4SP3", "bpol", "nB", "BCUTw.1l", "BCUTw.1h", "BCUTc.1l",
    "BCUTc.1h", "BCUTp.1l", "BCUTp.1h", "nBase", "ATSp1", "ATSp2", "ATSp3",
    "ATSp4", "ATSp5", "ATSm1", "ATSm2", "ATSm3", "ATSm4", "ATSm5", "ATSc1",
    "ATSc2", "ATSc3", "ATSc4", "ATSc5", "nAtom", "nAromBond", "naAromAtom",
    "apol", "ALogP", "ALogp2", "AMR", "nAcid", "nA", "nR", "nN", "nD", "nC",
    "nF", "nQ", "nE", "nG", "nH", "nI", "nP", "nL", "nK", "nM", "nS", "nT",
    "nY", "nV", "nW"
)

#' @noRd
#' @title All Chemical Descriptor Features available with rcdklibs v2.3
CDFeatures_v2.3 <- c(
    "Fsp3", "nSmallRings", "nAromRings", "nRingBlocks", "nAromBlocks",
    "nRings3", "nRings4", "nRings5", "nRings6", "nRings7", "nRings8", "nRings9",
    "tpsaEfficiency", "Zagreb", "XLogP", "WPATH", "WPOL", "Wlambda1.unity",
    "Wlambda2.unity", "Wlambda3.unity", "Wnu1.unity", "Wnu2.unity",
    "Wgamma1.unity", "Wgamma2.unity", "Wgamma3.unity", "Weta1.unity",
    "Weta2.unity", "Weta3.unity", "WT.unity", "WA.unity", "WV.unity",
    "WK.unity", "WG.unity", "WD.unity", "WTPT.1", "WTPT.2", "WTPT.3", "WTPT.4",
    "WTPT.5", "MW", "VAdjMat", "VABC", "TopoPSA", "LipinskiFailures", "nRotB",
    "topoShape", "geomShape", "PetitjeanNumber", "MOMI.X", "MOMI.Y", "MOMI.Z",
    "MOMI.XY", "MOMI.XZ", "MOMI.YZ", "MOMI.R", "MDEC.11", "MDEC.12", "MDEC.13",
    "MDEC.14", "MDEC.22", "MDEC.23", "MDEC.24", "MDEC.33", "MDEC.34", "MDEC.44",
    "MDEO.11", "MDEO.12", "MDEO.22", "MDEN.11", "MDEN.12", "MDEN.13", "MDEN.22",
    "MDEN.23", "MDEN.33", "MLogP", "nAtomLAC", "LOBMAX", "LOBMIN", "nAtomP",
    "nAtomLC", "khs.sLi", "khs.ssBe", "khs.ssssBe", "khs.ssBH", "khs.sssB",
    "khs.ssssB", "khs.sCH3", "khs.dCH2", "khs.ssCH2", "khs.tCH", "khs.dsCH",
    "khs.aaCH", "khs.sssCH", "khs.ddC", "khs.tsC", "khs.dssC", "khs.aasC",
    "khs.aaaC", "khs.ssssC", "khs.sNH3", "khs.sNH2", "khs.ssNH2", "khs.dNH",
    "khs.ssNH", "khs.aaNH", "khs.tN", "khs.sssNH", "khs.dsN", "khs.aaN",
    "khs.sssN", "khs.ddsN", "khs.aasN", "khs.ssssN", "khs.sOH", "khs.dO",
    "khs.ssO", "khs.aaO", "khs.sF", "khs.sSiH3", "khs.ssSiH2", "khs.sssSiH",
    "khs.ssssSi", "khs.sPH2", "khs.ssPH", "khs.sssP", "khs.dsssP", "khs.sssssP",
    "khs.sSH", "khs.dS", "khs.ssS", "khs.aaS", "khs.dssS", "khs.ddssS",
    "khs.sCl", "khs.sGeH3", "khs.ssGeH2", "khs.sssGeH", "khs.ssssGe",
    "khs.sAsH2", "khs.ssAsH", "khs.sssAs", "khs.sssdAs", "khs.sssssAs",
    "khs.sSeH", "khs.dSe", "khs.ssSe", "khs.aaSe", "khs.dssSe", "khs.ddssSe",
    "khs.sBr", "khs.sSnH3", "khs.ssSnH2", "khs.sssSnH", "khs.ssssSn", "khs.sI",
    "khs.sPbH3", "khs.ssPbH2", "khs.sssPbH", "khs.ssssPb", "Kier1", "Kier2",
    "Kier3", "HybRatio", "nHBDon", "nHBAcc", "GRAV.1", "GRAV.2", "GRAV.3",
    "GRAVH.1", "GRAVH.2", "GRAVH.3", "GRAV.4", "GRAV.5", "GRAV.6", "fragC",
    "FMF", "ECCEN", "PPSA.1", "PPSA.2", "PPSA.3", "PNSA.1", "PNSA.2", "PNSA.3",
    "DPSA.1", "DPSA.2", "DPSA.3", "FPSA.1", "FPSA.2", "FPSA.3", "FNSA.1",
    "FNSA.2", "FNSA.3", "WPSA.1", "WPSA.2", "WPSA.3", "WNSA.1", "WNSA.2",
    "WNSA.3", "RPCG", "RNCG", "RPCS", "RNCS", "THSA", "TPSA", "RHSA", "RPSA",
    "SP.0", "SP.1", "SP.2", "SP.3", "SP.4", "SP.5", "SP.6", "SP.7", "VP.0",
    "VP.1", "VP.2", "VP.3", "VP.4", "VP.5", "VP.6", "VP.7", "SPC.4", "SPC.5",
    "SPC.6", "VPC.4", "VPC.5", "VPC.6", "SC.3", "SC.4", "SC.5", "SC.6", "VC.3",
    "VC.4", "VC.5", "VC.6", "SCH.3", "SCH.4", "SCH.5", "SCH.6", "SCH.7",
    "VCH.3", "VCH.4", "VCH.5", "VCH.6", "VCH.7", "C1SP1", "C2SP1", "C1SP2",
    "C2SP2", "C3SP2", "C1SP3", "C2SP3", "C3SP3", "C4SP3", "bpol", "nB",
    "BCUTw.1l", "BCUTw.1h", "BCUTc.1l", "BCUTc.1h", "BCUTp.1l", "BCUTp.1h",
    "nBase", "ATSp1", "ATSp2", "ATSp3", "ATSp4", "ATSp5", "ATSm1", "ATSm2",
    "ATSm3", "ATSm4", "ATSm5", "ATSc1", "ATSc2", "ATSc3", "ATSc4", "ATSc5",
    "nAtom", "nAromBond", "naAromAtom", "apol", "ALogP", "ALogp2", "AMR",
    "nAcid"
)

#' @noRd
#' @title All Chemical Descriptor Features available with rcdklibs v2.9
CDFeatures_v2.9 <- c(
    "Fsp3", "nSmallRings", "nAromRings", "nRingBlocks", "nAromBlocks",
    "nRings3", "nRings4", "nRings5", "nRings6", "nRings7", "nRings8", "nRings9",
    "tpsaEfficiency", "Zagreb", "XLogP", "WPATH", "WPOL", "Wlambda1.unity",
    "Wlambda2.unity", "Wlambda3.unity", "Wnu1.unity", "Wnu2.unity",
    "Wgamma1.unity", "Wgamma2.unity", "Wgamma3.unity", "Weta1.unity",
    "Weta2.unity", "Weta3.unity", "WT.unity", "WA.unity", "WV.unity",
    "WK.unity", "WG.unity", "WD.unity", "WTPT.1", "WTPT.2", "WTPT.3", "WTPT.4",
    "WTPT.5", "MW", "VAdjMat", "VABC", "TopoPSA", "LipinskiFailures", "nRotB",
    "topoShape", "geomShape", "PetitjeanNumber", "MOMI.X", "MOMI.Y", "MOMI.Z",
    "MOMI.XY", "MOMI.XZ", "MOMI.YZ", "MOMI.R", "MDEC.11", "MDEC.12", "MDEC.13",
    "MDEC.14", "MDEC.22", "MDEC.23", "MDEC.24", "MDEC.33", "MDEC.34", "MDEC.44",
    "MDEO.11", "MDEO.12", "MDEO.22", "MDEN.11", "MDEN.12", "MDEN.13", "MDEN.22",
    "MDEN.23", "MDEN.33", "MLogP", "nAtomLAC", "LOBMAX", "LOBMIN", "nAtomP",
    "nAtomLC", "khs.sLi", "khs.ssBe", "khs.ssssBe", "khs.ssBH", "khs.sssB",
    "khs.ssssB", "khs.sCH3", "khs.dCH2", "khs.ssCH2", "khs.tCH", "khs.dsCH",
    "khs.aaCH", "khs.sssCH", "khs.ddC", "khs.tsC", "khs.dssC", "khs.aasC",
    "khs.aaaC", "khs.ssssC", "khs.sNH3", "khs.sNH2", "khs.ssNH2", "khs.dNH",
    "khs.ssNH", "khs.aaNH", "khs.tN", "khs.sssNH", "khs.dsN", "khs.aaN",
    "khs.sssN", "khs.ddsN", "khs.aasN", "khs.ssssN", "khs.sOH", "khs.dO",
    "khs.ssO", "khs.aaO", "khs.sF", "khs.sSiH3", "khs.ssSiH2", "khs.sssSiH",
    "khs.ssssSi", "khs.sPH2", "khs.ssPH", "khs.sssP", "khs.dsssP", "khs.sssssP",
    "khs.sSH", "khs.dS", "khs.ssS", "khs.aaS", "khs.dssS", "khs.ddssS",
    "khs.sCl", "khs.sGeH3", "khs.ssGeH2", "khs.sssGeH", "khs.ssssGe",
    "khs.sAsH2", "khs.ssAsH", "khs.sssAs", "khs.sssdAs", "khs.sssssAs",
    "khs.sSeH", "khs.dSe", "khs.ssSe", "khs.aaSe", "khs.dssSe", "khs.ddssSe",
    "khs.sBr", "khs.sSnH3", "khs.ssSnH2", "khs.sssSnH", "khs.ssssSn", "khs.sI",
    "khs.sPbH3", "khs.ssPbH2", "khs.sssPbH", "khs.ssssPb", "Kier1", "Kier2",
    "Kier3", "HybRatio", "nHBDon", "nHBAcc", "GRAV.1", "GRAV.2", "GRAV.3",
    "GRAVH.1", "GRAVH.2", "GRAVH.3", "GRAV.4", "GRAV.5", "GRAV.6", "fragC",
    "FMF", "ECCEN", "PPSA.1", "PPSA.2", "PPSA.3", "PNSA.1", "PNSA.2", "PNSA.3",
    "DPSA.1", "DPSA.2", "DPSA.3", "FPSA.1", "FPSA.2", "FPSA.3", "FNSA.1",
    "FNSA.2", "FNSA.3", "WPSA.1", "WPSA.2", "WPSA.3", "WNSA.1", "WNSA.2",
    "WNSA.3", "RPCG", "RNCG", "RPCS", "RNCS", "THSA", "TPSA", "RHSA", "RPSA",
    "SP.0", "SP.1", "SP.2", "SP.3", "SP.4", "SP.5", "SP.6", "SP.7", "VP.0",
    "VP.1", "VP.2", "VP.3", "VP.4", "VP.5", "VP.6", "VP.7", "SPC.4", "SPC.5",
    "SPC.6", "VPC.4", "VPC.5", "VPC.6", "SC.3", "SC.4", "SC.5", "SC.6", "VC.3",
    "VC.4", "VC.5", "VC.6", "SCH.3", "SCH.4", "SCH.5", "SCH.6", "SCH.7",
    "VCH.3", "VCH.4", "VCH.5", "VCH.6", "VCH.7", "C1SP1", "C2SP1", "C1SP2",
    "C2SP2", "C3SP2", "C1SP3", "C2SP3", "C3SP3", "C4SP3", "bpol", "nB",
    "BCUTw.1l", "BCUTw.1h", "BCUTc.1l", "BCUTc.1h", "BCUTp.1l", "BCUTp.1h",
    "nBase", "ATSp1", "ATSp2", "ATSp3", "ATSp4", "ATSp5", "ATSm1", "ATSm2",
    "ATSm3", "ATSm4", "ATSm5", "ATSc1", "ATSc2", "ATSc3", "ATSc4", "ATSc5",
    "nAtom", "nAromBond", "naAromAtom", "apol", "ALogP", "ALogp2", "AMR",
    "nAcid", "TAE0", "TAE1", "TAE2", "TAE3", "TAE4", "TAE5", "TAE6", "TAE7",
    "TAE8", "TAE9", "TAE10", "TAE11", "TAE12", "TAE13", "TAE14", "TAE15",
    "TAE16", "TAE17", "TAE18", "TAE19", "TAE20", "TAE21", "TAE22", "TAE23",
    "TAE24", "TAE25", "TAE26", "TAE27", "TAE28", "TAE29", "TAE30", "TAE31",
    "TAE32", "TAE33", "TAE34", "TAE35", "TAE36", "TAE37", "TAE38", "TAE39",
    "TAE40", "TAE41", "TAE42", "TAE43", "TAE44", "TAE45", "TAE46", "TAE47",
    "TAE48", "TAE49", "TAE50", "TAE51", "TAE52", "TAE53", "TAE54", "TAE55",
    "TAE56", "TAE57", "TAE58", "TAE59", "TAE60", "TAE61", "TAE62", "TAE63",
    "TAE64", "TAE65", "TAE66", "TAE67", "TAE68", "TAE69", "TAE70", "TAE71",
    "TAE72", "TAE73", "TAE74", "TAE75", "TAE76", "TAE77", "TAE78", "TAE79",
    "TAE80", "TAE81", "TAE82", "TAE83", "TAE84", "TAE85", "TAE86", "TAE87",
    "TAE88", "TAE89", "TAE90", "TAE91", "TAE92", "TAE93", "TAE94", "TAE95",
    "TAE96", "TAE97", "TAE98", "TAE99", "TAE100", "TAE101", "TAE102", "TAE103",
    "TAE104", "TAE105", "TAE106", "TAE107", "TAE108", "TAE109", "TAE110",
    "TAE111", "TAE112", "TAE113", "TAE114", "TAE115", "TAE116", "TAE117",
    "TAE118", "TAE119", "TAE120", "TAE121", "TAE122", "TAE123", "TAE124",
    "TAE125", "TAE126", "TAE127", "TAE128", "TAE129", "TAE130", "TAE131",
    "TAE132", "TAE133", "TAE134", "TAE135", "TAE136", "TAE137", "TAE138",
    "TAE139", "TAE140", "TAE141", "TAE142", "TAE143", "TAE144", "TAE145",
    "TAE146", "nA", "nR", "nN", "nD", "nC", "nF", "nQ", "nE", "nG", "nH", "nI",
    "nP", "nL", "nK", "nM", "nS", "nT", "nY", "nV", "nW"
)

#' @noRd
#' @keywords internal
#'
#' @title Wrapper for parallel apply operations
#'
#' @description
#' Helper function that provides a unified interface for parallel processing,
#' automatically choosing between cluster-based parallelization (for multiple
#' workers) and sequential processing (for single worker). This function handles
#' cluster creation, package loading, and cleanup automatically.
#'
#' @param NW
#' Number of workers. If > 1, creates a cluster; if <= 1, uses sequential
#' processing
#'
#' @param ITERABLE
#' Vector or list to apply the function over
#'
#' @param EXPORT
#' Character vector of object names to export to cluster workers
#'
#' @param ENVIR
#' Environment from which to export objects
#'
#' @param BENCHMARK
#' Logical indicating whether to return timing information
#'
#' @param FUN
#' Function to apply to each element of ITERABLE
#'
#' @param ... Additional arguments passed to FUN
#'
#' @details
#' When NW > 1:
#' - Creates a cluster with makeCluster()
#' - Exports specified objects to workers via clusterExport()
#' - Attaches the FastRet namespace to each worker via attach(asNamespace("FastRet"))
#' - Uses parLapply for parallel execution
#' - Automatically cleans up cluster on exit
#' - No stdout output from workers (uses nullfile())
#'
#' When NW <= 1:
#' - Falls back to sequential lapply()
#'
#' This approach ensures that parallel workers have access to all FastRet functions
#' without requiring explicit exports or complex namespace management.
#'
#' @return List containing the results of applying FUN to each element of ITERABLE.
#' If BENCHMARK=TRUE, includes timing attributes.
#'
parLapply2 <- function( NW,
                        ITERABLE,
                        FUN,
                        ...,
                        EXPORT = character(),
                        ENVIR = parent.frame(),
                        BENCHMARK = FALSE,
                        ATTACHPE = TRUE,
                        ATTACHPNS = TRUE) {
    timestamps <- rep(Sys.time(), 6)
    if (NW <= 1) {
        RETOBJ <- lapply(X = ITERABLE, FUN = FUN, ...)
        timestamps[6] <- Sys.time()
    } else {
        # Create Cluster
        cl <- makeCluster(NW)
        on.exit(stopCluster(cl), add = TRUE)
        timestamps[2] <- Sys.time()

        # Object Export
        clusterExport(cl, varlist = EXPORT, envir = ENVIR)
        timestamps[3] <- Sys.time()

        # Attach Package Environment
        if (ATTACHPE || ATTACHPNS) {
            expr <- parse(text = 'library("FastRet")')
            clusterCall(cl, eval, expr)
        }
        timestamps[4] <- Sys.time()

        # Attach Package Namespace
        if (ATTACHPNS) {
            expr <- parse(text = 'attach(asNamespace("FastRet"))')
            clusterCall(cl, eval, expr)
        }
        timestamps[5] <- Sys.time()

        # Function Execution
        RETOBJ <- parLapply(cl = cl, X = ITERABLE, fun = FUN, ...)
        timestamps[6] <- Sys.time()
    }
    if (BENCHMARK) {
        timings <- as.numeric(difftime(
            time1 = timestamps[-1],
            time2 = timestamps[-length(timestamps)],
            units = "secs"
        ))
        names(timings) <- c("init", "export", "pe", "pns", "call")
        attributes(RETOBJ)$timings <- timings
    }
    RETOBJ
}
