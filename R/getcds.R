# Main #####

#' @title Get Chemical Descriptors for a list of molecules
#' @description Calculate Chemical Descriptors for a list of molecules. Molecules can appear multiple times in the list.
#' @param df dataframe with two mandatory columns: "NAME" and "SMILES"
#' @param verbose 0: no output, 1: progress, 2: more progress and warnings
#' @param nw number of workers for parallel processing
#' @keywords public
#' @export
getCDs <- function(df = read_rp_xlsx(), verbose = 1, nw = 1) {

  # Configure logging behaviour
  if (verbose == 0) {
    catf <- function(...) invisible() # disable catf prints
  }
  if (verbose >= 2) {
    opts <- options(warn = 1)
    on.exit(options(opts), add = TRUE)
    suppressWarnings <- function(expr) expr # dont suppress warnings later on
  }

  # Return pregenerated results when mocking is enabled (e.g. while testing the GUI)
  if ("getCDs" %in% getOption("FastRet.mocks", c())) {
    catf("Mocking is enabled for 'getCDs'. Returning 'mockdata/RPCD.rds'.")
    return(readRDS(pkg_file("mockdata/RPCD.rds")))
  }

  catf("Obtaining chemical descriptors for a dataframe of dimension %d x %d", nrow(df), ncol(df))
  a <- Sys.time()
  if (nw > 1) {
    catf("Creating cluster for parallel processing")
    # Parallelization using forking is not possible. Probably because the `getCDsFor1Molecule` calls `rcdk::eval.desc` which calls java in the background which is probably not thread-safe.
    cl <- makeCluster(nw, outfile = if (verbose >= 2) "" else nullfile())
    clusterExport(cl, c("getCDsFor1Molecule", "get_cache_dir", "CDNames"))
    on.exit(stopCluster(cl), add = TRUE)
    catf("Calculating chemical descriptors using %d workers", nw)
    cds <- parLapply(cl, df$SMILES, getCDsFor1Molecule, verbose = verbose - 1)
    catf("Collecting results")
    cds <- do.call(rbind, cds)
  } else {
    catf("Calculating chemical descriptors")
    cds <- lapply(df$SMILES, getCDsFor1Molecule, verbose = verbose - 1)
    cds <- do.call(rbind, cds)
  }
  retdf <- cbind(df, cds)
  b <- Sys.time()
  catf("Finished calculating chemical descriptors in %s", format(b - a))
  invisible(retdf)
}

getCDsCache <- as.environment(list())

#' @title Get Chemical Descriptors for a single molecule
#' @description Helper function for [getCDs()]. Calculates chemical descriptors for a single molecule, specified as [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) string. This function should NOT be used directly. It is only exported so [getCDs()] can spawn background worker processes that are able to call this function.
#' @details Chemical descriptors in [getCDs()] are calculated individually for each molecule. This is due to the inconsistent ordering of output dataframes when a list of IAtomContainers is provided to `rcdk::eval.desc`. Although the input SMILES are set as rownames, they don't match the original input SMILES due to an unclear transformation, making mapping non-trivial. Calculating descriptors molecule by molecule also enables parallelization in [getCDs()].
#' @param smi SMILES string of the molecule
#' @param cache if TRUE, the results are cached in a directory `~/.cache/FastRet/getCDsFor1Molecule/` to speed up subsequent calls
#' @param verbose 0: no output, 1: show progress
#' @keywords internal
#' @export
getCDsFor1Molecule <- function(smi = "O=C(O)CCCCCCCCCO", cache = TRUE, verbose = 1) {
  if (verbose <- 0) catf <- function(...) invisible() # disable catf prints
  if (cache) {
    cache_dir <- get_cache_dir("getCDsFor1Molecule")
    cache_file <- paste0(digest::digest(smi), ".rds")
    cache_path <- file.path(cache_dir, cache_file)
    if (file.exists(cache_path)) {
      catf("Loading chemical descriptors for '%s' from cachefile '%s'", smi, cache_path)
      cds <- readRDS(cache_path)
    } else {
      cds <- getCDsFor1Molecule(smi, cache = FALSE)
      catf("Caching chemical descriptors for '%s' at '%s'", smi, cache_path)
      saveRDS(cds, cache_path)
    }
  } else {
    # To understand the following part, a rough understanding of rcdk is required. Therefore, it is recommended to read the following vignette: https://cran.r-project.org/web/packages/rcdk/vignettes/using-rcdk.html#molecular-descriptors
    catf("Calculating chemical descriptors for: '%s'", smi)
    obj <- rcdk::parse.smiles(smi)[[1]] # Externalptr to Java IAtomContainer objects
    rcdk::convert.implicit.to.explicit(obj) # Add hydrogen atoms if omitted in the SMILES
    rcdk::generate.2d.coordinates(obj) # Generate 2D coordinates
    cds <- suppressWarnings(rcdk::eval.desc(obj, CDNames, verbose = FALSE))
  }
  cds
}


# Constants #####

#' @title Analyze Chemical Descriptors Names
#' @description Analyze the chemical descriptors names and return a dataframe with the names and a boolean column indicating if all values are NA.
#' @details This function is used to analyze the chemical descriptors names and to identify which descriptors produce only NAs in the test datasets. The function is used to generate the CDNames object.
#' @keywords internal
#' @export
analyzeCDNames <- function() {
  df <- read_rp_xlsx()
  n <- nrow(df)
  descriptors <- rcdk::get.desc.names(type = "all")
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

#' @title Chemical Descriptors Names
#' @description This object contains the names of various chemical descriptors.
#' @details One descriptor can be associated with multiple features, e.g. the BCUT descriptor corresponds to the following features: BCUTw.1l, BCUTw.1h, BCUTc.1l, BCUTc.1h, BCUTp.1l, BCUTp.1h. Some descriptors produce warnings for certain molecules., e.g. "The AtomType null could not be found" or "Molecule must have 3D coordinates" and return NA in such cases. Descriptors that produce only NAs in our test datasets will be excluded. To see which descriptors produce only NAs, run `analyzeCDNames`. The "LongestAliphaticChain" descriptors sometimes even produces `Error: segfault from C stack overflow` error`, e.g. for SMILES `c1ccccc1C(Cl)(Cl)Cl` (== `rcdk::bpdata$SMILES[200]`) when using `OpenJDK Runtime Environment (build 11.0.23+9-post-Ubuntu-1ubuntu122.04.1)`.
#' @examples \dontrun{
#' CDNames <- rcdk::get.desc.names(type = "all")
#' skipPattern <- "(WHIM|VABC|MomentOfInertia|LengthOverBreadth|GravitationalIndex|CPSA|TaeAminoAcid|LongestAliphaticChain)"
#' CDNames <- CDNames[!grepl(skipPattern, CDNames)]
#' }
#' @seealso [analyzeCDNames()], [CDs]
#' @keywords internal
#' @export
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
CDs_v2.3 <- c("Fsp3", "nSmallRings", "nAromRings", "nRingBlocks", "nAromBlocks", "nRings3", "nRings4", "nRings5", "nRings6", "nRings7", "nRings8", "nRings9", "tpsaEfficiency", "Zagreb", "XLogP", "WPATH", "WPOL", "Wlambda1.unity", "Wlambda2.unity", "Wlambda3.unity", "Wnu1.unity", "Wnu2.unity", "Wgamma1.unity", "Wgamma2.unity", "Wgamma3.unity", "Weta1.unity", "Weta2.unity", "Weta3.unity", "WT.unity", "WA.unity", "WV.unity", "WK.unity", "WG.unity", "WD.unity", "WTPT.1", "WTPT.2", "WTPT.3", "WTPT.4", "WTPT.5", "MW", "VAdjMat", "VABC", "TopoPSA", "LipinskiFailures", "nRotB", "topoShape", "geomShape", "PetitjeanNumber", "MOMI.X", "MOMI.Y", "MOMI.Z", "MOMI.XY", "MOMI.XZ", "MOMI.YZ", "MOMI.R", "MDEC.11", "MDEC.12", "MDEC.13", "MDEC.14", "MDEC.22", "MDEC.23", "MDEC.24", "MDEC.33", "MDEC.34", "MDEC.44", "MDEO.11", "MDEO.12", "MDEO.22", "MDEN.11", "MDEN.12", "MDEN.13", "MDEN.22", "MDEN.23", "MDEN.33", "MLogP", "nAtomLAC", "LOBMAX", "LOBMIN", "nAtomP", "nAtomLC", "khs.sLi", "khs.ssBe", "khs.ssssBe", "khs.ssBH", "khs.sssB", "khs.ssssB", "khs.sCH3", "khs.dCH2", "khs.ssCH2", "khs.tCH", "khs.dsCH", "khs.aaCH", "khs.sssCH", "khs.ddC", "khs.tsC", "khs.dssC", "khs.aasC", "khs.aaaC", "khs.ssssC", "khs.sNH3", "khs.sNH2", "khs.ssNH2", "khs.dNH", "khs.ssNH", "khs.aaNH", "khs.tN", "khs.sssNH", "khs.dsN", "khs.aaN", "khs.sssN", "khs.ddsN", "khs.aasN", "khs.ssssN", "khs.sOH", "khs.dO", "khs.ssO", "khs.aaO", "khs.sF", "khs.sSiH3", "khs.ssSiH2", "khs.sssSiH", "khs.ssssSi", "khs.sPH2", "khs.ssPH", "khs.sssP", "khs.dsssP", "khs.sssssP", "khs.sSH", "khs.dS", "khs.ssS", "khs.aaS", "khs.dssS", "khs.ddssS", "khs.sCl", "khs.sGeH3", "khs.ssGeH2", "khs.sssGeH", "khs.ssssGe", "khs.sAsH2", "khs.ssAsH", "khs.sssAs", "khs.sssdAs", "khs.sssssAs", "khs.sSeH", "khs.dSe", "khs.ssSe", "khs.aaSe", "khs.dssSe", "khs.ddssSe", "khs.sBr", "khs.sSnH3", "khs.ssSnH2", "khs.sssSnH", "khs.ssssSn", "khs.sI", "khs.sPbH3", "khs.ssPbH2", "khs.sssPbH", "khs.ssssPb", "Kier1", "Kier2", "Kier3", "HybRatio", "nHBDon", "nHBAcc", "GRAV.1", "GRAV.2", "GRAV.3", "GRAVH.1", "GRAVH.2", "GRAVH.3", "GRAV.4", "GRAV.5", "GRAV.6", "fragC", "FMF", "ECCEN", "PPSA.1", "PPSA.2", "PPSA.3", "PNSA.1", "PNSA.2", "PNSA.3", "DPSA.1", "DPSA.2", "DPSA.3", "FPSA.1", "FPSA.2", "FPSA.3", "FNSA.1", "FNSA.2", "FNSA.3", "WPSA.1", "WPSA.2", "WPSA.3", "WNSA.1", "WNSA.2", "WNSA.3", "RPCG", "RNCG", "RPCS", "RNCS", "THSA", "TPSA", "RHSA", "RPSA", "SP.0", "SP.1", "SP.2", "SP.3", "SP.4", "SP.5", "SP.6", "SP.7", "VP.0", "VP.1", "VP.2", "VP.3", "VP.4", "VP.5", "VP.6", "VP.7", "SPC.4", "SPC.5", "SPC.6", "VPC.4", "VPC.5", "VPC.6", "SC.3", "SC.4", "SC.5", "SC.6", "VC.3", "VC.4", "VC.5", "VC.6", "SCH.3", "SCH.4", "SCH.5", "SCH.6", "SCH.7", "VCH.3", "VCH.4", "VCH.5", "VCH.6", "VCH.7", "C1SP1", "C2SP1", "C1SP2", "C2SP2", "C3SP2", "C1SP3", "C2SP3", "C3SP3", "C4SP3", "bpol", "nB", "BCUTw.1l", "BCUTw.1h", "BCUTc.1l", "BCUTc.1h", "BCUTp.1l", "BCUTp.1h", "nBase", "ATSp1", "ATSp2", "ATSp3", "ATSp4", "ATSp5", "ATSm1", "ATSm2", "ATSm3", "ATSm4", "ATSm5", "ATSc1", "ATSc2", "ATSc3", "ATSc4", "ATSc5", "nAtom", "nAromBond", "naAromAtom", "apol", "ALogP", "ALogp2", "AMR", "nAcid")
CDs_v2.9 <- c("Fsp3", "nSmallRings", "nAromRings", "nRingBlocks", "nAromBlocks", "nRings3", "nRings4", "nRings5", "nRings6", "nRings7", "nRings8", "nRings9", "tpsaEfficiency", "Zagreb", "XLogP", "WPATH", "WPOL", "Wlambda1.unity", "Wlambda2.unity", "Wlambda3.unity", "Wnu1.unity", "Wnu2.unity", "Wgamma1.unity", "Wgamma2.unity", "Wgamma3.unity", "Weta1.unity", "Weta2.unity", "Weta3.unity", "WT.unity", "WA.unity", "WV.unity", "WK.unity", "WG.unity", "WD.unity", "WTPT.1", "WTPT.2", "WTPT.3", "WTPT.4", "WTPT.5", "MW", "VAdjMat", "VABC", "TopoPSA", "LipinskiFailures", "nRotB", "topoShape", "geomShape", "PetitjeanNumber", "MOMI.X", "MOMI.Y", "MOMI.Z", "MOMI.XY", "MOMI.XZ", "MOMI.YZ", "MOMI.R", "MDEC.11", "MDEC.12", "MDEC.13", "MDEC.14", "MDEC.22", "MDEC.23", "MDEC.24", "MDEC.33", "MDEC.34", "MDEC.44", "MDEO.11", "MDEO.12", "MDEO.22", "MDEN.11", "MDEN.12", "MDEN.13", "MDEN.22", "MDEN.23", "MDEN.33", "MLogP", "nAtomLAC", "LOBMAX", "LOBMIN", "nAtomP", "nAtomLC", "khs.sLi", "khs.ssBe", "khs.ssssBe", "khs.ssBH", "khs.sssB", "khs.ssssB", "khs.sCH3", "khs.dCH2", "khs.ssCH2", "khs.tCH", "khs.dsCH", "khs.aaCH", "khs.sssCH", "khs.ddC", "khs.tsC", "khs.dssC", "khs.aasC", "khs.aaaC", "khs.ssssC", "khs.sNH3", "khs.sNH2", "khs.ssNH2", "khs.dNH", "khs.ssNH", "khs.aaNH", "khs.tN", "khs.sssNH", "khs.dsN", "khs.aaN", "khs.sssN", "khs.ddsN", "khs.aasN", "khs.ssssN", "khs.sOH", "khs.dO", "khs.ssO", "khs.aaO", "khs.sF", "khs.sSiH3", "khs.ssSiH2", "khs.sssSiH", "khs.ssssSi", "khs.sPH2", "khs.ssPH", "khs.sssP", "khs.dsssP", "khs.sssssP", "khs.sSH", "khs.dS", "khs.ssS", "khs.aaS", "khs.dssS", "khs.ddssS", "khs.sCl", "khs.sGeH3", "khs.ssGeH2", "khs.sssGeH", "khs.ssssGe", "khs.sAsH2", "khs.ssAsH", "khs.sssAs", "khs.sssdAs", "khs.sssssAs", "khs.sSeH", "khs.dSe", "khs.ssSe", "khs.aaSe", "khs.dssSe", "khs.ddssSe", "khs.sBr", "khs.sSnH3", "khs.ssSnH2", "khs.sssSnH", "khs.ssssSn", "khs.sI", "khs.sPbH3", "khs.ssPbH2", "khs.sssPbH", "khs.ssssPb", "Kier1", "Kier2", "Kier3", "HybRatio", "nHBDon", "nHBAcc", "GRAV.1", "GRAV.2", "GRAV.3", "GRAVH.1", "GRAVH.2", "GRAVH.3", "GRAV.4", "GRAV.5", "GRAV.6", "fragC", "FMF", "ECCEN", "PPSA.1", "PPSA.2", "PPSA.3", "PNSA.1", "PNSA.2", "PNSA.3", "DPSA.1", "DPSA.2", "DPSA.3", "FPSA.1", "FPSA.2", "FPSA.3", "FNSA.1", "FNSA.2", "FNSA.3", "WPSA.1", "WPSA.2", "WPSA.3", "WNSA.1", "WNSA.2", "WNSA.3", "RPCG", "RNCG", "RPCS", "RNCS", "THSA", "TPSA", "RHSA", "RPSA", "SP.0", "SP.1", "SP.2", "SP.3", "SP.4", "SP.5", "SP.6", "SP.7", "VP.0", "VP.1", "VP.2", "VP.3", "VP.4", "VP.5", "VP.6", "VP.7", "SPC.4", "SPC.5", "SPC.6", "VPC.4", "VPC.5", "VPC.6", "SC.3", "SC.4", "SC.5", "SC.6", "VC.3", "VC.4", "VC.5", "VC.6", "SCH.3", "SCH.4", "SCH.5", "SCH.6", "SCH.7", "VCH.3", "VCH.4", "VCH.5", "VCH.6", "VCH.7", "C1SP1", "C2SP1", "C1SP2", "C2SP2", "C3SP2", "C1SP3", "C2SP3", "C3SP3", "C4SP3", "bpol", "nB", "BCUTw.1l", "BCUTw.1h", "BCUTc.1l", "BCUTc.1h", "BCUTp.1l", "BCUTp.1h", "nBase", "ATSp1", "ATSp2", "ATSp3", "ATSp4", "ATSp5", "ATSm1", "ATSm2", "ATSm3", "ATSm4", "ATSm5", "ATSc1", "ATSc2", "ATSc3", "ATSc4", "ATSc5", "nAtom", "nAromBond", "naAromAtom", "apol", "ALogP", "ALogp2", "AMR", "nAcid", "TAE0", "TAE1", "TAE2", "TAE3", "TAE4", "TAE5", "TAE6", "TAE7", "TAE8", "TAE9", "TAE10", "TAE11", "TAE12", "TAE13", "TAE14", "TAE15", "TAE16", "TAE17", "TAE18", "TAE19", "TAE20", "TAE21", "TAE22", "TAE23", "TAE24", "TAE25", "TAE26", "TAE27", "TAE28", "TAE29", "TAE30", "TAE31", "TAE32", "TAE33", "TAE34", "TAE35", "TAE36", "TAE37", "TAE38",  "TAE39", "TAE40", "TAE41", "TAE42", "TAE43", "TAE44", "TAE45", "TAE46", "TAE47", "TAE48", "TAE49", "TAE50", "TAE51", "TAE52", "TAE53", "TAE54", "TAE55", "TAE56", "TAE57", "TAE58", "TAE59", "TAE60", "TAE61", "TAE62", "TAE63", "TAE64", "TAE65", "TAE66", "TAE67", "TAE68", "TAE69", "TAE70", "TAE71", "TAE72", "TAE73", "TAE74", "TAE75", "TAE76", "TAE77", "TAE78", "TAE79", "TAE80", "TAE81", "TAE82", "TAE83", "TAE84", "TAE85", "TAE86", "TAE87", "TAE88", "TAE89", "TAE90", "TAE91", "TAE92", "TAE93", "TAE94", "TAE95", "TAE96", "TAE97", "TAE98", "TAE99", "TAE100", "TAE101", "TAE102", "TAE103", "TAE104", "TAE105", "TAE106", "TAE107", "TAE108", "TAE109", "TAE110", "TAE111", "TAE112", "TAE113", "TAE114", "TAE115", "TAE116", "TAE117", "TAE118", "TAE119", "TAE120", "TAE121", "TAE122", "TAE123", "TAE124", "TAE125", "TAE126", "TAE127", "TAE128", "TAE129", "TAE130", "TAE131", "TAE132", "TAE133", "TAE134", "TAE135", "TAE136", "TAE137", "TAE138", "TAE139", "TAE140", "TAE141", "TAE142", "TAE143", "TAE144", "TAE145", "TAE146", "nA", "nR", "nN", "nD", "nC", "nF", "nQ", "nE", "nG", "nH", "nI", "nP", "nL", "nK", "nM", "nS", "nT", "nY", "nV", "nW")

#' @title Chemical Descriptors
#' @description This object contains the feature names of various chemical descriptors.
#' @keywords internal
#' @seealso [CDNames]
#' @export
CDs <- CDs_v2.9
