# Build a small xlsx of valid SMILES that are guaranteed NOT in the shipped
# descriptor cache, so the GUI's first train/predict on it exercises the slow
# (uncached) rCDK path and a second run is fast (served from CDs.sqlite).
#   Rscript misc/scripts/make-uncached-test-xlsx.R
# Output: misc/uncached-test-molecules.xlsx (NAME, RT, SMILES)

suppressMessages(pkgload::load_all(".", quiet = TRUE))

# Candidate drug-like molecules (not LC-MS metabolites, so unlikely to be cached).
cand <- c(
    Ibuprofen        = "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    Aspirin          = "CC(=O)Oc1ccccc1C(=O)O",
    Paracetamol      = "CC(=O)Nc1ccc(O)cc1",
    Naproxen         = "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    Diclofenac       = "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
    Ketoprofen       = "CC(C(=O)O)c1cccc(c1)C(=O)c1ccccc1",
    Indomethacin     = "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
    Warfarin         = "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
    Probenecid       = "CCCN(CCC)S(=O)(=O)c1ccc(cc1)C(=O)O",
    Furosemide       = "NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl",
    Sulfamethoxazole = "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1",
    Trimethoprim     = "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC",
    Propranolol      = "CC(C)NCC(O)COc1cccc2ccccc12",
    Metoprolol       = "COCCc1ccc(OCC(O)CNC(C)C)cc1",
    Amlodipine       = "CCOC(=O)C1=C(COCCN)NC(C)=C(C1c1ccccc1Cl)C(=O)OC",
    Verapamil        = "COc1ccc(CCN(C)CCCC(C#N)(C(C)C)c2ccc(OC)c(OC)c2)cc1OC",
    Fluoxetine       = "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",
    Sertraline       = "CNC1CCC(c2ccc(Cl)c(Cl)c2)c2ccccc21",
    Diazepam         = "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
    Carbamazepine    = "NC(=O)N1c2ccccc2C=Cc2ccccc21",
    Phenytoin        = "O=C1NC(=O)C(c2ccccc2)(c2ccccc2)N1",
    Lidocaine        = "CCN(CC)CC(=O)Nc1c(C)cccc1C",
    Atorvastatin     = "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
    Losartan         = "CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1",
    Ciprofloxacin    = "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
    Amoxicillin      = "CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O",
    Ampicillin       = "CC1(C)SC2C(NC(=O)C(N)c3ccccc3)C(=O)N2C1C(=O)O",
    Diphenhydramine  = "CN(C)CCOC(c1ccccc1)c1ccccc1",
    Loratadine       = "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc32)CC1",
    Cetirizine       = "OC(=O)COCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1",
    Omeprazole       = "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1",
    Simvastatin      = "CCC(C)(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12C",
    Quinine          = "C=CC1CN2CCC1CC2C(O)c1ccnc2ccc(OC)cc12",
    Estradiol        = "CC12CCC3c4ccc(O)cc4CCC3C1CCC2O",
    Gemfibrozil      = "CC(C)(CCCc1cccc(C)c1C)C(=O)O",
    Tolbutamide      = "CCCCNC(=O)NS(=O)(=O)c1ccc(C)cc1",
    Glibenclamide    = "COc1ccc(Cl)cc1C(=O)NCCc1ccc(S(=O)(=O)NC(=O)NC2CCCCC2)cc1",
    Nifedipine       = "COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]",
    Captopril        = "CC(CS)C(=O)N1CCCC1C(=O)O",
    Enalapril        = "CCOC(=O)C(CCc1ccccc1)NC(C)C(=O)N1CCCC1C(=O)O",
    Ketoconazole     = "CC(=O)N1CCN(c2ccc(OCC3COC(Cn4ccnc4)(c4ccc(Cl)cc4Cl)O3)cc2)CC1",
    Clozapine        = "CN1CCN(C2=Nc3cc(Cl)ccc3Nc3ccccc32)CC1",
    Haloperidol      = "OC1(c2ccc(Cl)cc2)CCN(CCCC(=O)c2ccc(F)cc2)CC1",
    Naproxen2        = "CC(C(=O)O)c1ccc2ccccc2c1",
    Salbutamol       = "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1",
    Theobromine_like = "Cn1cnc2c1c(=O)[nH]c(=O)n2C",
    Atrazine         = "CCNc1nc(Cl)nc(NC(C)C)n1",
    Bisphenol_A      = "CC(C)(c1ccc(O)cc1)c1ccc(O)cc1",
    Triclosan        = "Oc1cc(Cl)ccc1Oc1ccc(Cl)cc1Cl",
    Caffeine_analog  = "CCn1cnc2c1c(=O)n(CC)c(=O)n2C"
)

# Cached SMILES (literal keys in the shipped database).
con <- DBI::dbConnect(RSQLite::SQLite(), pkg_file("cachedata/CDs.sqlite", mustWork = TRUE))
cached <- DBI::dbGetQuery(con, "SELECT smiles FROM cds")$smiles
DBI::dbDisconnect(con)

keep_name <- character(0); keep_smi <- character(0)
for (i in seq_along(cand)) {
    smi <- unname(cand[i])
    ok <- tryCatch({
        obj <- rcdk::parse.smiles(smi)
        !is.null(obj[[1]])
    }, error = function(e) FALSE)
    if (!isTRUE(ok)) next
    canon <- tryCatch(as_canonical(smi), error = function(e) NA_character_)
    if (smi %in% cached || (!is.na(canon) && canon %in% cached)) next # already cached
    keep_name <- c(keep_name, names(cand)[i])
    keep_smi  <- c(keep_smi, smi)
}

n <- length(keep_smi)
set.seed(1)
out <- data.frame(
    NAME = keep_name,
    RT = round(stats::runif(n, 1.5, 14.5), 2),
    SMILES = keep_smi,
    stringsAsFactors = FALSE
)
path <- "misc/uncached-test-molecules.xlsx"
openxlsx::write.xlsx(out, path)
cat(sprintf("Wrote %s with %d uncached molecules (of %d candidates).\n", path, n, length(cand)))
cat(sprintf("Verified: 0 of these %d SMILES are in the shipped cache.\n", n))
