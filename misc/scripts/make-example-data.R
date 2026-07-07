# Regenerate the example datasets/models shipped with FastRet from the published
# Measurements_v10P.xlsx. Run from the package root after devtools::load_all().
#
#   Rscript misc/scripts/make-example-data.R
#
# Produces (all derived from DSID subsets of Measurements_v10P.xlsx, raw SMILES):
#   - inst/extdata/RP.xlsx           RP dataset, 458 x 4 (RT, SMILES, NAME, INCHIKEY)
#   - data/RP.rda                    lazy-loaded `RP` (== read_rp_xlsx())
#   - inst/extdata/RP_adj.xlsx       RP_Steep, 25 x 4 (RT, NAME, SMILES, INCHIKEY)
#   - inst/extdata/RP_lasso_model.rds  LASSO frm trained on the new RP
#
# Requires Java + rcdk. Rebuild inst/cachedata/CDs.sqlite first (updateCachedCDs())
# so training hits the cache; a fresh per-user cache dir is used here so it seeds
# from the freshly built shipped database.

devtools::load_all(".")

# Use a fresh writable cache so training seeds from the just-rebuilt shipped
# CDs.sqlite instead of a possibly stale per-user copy.
options(FastRet.cache_dir = file.path(tempdir(), "FastRet-make-example-cache"))

url <- "https://github.com/spang-lab/FastRet/releases/download/example-data/Measurements_v10P.xlsx"
v10p_path <- tempfile("Measurements_v10P", fileext = ".xlsx")
download.file(url, v10p_path, mode = "wb")
X <- openxlsx::read.xlsx(v10p_path, 1)

# 1. RP.xlsx (RT, SMILES, NAME, INCHIKEY) -- 458 rows
RP_new <- X[X$DSID == "RP", c("RT", "SMILES", "NAME", "INCHIKEY")]
RP_new <- RP_new[order(RP_new$RT), ]
rownames(RP_new) <- NULL
openxlsx::write.xlsx(RP_new, file.path("inst", "extdata", "RP.xlsx"), rowNames = FALSE)

# 2. data/RP.rda -- read back so RP == read_rp_xlsx()
RP <- read_rp_xlsx()
usethis::use_data(RP, overwrite = TRUE)

# 3. RP_adj.xlsx (RT, NAME, SMILES, INCHIKEY) -- 25 real steeper-gradient re-measurements
RP_adj <- X[X$DSID == "RP_Steep", c("RT", "NAME", "SMILES", "INCHIKEY")]
rownames(RP_adj) <- NULL
openxlsx::write.xlsx(RP_adj, file.path("inst", "extdata", "RP_adj.xlsx"), rowNames = FALSE)

# 4. RP_lasso_model.rds -- LASSO trained on the new RP
frm <- train_frm(RP, method = "lasso", verbose = 1, seed = 42)
saveRDS(frm, file.path("inst", "extdata", "RP_lasso_model.rds"))

# --- Sanity checks -----------------------------------------------------------
stopifnot(
    isTRUE(all.equal(read_rp_xlsx(), RP)),
    identical(dim(read_rpadj_xlsx()), c(25L, 4L)),
    nrow(RP) == 458L,
    all(read_rpadj_xlsx()$SMILES %in% RP$SMILES),
    all(read_rpadj_xlsx()$INCHIKEY %in% RP$INCHIKEY),
    nrow(frm$df) == 458L,
    "INCHIKEY" %in% colnames(frm$df)
)
cat("make-example-data.R: OK\n")
