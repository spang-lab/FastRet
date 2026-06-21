# One-off migration: convert the shipped in-memory descriptor cache
# inst/cachedata/CDs.rds into the on-disk SQLite cache inst/cachedata/CDs.sqlite
# used since v1.4.0. Run from the package root with the package loaded:
#   Rscript misc/scripts/migrate-cds-to-sqlite.R
# Going forward, updateCachedCDs() is the canonical regenerator.

suppressMessages(pkgload::load_all(".", quiet = TRUE))

rds <- "inst/cachedata/CDs.rds"
sqlite <- "inst/cachedata/CDs.sqlite"
stopifnot(file.exists(rds))

CDs <- readRDS(rds)
stopifnot(is.data.frame(CDs) || is.matrix(CDs))
CDs <- as.data.frame(CDs)

# Sanity: the curated CDFeatures must be a subset of what rcdklibs >= 2.9 yields
# (the amino-acid-count features nA..nW were added after v2.3; all CDFeatures are
# present from v2.9 on), and of the shipped cache's columns.
stopifnot(all(CDFeatures %in% CDFeatures_v2.9))
stopifnot(all(CDFeatures %in% colnames(CDs)))
stopifnot(!is.null(rownames(CDs)))

df <- data.frame(smiles = rownames(CDs), CDs[, CDFeatures, drop = FALSE],
                 check.names = FALSE, row.names = NULL, stringsAsFactors = FALSE)

if (file.exists(sqlite)) unlink(sqlite)
con <- DBI::dbConnect(RSQLite::SQLite(), sqlite)
on.exit(DBI::dbDisconnect(con), add = TRUE)
cols <- paste(sprintf('"%s" REAL', CDFeatures), collapse = ", ")
DBI::dbExecute(con, sprintf("CREATE TABLE cds (smiles TEXT PRIMARY KEY, %s);", cols))
DBI::dbWriteTable(con, "cds", df, append = TRUE)
DBI::dbExecute(con, "PRAGMA wal_checkpoint(TRUNCATE);") # ensure a single self-contained file

n <- DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM cds")$n
cat(sprintf("Wrote %s: %d rows, %d descriptor columns\n", sqlite, n, length(CDFeatures)))
stopifnot(n == nrow(df))
