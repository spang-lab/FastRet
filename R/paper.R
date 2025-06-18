# Main #####

reproduce_paper <- function() {
    dat <- get_norm_ds()

    # VALIDATE: number of metabolites in in-house RT library (~380)
    pdf.open("out/heatmap_met_cc.pdf")
    try(plot_met_cc_heatmap(dat))
    dev.off()

    pdf.open("out/venn_met_cc.pdf")
    try({
        plot_met_cc_venn(dat, HILIC = FALSE, HILIC_Retip = FALSE)
        plot_met_cc_venn(dat, HILIC = TRUE, HILIC_Retip = FALSE)
        plot_met_cc_venn(dat, HILIC = TRUE, HILIC_Retip = TRUE)
    })
    dev.off()

    # VALIDATE: model adjustment works with as few as 25 metabolites
    # VALIDATE: BRT MAE for RP (0.5), HILIC (0.86), RP-AX (X.XX)
    # VALIDATE: min-size for metabolites for adjustment
    # VALIDATE: number of metabolites in each dataset (Retip HILIC: 970, In-house HILIC: 364, Mixed mode: 412, RP: 40)
    # VALIDATE: RMSE, MAE, R2, PPB1M for all datasets and methods
    # VALIDATE: 200 candidate lambda values, 10-fold CV for Lasso
    # VALIDATE: xgboost hyperparameters (tree depth 2-5, nrounds 300-1000 step 100, eta 0.01/0.02), 10-fold CV
    # VALIDATE: selective measuring method as described is implemented and used (steps 1-5)
    # VALIDATE: k = 25 is used for medoid selection
    # VALIDATE: 25 medoids are used for further validation
    # VALIDATE: adjustment is performed using polynomial regression (linear, quadratic, cubic) on medoids
    # VALIDATE: 10-fold cross-validation is used for all performance evaluations
    # VALIDATE: BRTs outperform Lasso by 0.26 RMSE and 0.24 MAE on average
    # VALIDATE: PPB1M for RP (85.95%), RP-AX (46.37%)
    # VALIDATE: HILIC in-house (364 metabolites, 327 per fold), Retip (970, 873 per fold), in-house HILIC performs slightly better
    # VALIDATE: Figure 1 shows RMSE vs training set size for RP, saturation point
    # VALIDATE: Lasso saturates at ~25, BRTs at ~60 training samples
    # VALIDATE: For small sets Lasso > BRTs, for large sets BRTs > Lasso
    # VALIDATE: external validation set contains 25 metabolites not in training
    # VALIDATE: Figure 2 shows RT changes for all tested conditions
    # VALIDATE: original model MAE = 1 min on validation set
    # VALIDATE: adjusted model MAE = 0.65 min
    # VALIDATE: medoids-only model MAE = 1.06 min, original model trained on 401 metabolites
    # VALIDATE: adjusted models outperform medoids-only and original model
    # VALIDATE: Table values match actual MAE for each approach
    # VALIDATE: performance comparable to Retip paper
    # VALIDATE: in-house HILIC slightly better than Retip HILIC
    # VALIDATE: BRTs > Lasso for large datasets
    # VALIDATE: Lasso reliable for small training sets
    # VALIDATE: chemical descriptor quality is limiting factor (if possible to check)
}

# Saving #####

pdf.open <- function(path, width = 3.33, height = 3.33) {
    catf("Saving %s", path)
    if (!dir.exists("out")) dir.create("out", recursive = TRUE)
    pdf(path, width = width, height = height)
}

#' @noRd
#' @title Save a data.frame as an interactive standalone HTML file
#' @description
#' This function creates an interactive, sortable, filterable HTML table from
#' a data.frame and saves it as a standalone HTML file using the DT and
#' htmlwidgets packages.
#' @param x A data.frame or tibble to be rendered.
#' @param file Character. Output HTML file path (e.g. "table.html").
#' @param title Character. Optional HTML page title.
#' @param caption Character. Optional table caption.
#' @param pageLength Integer. Number of rows to show per page (default is 25).
#' @param row.names Logical. Whether to include row names. Default is FALSE.
#' @return Invisibly returns the DT widget object.
#' @examples
#' save_as_html(iris, tempfile(fileext = ".html"))
save_as_html <- function(x,
                         file,
                         digits = 2,
                         title = "Interactive Table",
                         caption = NULL,
                         pageLength = 999,
                         row.names = FALSE) {
    kbl <- knitr::kable(x, "html", digits, row.names, caption = caption)
    kbl <- kableExtra::kable_styling(kbl, fixed_thead = TRUE)
    maxchar <- apply(x, 2, function(col) max(nchar(col)))
    longcols <- which(maxchar > 30)
    for (col in longcols) kbl <- kableExtra::column_spec(
        kbl, col, extra_css = "word-break: break-word;"
    )
    html_page <- paste0(
        "<!DOCTYPE html>\n<html><head><meta charset='utf-8'>\n",
        "<title>", title, "</title>\n",
        "<style>",
        "body {font-family: sans-serif; padding: 2em;}",
        "thead th { background: #ccc; }",
        "table {table-layout: auto; width: 100%; border-collapse: collapse;}",
        "td, th {padding: 0.4em; border: 1px solid #ccc; vertical-align: top;}",
        "</style>\n",
        "</head><body>\n",
        as.character(kbl),
        "\n</body></html>"
    )
    writeLines(html_page, con = file)
    invisible(TRUE)
}

save_as_xlsx <- function(x, file, rowNames = FALSE) {
    openxlsx::write.xlsx(x, file, rowNames = rowNames)
}

# Analysis Helpers #####

#' @noRd
#' @title Get duplicates in a vector
#' @param x A vector to check for duplicates.
#' @return A list of indices for each duplicate value in the vector.
#' @examples
#' x <- c("a", "b", "c", "a", "d", "a", "b", "e")
#' get_duplicates(x)
#' ## c(0, 0, 0, 1, 0, 1, 2, 0)
#' ##   |        |     |  |
#' ##   |        |     |  duplicate of element 2 ("b")
#' ##   |        |     another duplicate of element 1 ("a")
#' ##   |        duplicate of element 1 ("a")
#' ##   no duplicate
get_duplicates <- function(x) {
    isdup <- duplicated(x)
    dupid <- which(isdup) # c(4, 6, 7)
    dupof <- vapply(
        dupid,
        function(i) which(x == x[i])[1],
        integer(1)
    ) # c(1, 1, 2)
    df <- data.frame(x = x, isdup = isdup)
    df$dupof <- 0
    df$dupof[dupid] <- dupof
    df$dupof
}

get_dupof <- function(x) {
    namdup <- get_duplicates(x$NAME)
    smidup <- get_duplicates(x$CANONICAL)
    smidup[smidup == namdup] <- 0
    gsub("\\b0\\b", "", paste(namdup, smidup))
}

get_canonical_smiles <- function(smiles) {
    molecules <- rcdk::parse.smiles(smiles)
    flavor <- rcdk::smiles.flavors("Canonical")
    sapply(molecules, rcdk::get.smiles, flavor)
}

# Runtime Helpers #####

#' @noRd
#' @title Measure runtime of an expression.
#' @description
#' Measures runtime the same way as [system.time()] does, but only prints the
#' elapsed time.
#' @param expr An R expression to measure.
#' @return The runtime in seconds.
#' @examples
#' rt(Sys.sleep(1))
rt <- function(expr) {
    gc(verbose = FALSE)
    time <- proc.time()
    expr
    new.time <- proc.time()
    unname(new.time - time)[3]
}

readRdsVerbose <- function(path) {
    catf("Reading '%s'", path)
    readRDS(path)
}

saveRdsVerbose <- function(x, path) {
    catf("Writing '%s'", path)
    saveRDS(x, path)
    catf("Done")
}

# Formatting Helpers #####

wrap_n <- function(x, n = 3) {
  vapply(x, function(s) {
    gsub(paste0("(.{", n, "})"), "\\1\n", s)
  }, character(1))
}

#' @noRd
#' @title Pretty Print Data Frame
pprint <- function(x, ...,
                   digits = 3, # base default is NULL
                   quote = FALSE, # base default is FALSE
                   right = FALSE, # base default is TRUE
                   row.names = TRUE, # base default is TRUE
                   max = NULL, # base default is NULL
                   maxrow = 10, # new argument
                   maxchar = 32) { # new argument
    nrows_total <- nrow(x)
    xhead <- head(x, maxrow)
    nrows <- nrow(xhead)
    for (k in colnames(xhead)) {
        v <- xhead[[k]]
        if (is.character(v)) {
            is_long <- which(nchar(v) > maxchar)
            if (length(is_long) == 0) next
            longs <- v[is_long]
            short <- paste0(substr(longs, 1, 25), "[...]")
            v[is_long] <- short
            xhead[[k]] <- v
        }
        if (is.list(v)) {
            xhead[[k]] <- NULL
            xhead[[k]] <- doc("SUBLIST", nrows)
        }
    }
    print(
        xhead,
        digits = digits,
        quote = quote,
        right = right,
        row.names = row.names,
        max = max
    )
    if (nrows < nrows_total && interactive() && isatty(stdout())) {
        cat(sprintf("\nOmitted rows %d-%d\n", nrows + 1, nrows_total))
        cat("Set 'maxrow' and/or 'maxchar' to configure printing\n\n")
    }
}

# Plotting Helpers #####

rect1 <- function(x, y, col) {
    rect(
        xleft = x - 0.5, ybottom = y - 0.5,
        xright = x + 0.5, ytop = y + 0.5,
        col = col, border = NULL, density = NA
    )
}

# Analysis #####

get_ds <- function() {
    list(
        HILIC = FastRet::HILIC,
        HILIC_Retip = FastRet::read_retip_hilic_data(),
        RP_AXMM = FastRet::RP_AXMM,
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
        RP_Val_T25_FR25_Steep = FastRet::RP_Val$T25_FR25_Steep
    )
}

#' @noRd
#' @title Get Normalized Datasets
#' @details Runtime: 1.26 seconds
#' @examples
#' get_norm_ds(read_cache = FALSE, update_cache = TRUE) # update cache
#' ds <- get_norm_ds()
#' str(ds, 1)
get_norm_ds <- function(read_cache = TRUE, update_cache = TRUE) {
    cachefile <- file.path(get_cache_dir("get_norm_ds"), "nds.rds")
    if (read_cache && file.exists(cachefile)) return(readRdsVerbose(cachefile))
    ds <- get_ds()
    ids <- names(ds)
    cols <- c("NUM", "RT", "NAME", "DUPOF", "SMILES", "CANONICAL")
    for (i in seq_along(ids)) {
        catf("Processing dataset %s (%d/%d)", ids[i], i, length(ids))
        ds[[i]]$NUM <- seq_len(nrow(ds[[i]]))
        ds[[i]]$CANONICAL <- get_canonical_smiles(ds[[i]]$SMILES) # slow part
        ds[[i]]$DUPOF <- get_dupof(ds[[i]])
        ds[[i]] <- as.data.frame(ds[[i]][, cols])
    }
    if (update_cache) saveRdsVerbose(ds, cachefile)
    invisible(ds)
}

# Plotting #####

#' @noRd
#' @title Heatmap of Metabolites per Column
#' @examples
#' dat = get_norm_ds()
#' plot_met_cc_heatmap(dat)
plot_met_cc_heatmap <- function(dat) {

    # Prapare data frame for plotting
    IDs <- unique(unlist(sapply(dat, function(df) df$ID)))
    mctbl <- sapply(dat, function(df) as.integer(IDs %in% df$ID))
    colnams <- sprintf("%s (%d)", names(dat), sapply(dat, nrow))
    rownames(mctbl) <- IDs
    colnames(mctbl) <- colnams
    nr <- nrow(mctbl)
    nc <- ncol(mctbl)

    # Sort metabolites
    in_Mod <- IDs %in% dat$RP_Steep$ID
    in_Val <- IDs %in% dat$RP_Val_Normal$ID
    in_Normal <- IDs %in% dat$RP_Normal$ID
    in_HILIC <- IDs %in% dat$HILIC$ID
    in_HILIC_Retip <- IDs %in% dat$HILIC_Retip$ID
    idx_sorted <- c(
        which( in_Mod &  in_Normal &  in_HILIC &   in_HILIC_Retip),
        which( in_Mod &  in_Normal &  in_HILIC &  !in_HILIC_Retip),
        which( in_Mod &  in_Normal & !in_HILIC &   in_HILIC_Retip),
        which( in_Mod &  in_Normal & !in_HILIC &  !in_HILIC_Retip),
        which(!in_Mod &  in_Normal &  in_HILIC &  in_HILIC_Retip),
        which(!in_Mod &  in_Normal &  in_HILIC & !in_HILIC_Retip),
        which(!in_Mod &  in_Normal & !in_HILIC &  in_HILIC_Retip),
        which(!in_Mod &  in_Normal & !in_HILIC & !in_HILIC_Retip),
        which( in_Val &               in_HILIC &  in_HILIC_Retip),
        which( in_Val &               in_HILIC & !in_HILIC_Retip),
        which( in_Val &              !in_HILIC &  in_HILIC_Retip),
        which( in_Val &              !in_HILIC & !in_HILIC_Retip),
        which(!in_Val & !in_Normal &  in_HILIC &  in_HILIC_Retip),
        which(!in_Val & !in_Normal &  in_HILIC & !in_HILIC_Retip),
        which(!in_Val & !in_Normal & !in_HILIC &  in_HILIC_Retip)
    )
    stopifnot(sort(idx_sorted) == seq_along(IDs)) # Sanity check
    mctbl <- mctbl[idx_sorted, ]

    # Do Plotting
    opar <- par(mar = c(1.1, 1.1, 6.1, 1.1))
    on.exit(par(opar), add = TRUE, after = FALSE)
    plot(
        x = 0, y = 0,
        xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5),
        xlab = "", ylab = "",
        xaxs = "i", yaxs = "i",
        type = "n", axes = FALSE
    )
    mtext("Metabolites", side = 2, line = 0.1, cex = 0.66)
    for (i in seq_len(nr)) {
        for (j in seq_len(nc)) {
            col <- if (mctbl[i, j] == 1) "#2ECC71" else "#F5F6FA"
            rect1(x = j, y = i, col = col)
        }
    }
    for (j in seq_len(nc)) {
        mtext(colnams[j], side = 3, at = j, las = 2, line = 0.1, cex = 0.5)
    }
}

plot_met_cc_venn <- function(dat, HILIC = TRUE, HILIC_Retip = TRUE) {
    opar <- par(cex = 0.5)
    on.exit(par(opar), add = TRUE, after = FALSE)
    IDs <- unique(unlist(sapply(dat, function(df) df$ID)))
    dat2 <- list(
        RP_Mod = dat$RP_Steep$ID,
        RP_Normal = dat$RP_Normal$ID,
        RP_Val = dat$RP_Val_Normal$ID,
        HILIC = dat$HILIC$ID,
        HILIC_Retip = dat$HILIC_Retip$ID
    )
    if (!HILIC) dat2$HILIC <- NULL
    if (!HILIC_Retip) dat2$HILIC_Retip <- NULL
    names(dat2) <- sprintf("%s (%d)", names(dat2), lengths(dat2))
    venn::venn(dat2, ilabels = "counts", zcolor = "bw")
}

plot_met_cc_upset  <- function(dat) {

}

# Tables #####

table_datasets <- function() {
    ds <- get_ds()
    ids <- names(ds)
    Column <- c("HILIC1", "HILIC2", "RPAXMM1", doc("RP1", 14))
    ColumnTypes <- unique(Column)
    ColumnTypesStr <- paste(ColumnTypes, collapse = ", ")
    Usage <- c(doc("Training", 4), doc("Adjustment", 6), doc("Validation", 7))
    tbl <- data.frame(ids, Column, Usage)
    dir.create("out", showWarnings = FALSE)
    message("")

    for (i in seq_along(ids)) {
        message(sprintf("Processing dataset %s (%d/%d)", ids[i], i, length(ids)))
        ds[[i]]$NUM <- seq_len(nrow(ds[[i]]))
        ds[[i]]$CANONICAL <- get_canonical_smiles(ds[[i]]$SMILES)
        ds[[i]]$DUPOF <- get_dupof(ds[[i]])
        ds[[i]] <- ds[[i]][, c("NUM", "RT", "NAME", "DUPOF", "SMILES", "CANONICAL")]
        path <- file.path("out", paste0(ids[i], ".html"))
        tbl[i, "ID"] <- sprintf("[%s](%s)", ids[i], path)
        tbl[i, "nMeas"] <- nrow(ds[[i]])
        x <- 1
        tbl[i, "nNames"] <- length(unique(ds[[i]]$NAME))
        tbl[i, "nSMILES"] <- length(unique(ds[[i]]$SMILES))
        tbl[i, "nMets"] <- length(unique(ds[[i]]$CANONICAL))
        save_as_html(ds[[i]], path, title = ids[i])
    }
    nColumns <- length(unique(Column))
    nDS <- length(ds)
    x <- 1
    nMeas <- sum(sapply(ds, nrow))
    nMets <- length(unique(unlist(sapply(ds, `[[`, "CANONICAL"))))
}
