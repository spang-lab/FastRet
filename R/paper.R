# Main #####

reproduce_paper <- function() {
    # [x] Implement function [make_metabolites_df()].
    # [x] Print all the tables in the a markdown document.
    # [ ] Convert the markdown document to PDF using the typst engine. <-- CONTINUE HERE
    # [ ] Add the intro text from the old Rmd.
    # [ ] Add the figures from the old Rmd.
    # [ ] Add some checks to the dataset objects during [make_datasets_list()] to ensure that all links are valid.

    Columns <- make_columns_df()
    Conditions <- make_conditions_df()
    Gradients <- make_gradients_df()
    Eluents <- make_eluents_df()
    Datasets_List <- make_datasets_list()
    Datasets <- make_datasets_df(Datasets_List)
    Measurements <- make_measurements_df(Datasets_List)
    Metabolites <- make_metabolites_df(Datasets_List)

    intro <- glue::glue(
        "To evaluate FastRet, we used metabolite retention time measurements from 4 different chromatography columns: RP, RP_AXMM, HILIC and HILIC_Retip. The RP column was used under seven different chromatographic conditions and with two distinct metabolite sets, resulting in a total of 17 datasets comprising 2546 measurements of 1256 metabolites. Table [Datasets](#table-datasets) summarizes the key properties of each dataset."
    )
    stopifnot(
        all(nrow(Columns) == 4),
        all(sort(Columns$ID) == sort(c("RP", "RP_AXMM", "HILIC", "HILIC_Retip"))),
        length(unique(Datasets$ConditionID[Datasets$ColumnID == "RP"])) == 7,
        nrow(Datasets) == 17,
        nrow(Measurements) == 2546,
        nrow(Metabolites) == 1256
    )

    while (doc.cur() != 0) doc.off()
    doc <- doc.new()
    sec("## Table: Columns"); tbl(Columns)
    sec("## Table: Conditions"); tbl(Conditions)
    sec("## Table: Gradients"); tbl(Gradients)
    sec("## Table: Eluents"); tbl(Eluents)
    sec("## Table: Datasets"); tbl(Datasets)
    # sec("## Table: Metabolites"); tbl(Metabolites)
    sec("## Datasets");
    for (id in names(Datasets_List)[2:4]) {
        sec(sprintf("### Table: %s", id))
        x <- Datasets_List[[id]]
        x$SMILES <- make_breakable(x$SMILES, n = 20)
        x$CANONICAL <- make_breakable(x$CANONICAL, n = 20)
        tbl(x)
    }
    doc.save()
    doc.fetch()
    doc.off()

    par(sprintf(
        "The full list of metabolites as canonical SMILES strings is provided in section [Metabolites](#metabolites). For each metabolite, the datasets in which it was measured, the number of measurements, the input name and the input SMILES string entered by the experimenter are listed. This information is also visualized in several ways: in Figure [Met-DS-Heatmap](#met-dsl-heatmap), a heatmap of measurement counts per metabolite and dataset is presented; in Figure [Met-DS-Venn](#met-dsl-venn), Venn diagrams illustrating metabolite overlap and uniqueness across datasets are shown; and in Figure [Met-DS-Upset](#met-dsl-upset), an UpSet plot of the same information is displayed."
    ))

    tbl(df = tblDS, caption = paste(sep = "\n",
        "**Table Datasets:**",
        "Overview of datasets used for retention time prediction.",
        "*ID:* Dataset identifier.",
        "*Column:* Chromatographic Column Identifier.",
        "*Usage:* Dataset purpose.",
        "*nMeas:* Number of unique measurements.",
        "*nNames:* Number of unique metasbolite input names.",
        "*nSMILES:* Number of unique input SMILES.",
        "*nMets:* Number of unique canonical SMILES.",
        "'Input' refers to names and SMILES as entered by the experimenter.",
        "'Canonical' refers to standardized SMILES derived from input SMILES."
    ))

    # VALIDATE: number of metabolites in in-house RT library (~380)
    pdf.open("out/heatmap_met_cc.pdf")
    try(plot_met_cc_heatmap(ds))
    dev.off()

    pdf.open("out/venn_met_cc.pdf")
    try({
        plot_met_cc_venn(ds, HILIC = FALSE, HILIC_Retip = FALSE)
        plot_met_cc_venn(ds, HILIC = TRUE, HILIC_Retip = FALSE)
        plot_met_cc_venn(ds, HILIC = TRUE, HILIC_Retip = TRUE)
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

# Analysis #####

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

# Helpers #####

make_breakable <- function(x,
                           n = NULL,
                           at_underscore = TRUE,
                           sep = c(" ", "_\u200b")[1]) {
    patt <- paste0("(.{", n, "})")
    repl <- paste0("\\1", sep)
    if (is.numeric(n)) x <- gsub(patt, repl, x)
    if (isTRUE(at_underscore)) x <- gsub("_", repl, x)
    x
}

#' @noRd
#' @title Pretty Print Data Frame
pp <- function(x, ...,
               digits = 3, # base default is NULL
               quote = FALSE, # base default is FALSE
               right = FALSE, # base default is TRUE
               row.names = NULL, # base default is TRUE
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
            xhead[[k]] <- rep("SUBLIST", nrows)
        }
    }
    if (is.null(row.names)) {
        rnams <- rownames(xhead)
        rseq <- as.character(seq_len(nrows))
        row.names <- if (identical(rnams, rseq)) FALSE else TRUE
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


# Plotting #####

#' @noRd
#' @title Heatmap of Metabolites per Column
#' @examples
#' ds = get_norm_ds()
#' plot_met_cc_heatmap(ds)
plot_met_cc_heatmap <- function(ds) {

    # Prapare data frame for plotting
    IDs <- unique(unlist(sapply(ds, function(df) df$ID)))
    mctbl <- sapply(ds, function(df) as.integer(IDs %in% df$ID))
    colnams <- sprintf("%s (%d)", names(ds), sapply(ds, nrow))
    rownames(mctbl) <- IDs
    colnames(mctbl) <- colnams
    nr <- nrow(mctbl)
    nc <- ncol(mctbl)

    # Sort metabolites
    in_Mod <- IDs %in% ds$RP_Steep$ID
    in_Val <- IDs %in% ds$RP_Val_Normal$ID
    in_Normal <- IDs %in% ds$RP_Normal$ID
    in_HILIC <- IDs %in% ds$HILIC$ID
    in_HILIC_Retip <- IDs %in% ds$HILIC_Retip$ID
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

plot_met_cc_venn <- function(ds, HILIC = TRUE, HILIC_Retip = TRUE) {
    opar <- par(cex = 0.5)
    on.exit(par(opar), add = TRUE, after = FALSE)
    IDs <- unique(unlist(sapply(ds, function(df) df$ID)))
    dat2 <- list(
        RP_Mod = ds$RP_Steep$ID,
        RP_Normal = ds$RP_Normal$ID,
        RP_Val = ds$RP_Val_Normal$ID,
        HILIC = ds$HILIC$ID,
        HILIC_Retip = ds$HILIC_Retip$ID
    )
    if (!HILIC) dat2$HILIC <- NULL
    if (!HILIC_Retip) dat2$HILIC_Retip <- NULL
    names(dat2) <- sprintf("%s (%d)", names(dat2), lengths(dat2))
    venn::venn(dat2, ilabels = "counts", zcolor = "bw")
}

plot_met_cc_upset  <- function(ds) {

}

rect1 <- function(x, y, col) {
    rect(
        xleft = x - 0.5, ybottom = y - 0.5,
        xright = x + 0.5, ytop = y + 0.5,
        col = col, border = NULL, density = NA
    )
}

pdf.open <- function(path, width = 3.33, height = 3.33) {
    catf("Saving %s", path)
    if (!dir.exists("out")) dir.create("out", recursive = TRUE)
    pdf(path, width = width, height = height)
}

#' @noRd
#' @title Save a data.frame as an interactive standalone HTML file
#' @description
#' This function creates an interactive, sortable, filterable HTML table from
#' a data.frame and saves it as a standalone HTML file using the ds and
#' htmlwidgets packages.
#' @param x A data.frame or tibble to be rendered.
#' @param file Character. Output HTML file path (e.g. "table.html").
#' @param title Character. Optional HTML page title.
#' @param caption Character. Optional table caption.
#' @param pageLength Integer. Number of rows to show per page (default is 25).
#' @param row.names Logical. Whether to include row names. Default is FALSE.
#' @return Invisibly returns the ds widget object.
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


