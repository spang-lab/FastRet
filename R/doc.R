# Document Environment #####

de <- local({
    rbuildignore <- system.file(".Rbuildignore", package = "FastRet")
    loaded_via_devtools <- file.exists(rbuildignore)
    if (loaded_via_devtools && !is.null(.GlobalEnv$de)) {
        de <- .GlobalEnv$de
    } else {
        de <- list(docs = list(), cur = 0, doc = NULL)
        de <- structure(as.environment(de), class = "de")
        if (loaded_via_devtools) .GlobalEnv$de <- de
    }
    de
})

doc.cur <- function() {
    de$cur
}

doc.new <- function(
    # Metadata
    title = NULL, # Document Title as string
    subtitle = NULL, # Document Subtitle as string
    author = NULL, # Document Authors as string
    # Debugging
    print_tbls = FALSE, # Print tables to stdout after creation?
    print_secs = TRUE, # Print sections to stdout after creation?
    show_figs = TRUE, # Show figures in grafical window after creation?
    # Pandoc
    work_dir = tempfile("dir"), # Path where files should be written to
    out_dir = ".", # Path where the output file should be written to
    out_file = "out", # Name of the output file without extension
    formats = "html", # Pandoc convert formats (html, pdf, docx and tex)
    defaults_yml = NULL, # Path to the `defaults.yml` file used by pandoc
    preamble_tex = NULL, # Path to the `preamble.tex` file used by pandoc
    reference_docx = NULL, # Path to the Word template file used by pandoc
    template_html = NULL, # Path to the HTML template file used by pandoc
    template_tex = NULL, # Path to the LaTeX template file used by pandoc
    open = FALSE, # Automatically open pandoc result after generation?
    # Counters
    content = c(""), # Document body, used by most other functions
    # Content
    blk_ctr = 2, # Block counter, used by most functions (1)
    sec_ctr = 0, # Section counter, used by `h()`
    fig_ctr = 0, # Figure counter, used by `fig()`
    tbl_ctr = 0, # Table counter, used by `tbl()`
    lst_ctr = 0, # List counter, used e.g. by `lst()`
    lst_sym = list(), # Symbol of every started list, e.g. "-" or "1."
    lst_ws = list(), # Indentation of every started list, e.g. "  "
    lst_ws0 = "", # Sum of lists indents with space for current `lst_sym`
    lst_ws1 = "" # Sum of list indents without space for current `lst_sym`
    #
    # (1) We consider every piece of text a "block", e.g. a paragraph, a
    # heading, a list item, etc. We start at block 2, because block 1 is
    # reserved for the yaml header.
    #
) {
    #
    # FUNCTION BODY
    #
    de$doc <- structure(as.list(environment()), class = "doc")
    de$docs[[de$cur + 1]] <- de$doc
    de$cur <- de$cur + 1
    invisible(de$doc)
}

doc.off <- function() {
    if (de$cur == 0) stop("Cannot shut down document 0.", call. = FALSE)
    de$docs[[de$cur]] <- NULL
    de$cur <- de$cur - 1
    de$doc <- if (de$cur == 0) NULL else de$docs[[de$cur]]
    de$cur
}

doc.write <- function(x, collapse = "\n", end = "\n") {
    if (de$cur == 0) stop("No open document. Call `doc.new()` first.", call. = FALSE)
    x <- paste(x, collapse = "\n")
    x <- paste0(x, end)
    blk <- de$doc$blk_ctr
    de$doc$content[[de$doc$blk_ctr]] <- x
    de$doc$blk_ctr <- de$doc$blk_ctr + 1
}

doc.get <- function() {
    if (de$cur == 0) stop("No open document. Call `doc.new()` first.", call. = FALSE)
    de$docs[[de$cur]]
}

doc.show <- function() {
    tty <- interactive() && isatty(stdout())
    if (tty) cat("\033[34m")
    for (blk in de$doc$content) cat(blk)
    if (tty) cat("\033[0m")
}

doc.save <- function() {
    if (de$cur == 0) stop("No open document. Call `doc.new()` first.", call. = FALSE)
    doc <- de$doc
    dir.create(de$doc$work_dir, recursive = TRUE, showWarnings = FALSE)
    md_path <- file.path(de$doc$work_dir, paste0(de$doc$out_file, ".md"))
    writeLines(de$doc$content, md_path, sep = "")
    catf("Document saved to %s", md_path)
    if (de$doc$open) utils::browseURL(md_path)
    invisible(md_path)
}

doc.fetch <- function(overwrite = TRUE) {
    if (de$cur == 0) stop("No open document. Call `doc.new()` first.", call. = FALSE)
    src <- file.path(de$doc$work_dir, paste0(de$doc$out_file, ".md"))
    dst <- file.path(de$doc$out_dir, paste0(de$doc$out_file, ".md"))
    dir.create(de$doc$out_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(from=src, to=dst, overwrite = overwrite)
    catf("Document copied to %s", dst)
}

# S3 Methods #####

#' @export
print.de <- function(x, ...) print(ls.str(x))

#' @export
print.doc <- function(x, ...) {
    print(ls.str(x))
}

# Headers #####

#' @title Print YAML Header
#' @noRd
#' @keywords md
#' @description
#' Markdown documents converted via pandoc can contain a YAML header containing
#' meta information about the document, such as title, author and generation
#' date. This function creates such a YAML header.
#' @param title The title of the document.
#' @param author The author of the document.
#' @param date The date of the document. Defaults to the current date in "8. August, 2025" format.
#' @param subtitle An optional subtitle for the document.
#' @param ... Reserved for future use, currently ignored.
#' @return NULL. Called for side effect of printing to stdout.
#' @examples
#' meta() # Uses title and subtitle from `opts`
#' meta(op = list(title = "My Report", subtitle = "A short description"))
meta <- function(title,
                 author,
                 date = format(Sys.time(), "%d. %B, %Y"),
                 subtitle = NULL,
                 ...) {
    if (de$cur == 0) stop("No open document. Call `doc.new()` first.", call. = FALSE)
    lines <- c(
        paste("---"),
        paste("title:", title),
        paste("author:", author),
        paste("date:", format(Sys.time(), "%d. %B, %Y")),
        if (!is.null(subtitle)) paste("subtitle:", subtitle),
        paste("---\n")
    )
    doc.write(lines)
}

#' @noRd
#' @keywords md
#' @title Print Markdown Section Heading
#' @description
#' Print a markdown heading with the specified level.
#' @param x Heading text.
#' @param n Level of the heading. Can be ommitted if x starts with hashes.
#' @examples \dontrun{
#' sec("# Hello World") # prints '# Hello World'
#' sec("## Hello World") # prints '## Hello World'
#' sec("Hello World", 2) # prints '## Hello World'
#' }
sec <- function(x, n = NULL) {
    if (de$cur == 0) stop("No open document. Call `doc.new()` first.", call. = FALSE)
    if (is.null(n)) {
        n <- attr(regexpr("^#+", x), "match.length")
        if (n == -1) stop("Please provide a valid heading level.", call. = FALSE)
    } else {
        hashs <- paste(rep("#", n), collapse = "")
        x <- paste(hashs, x)
    }
    de$doc$sec_ctr <- n
    if (!endsWith("\n", x)) x <- paste0(x, "\n")
    if (de$doc$print_secs) catf(trimws(x))
    doc.write(x)
}


# Others #####

#' @noRd
#' @keywords md
#' @title Print a Dataframe as Markdown Table
#' @description Prints a dataframe as table using [pander::pandoc.table()].
#' @param df The data frame to be converted into a table.
#' @param caption The caption for the table.
#' @param justify The justification of the table ("left", "right", "center").
#' @param row.names A logical value indicating whether to include row names.
#' @param col.names The column names for the table. If NULL, the column names of df are used.
#' @param widths The maximum width of the cells in the table.
#' @param split.tables The maximum number of rows in a table before it is split into multiple tables.
#' @param use.hyphening A logical value indicating whether to use hyphenation in the table.
#' @param keep.line.breaks A logical value indicating whether to keep line breaks in the table.
#' @param ... Additional arguments passed to pandoc.table.
#' @return A table created using pandoc.table.
#' @examples
#' tbl(data.frame(a = 1:5, b = 6:10), caption = "Example table")
tbl <- function(df,
                id = NULL,
                caption = "",
                justify = "left",
                row.names = FALSE,
                col.names = NULL,
                widths = pander::panderOptions("table.split.cells"),
                split.tables = Inf,
                use.hyphening = FALSE,
                keep.line.breaks = TRUE,
                ...) {
    colnames(df) <- make_breakable(colnames(df))
    if (is.character(id)) doc.write(sprintf("<a id = %s></a>", id))
    doc.write(pander::pandoc.table.return(
        df,
        caption = caption,
        row.names = row.names,
        col.names = col.names %||% colnames(df),
        justify = justify,
        split.cells = widths,
        split.tables = split.tables,
        use.hyphening = use.hyphening,
        keep.line.breaks = keep.line.breaks,
        ...
    ))
}

#' @noRd
#' @keywords md
#' @title Save figure and print its caption
#' @description
#' Saves the plot produced by expr as SVG and prints a Markdown image tag.
#' @param expr An expression that produces a plot.
#' @param name The name of the output file (without extension).
#' @param caption The caption for the image in the Markdown tag.
#' @param width The width of the output images in inches. Default is 7.
#' @param height The height of the output images in inches. Default is 5.
#' @title Save figure and print its caption
#' @return A Markdown image tag for the SVG file.
#' @examples \dontrun{
#' # Create a plot object
#' map <- ggplot2::aes(x = mpg, y = disp)
#' expr <- ggplot2::plot(mtcars, map)
#' expr <- expr + ggplot2::geom_point()
#' # Save the plot and generate a Markdown tag
#' fig("my_plot", expr, "This is a scatter plot of mtcars data.")
#' }
fig <- function(expr, name, caption = "", width = 7, height = 5) {
    svg_abs <- file.path(de$doc$work_dir, "figures", paste0(name, ".svg"))
    svg_rel <- paste0("figures/", name, ".svg")
    dir.create(dirname(svg_abs), recursive = TRUE, showWarnings = FALSE)
    svg(filename = svg_abs, width = width, height = height)
    tryCatch(expr, finally = dev.off())
    doc.write(sprintf("![%s](%s)\n", caption, svg_rel))
}

#' @noRd
#' @keywords md
#' @title Print a paragraph
#' @description Print a markdown paragraph.
#' @param x Text of the paragraph.
#' @examples
#' para("Dies ist ein deutscher Absatz.", "This is an English paragraph.")
para <- function(x) {
    doc.write(x)
    doc.write("")
}

# Listings #####

#' @noRd
#' @keywords md
#' @title Tags for writing Markdown Lists
#' @description
#' [lst()] starts a new list (can be either unordered or ordered),
#' [ulst()] starts a new unordered list,
#' [olsl()] starts a new ordered list,
#' [elst()] ends the current list and
#' [item()] adds an item to the list.
#' See the 'Examples' section for further details.
#' @examples
#' local({
#'      doc.new()
#'      on.exit(doc.off())
#'      para("Paragaph before")
#'      ulst()
#'      item("Unordered 1")
#'      item("Unordered 2")
#'      olst()
#'      item("Ordered 2.1")
#'      item("Ordered 2.2")
#'      ulst()
#'      item("Unordered 2.2.1")
#'      elst()
#'      elst()
#'      elst()
#'      para("Paragaph after")
#'      doc.show()
#' })
lst <- function(typ = "ul") {
    x <- de$doc$lst_ctr + 1
    de$doc$lst_ctr <- x
    de$doc$lst_sym[x] <- if (typ == "ul") "- " else "1. "
    de$doc$lst_ws[x] <- if (typ == "ul") "  " else "   "
    de$doc$lst_ws0 <- paste(de$doc$lst_ws[-x], collapse = "")
    de$doc$lst_ws1 <- paste(de$doc$lst_ws, collapse = "")
}

#' @noRd
#' @rdName lst
ulst <- function() lst("ul")

#' @noRd
#' @rdName lst
olst <- function() lst("ol")

#' @noRd
#' @rdname lst
elst <- function() {
    x <- de$doc$lst_ctr
    if (x == 0) stop("Not in list. Use lst() or ol() first.", call. = FALSE)
    de$doc$lst_sym[[x]] <- NULL
    de$doc$lst_ws[x] <- NULL
    x <- x - 1
    de$doc$lst_ctr <- x
    de$doc$lst_ws0 <- paste(de$doc$lst_ws[-x], collapse = "")
    de$doc$lst_ws1 <- paste(de$doc$lst_ws, collapse = "")
    if (x == 0) doc.write("") # End nested list with empty line
}

#' @noRd
#' @param x Text of the list item.
#' @rdname lst
item <- function(x) {
    if (de$doc$lst_ctr == 0) stop("Not in list. Use lst() or ol() first.", call. = FALSE)
    sym <- de$doc$lst_sym[[de$doc$lst_ctr]]
    lines <- unlist(strsplit(x, "\n"))
    ws <- rep(de$doc$lst_ws1, length(lines))
    ws[1] <- paste0(de$doc$lst_ws0, sym)
    indented_lines <- paste(ws, lines, sep = "")
    doc.write(indented_lines)
}

# Testing #####

test_doc_gen <- function() {
    while (de$cur >= 1) doc.off()
    doc.new()
    meta(title = "Test Document",
         author = "John Doe",
         date = Sys.Date(),
         subtitle = "This is a test document")

    sec("# Introduction")

    para("This is a test paragraph in the introduction section.")
    ulst()
    item("First item in the list")
    item("Second item in the list")
    olst()
    item("First item in the ordered list")
    item("Second item in the ordered list")
    elst()
    elst()

    sec("# Results")

    sec("## The Iris Dataset")
    tbl(df = iris,
        name = "table-iris",
        caption = "Iris Dataset")

    sec("## Scatter Plot of Sepal Length vs. Sepal Width")
    fig(expr = plot(iris$Sepal.Length, iris$Sepal.Width),
        name = "figure-iris",
        caption = "Scatter plot of Sepal Length vs. Sepal Width"
    )

    sec("# Conclusion")

    para("This is a test paragraph in the conclusion section.")
    doc.write("This is a test document generated by the `test_doc_gen()` function.")


    # doc.show()
    doc.save()
    doc.fetch()
    doc.off()
}
