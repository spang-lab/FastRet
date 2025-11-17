# Imports (Private) #####
#
#' @import graphics
#' @import grDevices
#' @import parallel
#' @import stats
#' @import utils
#
# IMPORTANT: Only ever import packages shipped with R, as these are loaded
# anyways. Loading 3rd party packages can significantly slow down package
# loading. E.g. glmnet or xgboost both take almost 1s to load.

# Public #####

#' @export
#' @title Get package file
#' @description Returns the path to a file within the FastRet package.
#' @param path The path to the file within the package.
#' @param mustWork If TRUE, an error is thrown if the file does not exist.
#' @return The path to the file.
#' @keywords internal
#' @examples
#' path <- pkg_file("extdata/RP.xlsx")
pkg_file <- function(path, mustWork = FALSE) {
    system.file(path, package = "FastRet", mustWork = mustWork)
}

#' @export
#' @keywords internal
#'
#' @title now
#'
#' @description
#' Returns the current system time formatted according to the provided format
#' string.
#'
#' @param format
#' A string representing the desired time format. Default is "%Y-%m-%d
#' %H:%M:%OS2".
#'
#' @return
#' A string representing the current system time in the specified format.
#'
#' @examples
#' now()            # e.g. "2024-06-12 16:09:32.41"
#' now("%H:%M:%S")  # e.g. "16:09:32"
#'
now <- function(format = "%Y-%m-%d %H:%M:%OS2") {
    format(Sys.time(), format)
}

#' @export
#' @title catf function
#' @description Prints a formatted string with optional prefix and end strings.
#' @param ... Arguments to be passed to sprintf for string formatting.
#' @param prefix A function returning a string to be used as the prefix. Default is a timestamp.
#' @param end A string to be used as the end of the message. Default is a newline character.
#' @return No return value. This function is called for its side effect of printing a message.
#' @keywords internal
#' @examples
#' catf("Hello, %s!", "world")
#' catf("Goodbye", prefix = NULL, end = "!\n")
catf <- function(...,
                 prefix = .Options$FastRet.catf.prefix,
                 end = .Options$FastRet.catf.end) {
    prefixstr <- if (is.null(prefix)) sprintf("%s%s%s ", GREY, now(), RESET) else prefix()
    endstr <- if (is.null(end)) "\n" else end
    middlestr <- sprintf(...)
    msg <- sprintf("%s%s%s", prefixstr, middlestr, endstr)
    cat(msg)
}

#' @export
#' @keywords internal
#'
#' @title Collect elements from a list of lists
#'
#' @description
#' Takes a list of lists where each inner list has the same names. It returns a
#' list where each element corresponds to a name of the inner list that is
#' extracted from each inner list. Especially useful for collecting results from
#' lapply.
#'
#' @param xx
#' A list of lists where each inner list has the same names.
#'
#' @return
#' A list where each element corresponds to a name of the inner list that is
#' extracted from each inner list.
#'
#' @examples
#' xx <- lapply(1:3, function(i) list(a = i, b = i^2, c = i^3))
#' ret <- collect(xx)
#'
collect <- function(xx) {
    ns <- names(xx[[1]])
    ret <- lapply(ns, function(n) sapply(xx, function(x) x[[n]]))
    `names<-`(ret, ns)
}

# Misc (Private) #####

#' @noRd
#' @title Find Random Positions of x in y
#' @description
#' Like `match(x, y)`, but if `y` contains multiple occurrences of
#' `x[i]`, a random occurrence is returned (instead of always the first,
#' like `match` does).
#' @param x Vector of values to match.
#' @param y Vector of values to match against.
#' @param seed Random seed for reproducibility.
match_random <- function(x, y, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    perm <- sample.int(length(y)) # shuffle indices of y
    yshuf <- y[perm]
    pos <- match(x, yshuf) # first match in shuffled y -> random occurrence
    perm[pos] # map back to original indices
}

#' @noRd
#' @title Create automatically named List
#' @description
#' Like normal `list()`, except that unnamed elements are automatically named according to their symbol
#'
#' COPIED OVER FROM TOSCUTIL. Can be replaced with original toscutil version as
#' soon as all NAMESPACE imports from toscutil have been removed (right now
#' loading toscutil takes 22ms, and we want to avoid that).
#'
#' @param ... List elements
#' @return Object of type `list` with names attribute set
#' @seealso [list()]
#' @keywords base
#' @examples
#' a <- 1:10
#' b <- "helloworld"
#' l1 <- list(a, b)
#' names(l1) <- c("a", "b")
#' l2 <- named(a, b)
#' stopifnot(identical(l1, l2))
#' l3 <- list(z = a, b = b)
#' l4 <- named(z = a, b)
#' stopifnot(identical(l3, l4))
named <- function(...) {
    .symbols <- as.character(substitute(list(...)))[-1]
    .elems <- list(...)
    .idx <- if (is.null(names(.elems))) {
        rep(TRUE, length(.elems))
    } else {
        names(.elems) == ""
    }
    names(.elems)[.idx] <- .symbols[.idx]
    .elems
}

#' @noRd
#' @title Null Function
#' @description A function that always returns invisibly NULL, ignoring all arguments.
null <- function(...) {
    invisible(NULL)
}

#' @noRd
#' @title Not In Operator
#' @description Inverse of the %in% operator.
#' @param x Vector of values to test
#' @param y Vector of values to test against
#' @return Logical vector indicating which elements of x are not in y
`%notin%` <- function(x, y) {
    !(x %in% y)
}

`%||%` <- function(a, b) {
    if (!is.null(a)) a else b
}

#' @noRd
#' @title Convert Vector to String
#' @description Converts a vector to a string representation.
#' @param x Vector
#' @return String representation of the vector
#' @examples
#' as_str(x = c(1, 2, 3))                   # "1, 2, 3"
#' as_str(x = c(a = 1, b = 2, c = "Hello")) # "a=1, b=2, c=Hello"
#' as_str(x = c(a = 1, 2, c = "Hello"))     # "a=1, 2, c=Hello"
as_str <- function(x) {
    vals <- as.character(x)
    nams <- names(x) %||% rep("", length(x))
    sep <- ifelse(nams == "", "", "=")
    paste(nams, sep, vals, sep = "", collapse = ", ")
}

#' @noRd
#' @title Canonicalize SMILES
#' @description Convert SMILES to canonical form.
#' @param smiles Character vector of SMILES.
#' @return Character vector of canonical SMILES.
#' @examples
#' as_canonical(c("CCO", "C(C)O"))
as_canonical <- function(smiles) {
    canons <- rep("INVALID", length(smiles))
    is_valid <- is_valid_smiles(smiles)
    molecules <- rcdk::parse.smiles(smiles[is_valid])
    flavor <- rcdk::smiles.flavors("Canonical")
    canons[is_valid] <- sapply(molecules, rcdk::get.smiles, flavor)
    canons
}

#' @noRd
#' @title Validate SMILES
#' @description Check if SMILES strings are valid.
#' @details Attention: this documentation has been generated automatically by
#' GitHub Copilot and has NOT been reviewed.
#' @param smiles Character vector of SMILES.
#' @return Logical vector of validity.
#' @examples
#' is_valid_smiles(c("CCO", "invalid"))
is_valid_smiles <- function(smiles) {
    molecules <- suppressWarnings(rcdk::parse.smiles(smiles))
    !sapply(molecules, is.null)
}

#' @noRd
#' @description Calculate the number of workers for parallel processing
#' @param mult A multiplier for the number of cores. Default is 0.5.
#' @param nmax The maximum number of workers allowed. Default is 16.
#' @return The number of workers to be used for parallel processing.
#' @examples
#' get_n_workers(0.5)      # returns  2 on a  4-core machine
#' get_n_workers(1.0)      # returns  4 on a  4-core machine
#' get_n_workers(2.0)      # returns  8 on a  4-core machine
#' get_n_workers(1.0,  4)  # returns  4 on a 32-core machine
#' get_n_workers(1.0, 16)  # returns 16 on a 32-core machine
#' get_n_workers(1.0, 64)  # returns 32 on a 32-core machine
get_n_workers <- function(mult = 0.5, nmax = 16) {
    n <- parallel::detectCores()
    nmul <- ceiling(n * mult)
    min(nmul, nmax)
}

# Development Helpers (Private) #####

load_all <- function() {
    x <- Sys.time()
    catf("Calling: pkgload::load_all(reset = TRUE)")
    pkgload::load_all(reset = TRUE, quiet = TRUE)
    catf("Calling: pkgload_env$insert_global_shims(force = TRUE)")
    pkgload_env <- environment(pkgload::load_all)
    pkgload_env$insert_global_shims(force = TRUE)
    diff <- Sys.time() - x
    catf("Elapsed: %s", format(diff))
}

document <- function() {
    x <- Sys.time()
    catf("Calling: devtools::document(quiet = TRUE)")
    devtools::document(quiet = TRUE)
    catf("Calling: pkgload_env$insert_global_shims(force = TRUE)")
    pkgload_env <- environment(pkgload::load_all)
    pkgload_env$insert_global_shims(force = TRUE)
    diff <- Sys.time() - x
    catf("Elapsed: %s", format(diff))
}

serve_docs <- function() {
    servr::httd("docs")
}

make_docs <- function(reload = TRUE) {
    if (reload) {
        pkgdown::build_site(lazy = TRUE)
        vignette_files <- dir(pkg_file("vignettes"), recursive = TRUE, full.names = TRUE)
        timestamps_old <- file.mtime(vignette_files)
        while (TRUE) {
            Sys.sleep(0.1)
            timestamps_new <- file.mtime(vignette_files)
            if (!identical(timestamps_old, timestamps_new)) {
                pkgdown::build_site(lazy = TRUE)
                timestamps_old <- timestamps_new
            }
        }
    } else {
        pkgdown::build_site(lazy = TRUE)
    }
}

#' @noRd
#' @description
#' To enable either call add "eval(.Options$FastRet.onFuncEntry)" at the
#' beginning of each function or call `trace(f,
#' quote(eval(.Options$FastRet.onFuncEntry)))` for each function `f` you want to
#' trace.
enable_function_tracing <- function() {
    onFuncEntry <- quote({
        funcname <- as.character(sys.call(-2))
        cat(sprintf("Start: %s\n", funcname))
        on.exit(cat(sprintf("Exit: %s\n", funcname)), add = TRUE)
    })
    options(FastRet.onFuncEntry = onFuncEntry)
}

# Caret Replacements (Private) #####

#' @noRd
#' @title Create cross-validation folds
#' @param y Vector of indices/numeric values for splitting.
#' @param k
#' Number of folds. Default 10. If `k` is larger equal `length(y)`, k will be
#' silently set to `length(y)`.
#' @return A list of integer vectors with fold indices.
createFolds <- function(y, k = 10) {
    if (k >= length(y)) k <- length(y)
    foldnrs <- sample(rep(1:k, length.out = length(y)))
    folds <- split(seq_along(y), foldnrs)
    names(folds) <- paste0("Fold", seq_len(k))
    folds
}

#' @noRd
#' @title nearZeroVar
#' @description
#' Lightweight replacement for `caret::nearZeroVar`. Identifies predictors with
#' (a) a very large frequency ratio of the most common value to the second most
#' common value and (b) a low percent of distinct values. (a) and (b) must both
#' be true at the same time to consider a variable as near-zero-variance. This
#' mirrors the behavior used from `caret::nearZeroVar`.
#' @param X A data.frame.
#' @param freqCut Frequency ratio cutoff (default 95/5).
#' @param uniqueCut Percent unique cutoff (default 10).
#' @return
#' An integer vector with column indices of "near-zero-variance" predictors.
nearZeroVar <- function(X, freqCut = 95/5, uniqueCut = 10) {
  X <- as.data.frame(X, stringsAsFactors = FALSE)
  which(vapply(X, hasNearZeroVar, logical(1), freqCut, uniqueCut, USE.NAMES = FALSE))
}

#' @noRd
#' @title Check whether a predictor has near-zero variance
#' @description
#' Checks whether a predictor has near-zero variance based on frequency ratio
#' and percent unique values. Mirrors the behavior of `caret::nearZeroVar`.
#' @param x A vector of predictor values.
#' @param freqCut Frequency ratio cutoff (default 95/5).
#' @param uniqueCut Percent unique cutoff (default 10).
#' @return
#' A logical value indicating whether the predictor has near-zero variance.
hasNearZeroVar <- function(x, freqCut = 95/5, uniqueCut = 10) {
    n <- length(x)
    x <- x[!is.na(x)]
    if (!length(x)) return(TRUE)
    tab <- table(x)
    nUniq <- length(tab)
    if (nUniq <= 1) return(TRUE)
    top <- sort(tab, TRUE)
    freqRatio <- if (length(top) == 1) Inf else top[1] / top[2]
    pctUnique <- 100 * nUniq / n
    freqRatio > freqCut & pctUnique < uniqueCut || (is.numeric(x) && sd(x) == 0)
}

# Colors (Private) #####

GREY <- "\033[1;30m"
RED <- "\033[1;31m"
GREEN <- "\033[1;32m"
YELLOW <- "\033[1;33m"
BLUE <- "\033[1;34m"
RESET <- "\033[0m"
