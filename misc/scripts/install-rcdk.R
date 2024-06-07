# Usage: Rscript install-rcdk.R [--verbose]

args <- commandArgs(trailingOnly = TRUE)
verbose <- "--verbose" %in% args

options(Ncpus = parallel::detectCores())
pkgs <- c('rcdk', 'rcdklibs', 'rJava')
installed <- pkgs %in% rownames(installed.packages())
exitcode <- 0
ok <- "\033[32mok\033[39m\n" # green
failed <- "\033[31mfailed\033[39m\n" # red

with_sink <- function(expr) {
    if (!verbose) {
        con <- file(nullfile(), open = "wt")
        on.exit(close(con), add = TRUE, after = FALSE)
        sink(con)
        on.exit(sink(), add = TRUE, after = FALSE)
        sink(con, type = "message")
        on.exit(sink(type = "message"), add = TRUE, after = FALSE)
    }
    expr
}

cat('Removing rJava, rcdklibs and rcdk ... ')
for (lib in .libPaths()) {
    with_sink(try(remove.packages(pkgs[installed], lib = .libPaths())))
}
importable <- with_sink(suppressWarnings(sapply(pkgs, library, character.only = TRUE, logical.return = TRUE, quietly = TRUE)))
if (any(importable)) {
    cat(failed)
    exitcode <- 1
} else {
    cat(ok)
}

cat('Installing rJava, rcdklibs and rcdk ... ')
with_sink({
    install.packages('rJava',    type = 'source', repos = 'https://cloud.r-project.org/', quiet = TRUE)
    install.packages('rcdklibs', type = 'source', repos = 'https://cloud.r-project.org/', quiet = TRUE)
    install.packages('rcdk',     type = 'source', repos = 'https://cloud.r-project.org/', quiet = TRUE)
})
importable <- sapply(pkgs, library, character.only = TRUE, logical.return = TRUE, quietly = TRUE)
if (all(importable)) {
    cat(ok)
} else {
    cat(failed)
    exitcode <- 1
}

quit(status = exitcode, save = 'no')
