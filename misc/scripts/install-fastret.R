#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
help_idx <- match("--help", args)
show_help <- !is.na(help_idx)
if (show_help) {
    cat("Usage: install-fastret.R [--branch BRANCH] [--verbose] [--local]\n")
    cat("Options:\n")
    cat("  --branch BRANCH  Install a specific branch of fastret (default: main)\n")
    cat("  --verbose        Show verbose output\n")
    cat("  --local          Install fastret from local source\n")
    q(save = "no", status = 0)
}
branch_idx <- match("--branch", args)
branch <- if (is.na(branch_idx)) "main" else args[branch_idx + 1]
verbose_idx <- match("--verbose", args)
verbose <- !is.na(verbose_idx)
local_idx <- match("--local", args)
local <- !is.na(local_idx)
if (!verbose) {
    sink("/dev/null")
    on.exit(sink())
}
options(Ncpus = parallel::detectCores())
if (local) {
    devtools::install(".", dependencies = TRUE, upgrade = 'always')
} else {
    remotes::install_github("spang-lab/fastret", dependencies = TRUE, upgrade = 'always', ref = branch)
}
