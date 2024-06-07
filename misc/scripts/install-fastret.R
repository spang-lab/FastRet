#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
branch_idx <- match("--branch", args)
branch <- if (is.na(branch_idx)) "main" else args[branch_idx + 1]
verbose_idx <- match("--verbose", args)
verbose <- !is.na(verbose_idx)
if (!verbose) {
    sink("/dev/null")
    on.exit(sink())
}
options(Ncpus = parallel::detectCores())
remotes::install_github("spang-lab/fastret", dependencies = TRUE, upgrade = 'always', ref = branch)
