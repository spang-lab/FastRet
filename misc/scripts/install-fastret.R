#!/usr/bin/Rscript
options(Ncpus = parallel::detectCores())
remotes::install_github("spang-lab/fastret", dependencies = TRUE, upgrade = 'always')
