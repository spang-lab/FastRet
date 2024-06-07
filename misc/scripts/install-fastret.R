#!/usr/bin/Rscript
options(Ncpus = parallel::detectCores())
devtools::install_github("spang-lab/fastret", dependencies = TRUE, upgrade = 'always')
