#!/usr/bin/Rscript
options(Ncpus = parallel::detectCores())
install.packages(c('devtools', 'languageserver'), repos = 'https://cloud.r-project.org/')
