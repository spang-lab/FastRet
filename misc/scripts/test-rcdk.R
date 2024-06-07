library(rcdk, quietly = TRUE)

ok <- "\033[32mok\033[39m\n" # green
failed <- "\033[31mfailed\033[39m\n" # red

cat('Testing rcdk ... ')
data(bpdata)
mols <- parse.smiles(bpdata$SMILES)
nams  <- get.desc.names("all")
X <- suppressWarnings(eval.desc(mols[1:10], nams))
exitcode <- if (nrow(X) == 10 && !all(is.na(X))) 0 else 1
cat(if (exitcode == 0) ok else failed)
quit(status = exitcode, save = 'no')
