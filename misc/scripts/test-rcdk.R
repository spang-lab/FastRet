# Parse Args
args <- commandArgs(trailingOnly = TRUE)
help <- "--help" %in% args
if (help) {
    cat("Usage: Rscript test-rcdk.R [--verbose] [--progress] [--descriptors x:y] [--smiles a:b]\n")
    cat("Example: Rscript test-rcdk.R --verbose --descriptors 1:19 --smiles 190:210\n")
    cat("Ranges: x: 1..277, y: x..277, a: a..53, b: a..277\n")
    cat("\n")
    quit(status = 0)
}
verbose <- "--verbose" %in% args
progress <- "--progress" %in% args
options(warn = if (verbose) 1 else -1)
desc_idx <- match("--descriptors", args)
desc_str <- if (is.na(desc_idx)) "1:53" else args[desc_idx + 1]
desc_ids <- eval(parse(text = desc_str))
smi_idx <- match("--smiles", args)
smi_str <- if (is.na(smi_idx)) "1:277" else args[smi_idx + 1]
smi_ids <- eval(parse(text = smi_str))
if (!all(smi_ids) %in% 1:277) stop("Invalid SMILES range")
if (!all(desc_ids) %in% 1:53) stop("Invalid descriptors range")

# Load and configure rCDK
options("java.parameters" = c("-Xmx4000m")) # Allocate a large amount of memory for the JVM, as described here: https://cdk-r.github.io/cdkr/articles/rcdk.html#parsing-smiles
library(rcdk, quietly = TRUE)
data(bpdata)

# Test rCDK
cat('Testing rcdk ...\n')
cat("Descriptor IDs:", desc_str, "\n")
cat("SMILES IDs:", smi_str, "\n")
smis <- bpdata$SMILES[smi_ids]
descs  <- get.desc.names("all")[desc_ids]
# The following descriptor produces a `Error: segfault from C stack overflow` error, at least for `OpenJDK Runtime Environment (build 11.0.23+9-post-Ubuntu-1ubuntu122.04.1)`
# > descs[20]
# [1] "org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor"
n <- length(smis)
X <- vector("list", n)
for (i in seq_along(smis)) {
    smi <- smis[[i]]
    mol <- parse.smiles(smi)
    if (progress) cat(sprintf("%d/%d %s\n", i, n, smi))
    X[[i]] <- eval.desc(mol, descs)
    rJava::.jcall("java/lang/System","V","gc") # Free memory as described here: https://cdk-r.github.io/cdkr/articles/rcdk.html#parsing-smiles
    gc()
}
X <- do.call(rbind, X)
exitcode <- if (nrow(X) == n && !all(is.na(X))) 0 else 1
ok <- "\033[32mok\033[39m\n" # green
failed <- "\033[31mfailed\033[39m\n" # red
cat("Test result:", if (exitcode == 0) ok else failed)
quit(status = exitcode, save = 'no')
