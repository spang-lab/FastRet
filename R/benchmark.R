# r$> G_summary[order(G_summary$method, G_summary$do_cv), ]
# method do_cv rt_mean      rt_sd
# 1    BRT FALSE  0.9611 0.32503007
# 3    BRT  TRUE  2.7523 0.33611706
# 2  Lasso FALSE  0.8300 0.04980406
# 4  Lasso  TRUE  2.6365 0.19202445
benchmark_runtime <- function(seed = 1, cache = TRUE, do_cv = TRUE) {

    # Init Grid
    G <- expand.grid(
        method = c("Lasso", "BRT"),
        do_cv = c(TRUE, FALSE),
        seed = 1:100,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )

    # Benchmark
    n <- nrow(G)
    a <- Sys.time()
    for (i in seq_len(n)) {
        G[i, "rt"] <- measure_runtime(G$method[i], G$do_cv[i], G$seed[i])
        rttot <- as.numeric(Sys.time() - a, units = "mins")
        logf("%d/%d: rt=%.3f secs (total=%.1f mins)", i, n, G$rt[i], rttot)
    }
    G

    # Summarize
    S <- aggregate(rt ~ method + do_cv, data = G, FUN = function(x) c(mean = mean(x), sd = sd(x)))
    S <- do.call(data.frame, S)
    colnames(S) <- c("method", "do_cv", "rt_mean", "rt_sd")
    S[order(S$method, S$do_cv), ]
    S

}

measure_runtime <- function(method, do_cv, seed) {
    nw <- if (do_cv) 5 else 1
    runtime <- system.time(
        train_frm(RP, method, verbose=0, seed=seed, do_cv=do_cv, nw=5)
    )
    runtime[["elapsed"]]
}
