# Redirect the on-disk chemical-descriptor cache to a temporary directory for the
# whole test suite and delete it once all tests have finished, so the package
# never writes to (nor leaves files in) the user's real cache during testing or
# `R CMD check`. The option covers the main process; the environment variable is
# inherited by parallel cluster / future worker processes spawned during tests.
local({
    cache_dir <- file.path(tempdir(), "FastRet-test-cache")
    options(FastRet.cache_dir = cache_dir)
    Sys.setenv(FASTRET_CACHE_DIR = cache_dir)
    withr::defer(unlink(cache_dir, recursive = TRUE, force = TRUE), teardown_env())
})
