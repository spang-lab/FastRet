library(testthat)

test_that("train_frm works if `method == \"GBTree\"`", {
    model1 <- train_frm(df=RP[1:20, ], method="gbtree", nfolds=2, nw=1, verbose=0, seed=42)
    model2 <- train_frm(df=RP[1:20, ], method="gbtree", nfolds=2, nw=2, verbose=0, seed=42)
    expect_equal(names(model1), c("model", "df", "cv", "seed", "version", "args"))
    model2$args$nw <- 1 # Set to same value as in model1 before comparison
    # The per-fold CV sub-models are trained in the main process when nw = 1 but
    # in parallel workers when nw = 2. Under pkgload::load_all() those workers
    # load the *installed* FastRet, so they stamp the installed package version
    # into $version instead of the checked-out one, and the comparison would fail
    # on that field alone (it matches under R CMD check, which installs the
    # package first). Normalise this environment-dependent field, as we do for nw.
    for (f in names(model2$cv$models)) {
        model2$cv$models[[f]]$version <- model1$cv$models[[f]]$version
    }
    expect_true(all.equal(model1, model2)) # (1)
    # (1) we can't use expect_equal because the xgboost objects contain pointers
    # to the underlying C structs which causes expect_equal to fail
})

test_that("train_frm grid-searches xgboost params when `tune_grid` is supplied", {
    grid <- expand.grid(max_depth = c(3, 4), eta = c(0.1, 0.2), gamma = 1,
                        colsample_bytree = 0.8, subsample = 0.7, min_child_weight = 3)
    m <- train_frm(RP[1:40, ], method = "brt", tune_grid = grid, nfolds = 3,
                   do_cv = FALSE, verbose = 0, seed = 1)
    expect_s3_class(m$model, "xgb.Booster")
    expect_length(predict(m, RP[1:40, ]), 40)
    # tune_grid is ignored (with a warning) for non-BRT methods
    expect_warning(train_frm(RP[1:40, ], method = "lasso", tune_grid = grid,
                             do_cv = FALSE, verbose = 0, seed = 1))
})

test_that("train_frm accepts a grid-size keyword for tune_grid", {
    m <- train_frm(RP[1:40, ], method = "brt", tune_grid = "small", nfolds = 3,
                   do_cv = FALSE, verbose = 0, seed = 1)
    expect_s3_class(m$model, "xgb.Booster")
    expect_length(predict(m, RP[1:40, ]), 40)
})
