options(
    shiny.autoreload = TRUE,
    FastRet.mocks = c(
        # Shiny mocks
        "inpFRM",       # inits RV$ubInpFRM
        "inpDf",        # inits RV$inpDf
        "adjDf",        # inits RV$adjDf
        # "btnTrain",     # triggers SE$ABH$btnTrain
        "cluster_calc", # TODO
        "tiPredSmiles", # TODO
        # Functions mocks
        # "getCDs",              # mocks getCDs
        # "preprocess_data",     # mocks preprocess_data
        # "train_frm",             # mocks train_frm
        # "selective_measuring", # mocks selective_measuring
        NULL
    ),
    # FastRet.UI.startMode = "Train new Model",
    FastRet.UI.startMode = "Predict Retention Times",
    # FastRet.UI.startMode = "Selective Measuring",
    # FastRet.UI.startMode = "Adjust existing Model",
    warn = 1
)
catf("Initializing cluster")
# future::plan("sequential")
future::plan("multicore", workers = 3)
catf("Done")
fastret_app(port = 8080, host = "0.0.0.0", reload = FALSE, nsw = 1)
