# Main (Public) #####

#' @noRd
#' @keywords internal
#'
#' @title Server function for the FastRet GUI
#'
#' @description
#' This function initializes the server-side of the FastRet GUI. It sets up
#' shiny::reactive values, initializes handlers for various events and tasks,
#' and sets up observers and outputs.
#'
#' @param input
#' List of input values from the Shiny application.
#'
#' @param output
#' List of output values to be sent to the Shiny application.
#'
#' @param session
#' The session object passed to function given to shinyServer.
#'
#' @param nsw
#' The number of subprocesses each worker is allowed to start. The higher this
#' number, the faster individual tasks like model fitting can be processed.
#' Providing a value of 1 will disable parallel processing.
#'
#' @return
#' No return value. The function is used for its side effect of setting up the
#' server-side of the Shiny application.
#'
fastret_server <- function(input, output, session, nsw = 1) {
    catf("Start: fastret_server (session$token == %s)", session$token)
    if (!is.numeric(nsw) || nsw < 1) stop("nsw must be a positive integer")
    SE <- environment()
    SE$RV <- shiny::reactiveValues()
    init_log_dir(SE)
    init_extended_tasks(SE)
    init_extended_task_handlers(SE)
    init_reactives(SE)
    init_action_button_handlers(SE)
    init_input_handlers(SE)
    init_download_handlers(SE)
    init_upload_handlers(SE)
    init_special_event_handlers(SE)
    init_observers(SE)
    init_table_output_handler(SE)
    init_plot_output_handler(SE)
    init_outputs(SE)
    init_download_handlers(SE)
    catf("Exit: fastret_server")
}

# Inits (Private) #####

init_extended_tasks <- function(SE) {
    catf("Start: init_extended_tasks")
    btns <-  c("btnTrain",   "btnSM",               "btnPred",     "btnAdj")
    funcs <- c("train_frm",  "selective_measuring", "predict", "adjust_frm")
    logfiles <- file.path(SE$logdir, paste0(btns, ".log"))
    lapply(seq_along(btns), function(i) {
        SE$ET[[btns[i]]] <- extendedTask(func = funcs[i], logfile = logfiles[i])
    })
    # catf("Exit: init_extended_tasks")
}

init_extended_task_handlers <- function(SE) {
    catf("Start: init_extended_task_handlers")
    SE$ETH <- list()
    SE$ETH$btnTrain <- extendedTaskHandler("btnTrain",
        onSuccess = function(SE) {
            catf("Upating: SE$RV$trainedFRM, SE$RV$tblTrainResults")
            frm <- SE$ET$btnTrain$result()
            cds <- frm$df[!colnames(frm$df) %in% c("NAME", "RT", "SMILES")]
            cds <- round(cds, 2)
            RT <- round(frm$df$RT, 2)
            NAME <- frm$df$NAME
            SMILES <- frm$df$SMILES
            RT_PREDICTED <- round(predict(frm), 2)
            RT_PREDICTED_CV <- if (!is.null(frm$cv)) round(frm$cv$preds, 2) else rep(NA, nrow(frm$df))
            SE$RV$trainedFRM <- frm
            SE$RV$tblTrainResults <- cbind(RT, RT_PREDICTED, RT_PREDICTED_CV, NAME, SMILES, cds)
            catf("Showing buttons: dbSaveModel, dbSavePredictorSet")
            shinyjs::show("dbSaveModel")
            shinyjs::show("dbSavePredictorSet")
        },
        onRunning = function(SE) {
            catf("Hiding: dbSavePredictorSet, dbSaveModel")
            shinyjs::hide("dbSavePredictorSet")
            shinyjs::hide("dbSaveModel")
            # catf("Showing: abShowLogs") # TODO: enable logging output
            # shinyjs::show("abShowLogs")
            catf("Setting to NULL: SE$RV$trainedFRM, SE$RV$tblTrainResults")
            SE$RV$trainedFRM <- NULL
            SE$RV$tblTrainResults <- NULL
        }
    )
    SE$ETH$btnSM <- extendedTaskHandler("btnSM",
        onSuccess = function(SE) {
            smobj <- SE$ET$btnSM$result()
            mtbl <- smobj$clustering[smobj$clustering$IS_MEDOID, c("RT", "NAME", "CLUSTER", "SMILES")]
            mtbl <- mtbl[order(mtbl$CLUSTER), ]
            mtbl <- `colnames<-`(mtbl, c("RT", "NAME", "MEDOID", "SMILES"))
            ctbl <- smobj$clustering[, c("RT", "NAME", "CLUSTER", "IS_MEDOID", "SMILES")]
            SE$RV$cluster_calc <- smobj
            SE$RV$tblMedoids <- mtbl
            SE$RV$tblClustering <- ctbl
            shinyjs::show("dbSaveCluster")
        }
    )
    SE$ETH$btnPred <- extendedTaskHandler("btnPred",
        onSuccess = function(SE) {
            catf("Updating: SE$RV$tblPredResults")
            df <- SE$R$predDfCombined()
            df$PREDICTED_RT <- sprintf("%.2f", SE$ET$btnPred$result())
            SE$RV$tblPredResults <- df
            shinyjs::show("dbSavePred")
        },
        onRunning = function(SE) {
            catf("Hiding button: dbSavePred")
            shinyjs::hide("dbSavePred")
        }
    )
    SE$ETH$btnAdj <- extendedTaskHandler("btnAdj",
        onSuccess = function(SE) {
            shinyjs::show("dbSaveAdjModel")
        },
        onRunning = function(SE) {
            shinyjs::hide("dbSaveAdjModel")
        }
    )
    # catf("Exit: init_extended_task_handlers")
}

init_reactives <- function(SE) {
    SE$R <- list()
    SE$R$predDfCombined <- shiny::reactive({
        smiles <- SE$RV$predSmiles # data.frame(NAME = "Input SMILES", SMILES = smiles)
        df <- SE$RV$predDf # data.frame(NAME = names, SMILES = smiles)
        if (is.null(smiles) && is.null(df)) return(NULL)
        if (is.null(smiles)) return(df)
        if (is.null(df)) return(smiles)
        row1 <- structure(rep(NA, ncol(df)), .Names = colnames(df))
        row1[names(smiles)] <- smiles
        dfCombined <- rbind(row1, df)
        return(dfCombined)
    })
}

init_action_button_handlers <- function(SE) {
    catf("Start: init_action_button_handlers")
    SE$ABH <- list()
    SE$ABH$btnTrain <- function(SE) {
        if (is.null(SE$RV$inpDf)) stop("Please upload a excel sheet with the required data first")
        SE$ET$btnTrain$invoke( # takes same argument as [train_frm()]
            df = SE$RV$inpDf,
            method = c("lasso", "gbtree")[as.numeric(SE$input$rbMethod)],
            verbose = 1,
            nw = SE$nsw
        )
    }
    SE$ABH$btnPred <- function(SE) {
        frm <- SE$RV$inpFRM
        df <- SE$R$predDfCombined() # data.frame(NAME, SMILES)
        if (is.null(frm) || is.null(frm$model)) stop("Please upload a valid model first")
        if (is.null(df) || nrow(df) == 0) stop("Please enter a valid SMILES string first or upload a list of SMILES as xlsx")
        SE$ET$btnPred$invoke( # takes same argument as [predict()]
            object = frm,
            df = df,
            verbose = 1,
            adjust = NULL # i.e. adjust predictions if `object$adj` is not NULL
        )
    }
    SE$ABH$btnSM <- function(SE) {
        if (is.null(SE$RV$inpDf)) stop("Please upload a excel sheet with the required data first")
        SE$ET$btnSM$invoke( # takes same argument as [selective_measuring()]
            SE$RV$inpDf,
            k_cluster = SE$input$niK,
            verbose = 1
        )
    }
    SE$ABH$btnAdj <- function(SE) {
        if (!inherits(SE$RV$inpFRM, "frm")) stop("Please upload a valid, pretrained model first")
        if (is.null(SE$RV$adjDf)) stop("Please upload valid data for prediction adjustment first")
        SE$ET$btnAdj$invoke( # takes same argument as [adjust_frm()]
            frm = SE$RV$inpFRM,
            new_data = SE$RV$adjDf,
            predictors = as.numeric(SE$input$ciPredictors)
        )
    }
    # catf("Exit: init_action_button_handlers")
}

init_input_handlers <- function(SE) {
    catf("Start: init_input_handlers")
    SE$IPH <- list()
    SE$IPH$tiPredSmiles <- function(SE) {
        tryCatch({
            smiles <- SE$input$tiPredSmiles
            if (is.null(smiles) || nchar(trimws(smiles)) == 0) {
                catf("tiPredSmiles is empty. Setting SE$RV$predSmiles to NULL")
                SE$RV$predSmiles <- NULL
                SE$output$toPredSmilesError <- NULL
            } else {
                catf("Validating SMILES string")
                # Validate SMILES by trying to parse it with rCDK
                tmp <- try(silent = TRUE, {
                    obj <- rcdk::parse.smiles(smiles)[[1]]
                    rcdk::convert.implicit.to.explicit(obj)
                    TRUE
                })
                if (inherits(tmp, "try-error")) {
                    catf("Validation failed. Displaying 'Error: SMILES string is invalid'")
                    SE$RV$predSmiles <- NULL
                    SE$output$toPredSmilesError <- shiny::renderText("Error: SMILES string is invalid")
                } else {
                    catf("Validation successful. Updating SE$RV$predSmiles")
                    SE$RV$predSmiles <- data.frame(NAME = "Input SMILES", SMILES = smiles)
                    SE$output$toPredSmilesError <- NULL
                }
            }
        },
        error = function(e) {
            catf("Validation failed. Displaying error message")
            SE$RV$predSmiles <- NULL
            SE$output$toPredXlsxError <- shiny::renderText(sprintf("Error: %s\n", e$message))
        })
    }
    # catf("Exit: init_input_handlers")
}

init_download_handlers <- function(SE) {
    catf("Start: init_download_handlers")
    SE$DLH <- list()
    SE$DLH$dbSavePredictorSet <- shiny::downloadHandler(
        filename = function() {
            paste("predictor_set_", Sys.Date(), ".xlsx", sep = "")
        },
        content = function(file) {
            frm <- SE$RV$trainedFRM
            openxlsx::write.xlsx(frm$df, file, rowNames = FALSE)
        }
    )
    SE$DLH$dbSaveModel <- shiny::downloadHandler(
        filename = function() {
            frm <- SE$RV$trainedFRM
            mtype <- if (inherits(frm$model, "xgb.Booster")) "xgboost" else "lasso"
            sprintf("fastret-%s-model-%s.rds", mtype, Sys.Date())
        },
        content = function(file) {
            saveRDS(SE$RV$trainedFRM, file)
        }
    )
    SE$DLH$dbSaveCluster <- shiny::downloadHandler(
        filename = function() {
            sprintf("fastret-k-%s-clustering.xlsx", SE$input$niK)
        },
        content = function(file) {
            openxlsx::write.xlsx(SE$RV$cluster_calc$clustering, file, rowNames = FALSE)
        }
    )
    SE$DLH$dbSavePred <- shiny::downloadHandler(
        filename = function() {
            "predictions.xlsx"
        },
        content = function(file) {
            openxlsx::write.xlsx(SE$RV$tblPredResults, file, rowNames = FALSE)
        }
    )
    SE$DLH$dbSaveAdjModel <- shiny::downloadHandler(
        filename = function() {
            frm <- SE$ET$btnAdj$result()
            mtype <- if (inherits(frm$model, "xgb.Booster")) "xgboost" else "lasso"
            sprintf("fastret-%s-model-adjusted-%s.rds", mtype, Sys.Date())
        },
        content = function(file) {
            saveRDS(SE$ET$btnAdj$result(), file)
        }
    )
    # catf("Exit: init_download_handlers")
}

init_upload_handlers <- function(SE) {
    catf("Start: init_upload_handlers")
    SE$ULH <- list()
    SE$ULH$ubInpFRM <- function(SE) {
        tryCatch({
            rds <- SE$input$ubInpFRM$datapath
            catf("Reading and validating %s", rds)
            ubInpFRM <- readRDS(rds)
            ubInpFRM <- validate_inputmodel(ubInpFRM)
            catf("Validation successful. Updating: SE$RV$inpFRM and SE$output$toInpFRMError.")
            SE$RV$inpFRM <- ubInpFRM
            SE$output$toInpFRMError <- shiny::renderText("")
        },
        error = function(e) {
            catf("Validation failed. Updating: SE$output$toInpFRMError and SE$RV$inpFRM.")
            SE$RV$inpFRM <- NULL
            SE$output$toInpFRMError <- shiny::renderText(paste("Error:", e$message))
        })
    }
    SE$ULH$ubInpXlsx <- function(SE) {
        tryCatch({
            xlsx <- SE$input$ubInpXlsx$datapath
            catf("Reading and validating %s", xlsx)
            inpDf <- openxlsx::read.xlsx(xlsx, sheet = 1)
            inpDf <- validate_inputdata(inpDf, min_cds = 0, stop_on_unknown = FALSE)
            catf("Validation successful. Updating: SE$RV$inpDf and SE$output$toInpXlsxError.")
            SE$RV$inpDf <- inpDf
            SE$output$toInpXlsxError <- shiny::renderText("")
        },
        error = function(e) {
            catf("Validation failed. Updating: SE$output$toInpXlsxError and SE$RV$inpDf.")
            SE$RV$inpDf <- NULL
            SE$output$toInpXlsxError <- shiny::renderText(paste("Error:", e$message))
        })
    }
    SE$ULH$ubPredXlsx <- function(SE) {
        tryCatch({
            xlsx <- SE$input$ubPredXlsx$datapath
            catf("Reading and validating %s", xlsx)
            pred_df <- openxlsx::read.xlsx(xlsx, sheet = 1)
            pred_df <- validate_inputdata(pred_df, require = c("NAME", "SMILES"), min_cds = 0, stop_on_unknown = FALSE)
            catf("Validation successful. Updating: SE$RV$predDf and SE$output$toPredXlsxError.")
            SE$RV$predDf <- pred_df
            SE$output$toPredXlsxError <- NULL
        },
        error = function(e) {
            catf("Validation failed. Updating: SE$output$toPredXlsxError and SE$RV$pred_df.")
            SE$RV$predDf <- NULL
            SE$output$toPredXlsxError <- shiny::renderText(sprintf("Error: %s\n", e$message))
        })
    }
    SE$ULH$ubAdjXlsx <- function(SE) {
        tryCatch({
            xlsx <- SE$input$ubAdjXlsx$datapath
            catf("Reading and validating %s", xlsx)
            adjDf <- openxlsx::read.xlsx(xlsx, sheet = 1)
            adjDf <- validate_inputdata(adjDf, min_cds = 0, stop_on_unknown = FALSE)
            catf("Validation successful. Updating: SE$RV$adjDf and SE$output$toAdjXlsxError.")
            SE$RV$adjDf <- adjDf
            SE$output$toAdjXlsxError <- shiny::renderText("")
        },
        error = function(e) {
            catf("Validation failed. Updating: SE$output$toAdjXlsxError and SE$RV$adjDf.")
            SE$RV$adjDf <- NULL
            SE$output$toAdjXlsxError <- shiny::renderText(paste("Error:", e$message))
        })
    }
    # catf("Exit: init_upload_handlers")
}

init_special_event_handlers <- function(SE) {
    catf("Start: init_special_event_handlers")
    SE$SEH <- list()
    SE$SEH$SessionEnded <- function(SE) {
        catf("Ended session %s", SE$session$token)
    }
    # catf("Exit: init_special_event_handlers")
}

init_observers <- function(SE) {
    catf("Start: init_observers")

    # Special Events
    shinyhelper::observe_helpers()
    shiny::onSessionEnded(function() SE$SEH$SessionEnded(SE))

    # Upload Button Handler
    lapply(names(SE$ULH), function(x) {
        shiny::observeEvent(SE$input[[x]], {
            catf("Start: SE$ULH$%s", x)
            SE$ULH[[x]](SE)
            catf("Exit: SE$ULH$%s", x)
        })
    })

    # Action Button Handler
    lapply(names(SE$ABH), function(x) {
        shiny::observeEvent(SE$input[[x]], {
            catf("Start: SE$ABH$%s", x)
            withShowError(SE$ABH[[x]](SE))
            catf("Exit: SE$ABH$%s", x)
        })
    })

    # Input Widget Handler
    lapply(names(SE$IPH), function(x) {
        shiny::observeEvent(SE$input[[x]], {
            catf("Start: SE$IPH$%s", x)
            SE$IPH[[x]](SE)
            catf("Exit: SE$IPH$%s", x)
        })
    })

    # Extended Task Handler
    lapply(names(SE$ETH), function(x) {
        shiny::observeEvent(
            eventExpr = SE$ET[[x]]$status(),
            handlerExpr = {
                catf("Start: SE$ETH$%s", x)
                SE$ETH[[x]](SE)
                catf("Exit: SE$ETH$%s", x)
            },
            ignoreInit = TRUE
        )
        bslib::bind_task_button(SE$ET[[x]], x) # (1)
        # (1) Binds the given extended task to the input_task_button with ID
        # `x`, i.e. it gets disabled while the extended task is running. Note:
        # for the above to work, the extended task, the extended task handler
        # and the corresponding input_task_button must use the same ID.
    })

    # Download Buttons Handler
    # (must be stored inside `output` because that's where
    # `shiny::downloadButton` looks for them)
    lapply(names(SE$DLH), function(x) {
        SE$output[[x]] <- SE$DLH[[x]]
    })

    catf("Exit: init_observers")
}

init_plot_output_handler <- function(SE) {
    catf("Start: init_plot_output_handler")
    SE$POH <- list()
    SE$POH$poTrainPerfCV <- function(SE) {
        catf("Start: SE$POH$poTrainPerfCV")
        frm <- SE$RV$trainedFRM
        if (!is.null(frm) && !is.null(frm$cv)) {
            catf("Rendering: SE$output$poTrainPerfCV")
            plot_frm(frm, type = "scatter.cv")
        } else {
            catf("FRM is NULL. Clearing SE$output$poTrainPerfCV")
            NULL
        }
        catf("End: SE$POH$poTrainPerfCV")
    }
    SE$POH$poTrainPerf <- function(SE) {
        catf("Start: SE$POH$poTrainPerf")
        frm <- SE$RV$trainedFRM
        if (!is.null(frm)) {
            catf("Rendering: SE$output$poTrainPerf")
            plot_frm(frm, type = "scatter.train")
        } else {
            catf("FRM is NULL. Clearing SE$output$poTrainPerf")
            NULL
        }
        catf("End: SE$POH$poTrainPerf")
    }
    SE$POH$poAdjPerfCV <- function(SE) {
        catf("Start: SE$POH$poAdjPerfCV")
        if (SE$ET$btnAdj$status() != "error") {
            catf("Rendering: SE$output$poTrainPerf")
            plot_frm(SE$ET$btnAdj$result(), type = "scatter.cv.adj")
        } else {
            catf("FRM is NULL. Clearing SE$output$poTrainPerf")
            NULL
        }
    }
    SE$POH$poAdjPerf <- function(SE) {
        catf("Start: SE$POH$poAdjPerf")
        if (SE$ET$btnAdj$status() != "error") {
            catf("Rendering: SE$output$poAdjPerf")
            plot_frm(SE$ET$btnAdj$result(), type = "scatter.train.adj")
        } else {
            catf("FRM is NULL. Clearing SE$output$poAdjPerf")
            NULL
        }
        catf("End: SE$POH$poAdjPerf")
    }
    # catf("Exit: init_plot_output_handler")
}

init_table_output_handler <- function(SE) {
    catf("Start: init_table_output_handler")
    SE$TBH <- list()
    SE$TBH$tblTrainResults <- function(SE) {
        if (is.null(SE$RV$tblTrainResults)) {
            catf("Clearing: SE$output$tblTrainResults")
            SE$output$tblTrainResults <- NULL
        } else {
            catf("Rendering: SE$output$tblTrainResults")
            SE$output$tblTrainResults <- renderTbl(SE$RV$tblTrainResults)
        }
    }
    SE$TBH$tblPredResults <- function(SE) {
        if (is.null(SE$RV$tblPredResults)) {
            catf("Clearing: SE$output$tblPredResults")
            SE$output$tblPredResults <- renderTbl(NULL)
        } else {
            catf("Rendering: SE$output$tblPredResults")
            SE$output$tblPredResults <- renderTbl(SE$RV$tblPredResults)
        }
    }
    SE$TBH$tblMedoids <- function(SE) {
        if (is.null(SE$RV$tblMedoids)) {
            catf("Clearing: SE$output$tblMedoids")
            SE$output$tblMedoids <- NULL
        } else {
            catf("Rendering: SE$output$tblMedoids")
            SE$output$tblMedoids <- renderTbl(SE$RV$tblMedoids)
        }
    }
    SE$TBH$tblClustering <- function(SE) {
        if (is.null(SE$RV$tblClustering)) {
            catf("Clearing: SE$output$tblClustering")
            SE$output$tblClustering <- NULL
        } else {
            catf("Rendering: SE$output$tblClustering")
            SE$output$tblClustering <- renderTbl(SE$RV$tblClustering)
        }
    }
    # catf("Exit: init_table_output_handler")
}

init_outputs <- function(SE) {
    catf("Start: init_outputs")
    # UI outputs
    SE$output$ui_train_results   <- shiny::renderUI(ui_train_results(SE))
    SE$output$ui_sm_results      <- shiny::renderUI(ui_sm_results(SE))
    SE$output$ui_predict_results <- shiny::renderUI(ui_predict_results(SE))
    SE$output$ui_adjust_results  <- shiny::renderUI(ui_adjust_results(SE))
    # Plot Outputs
    SE$output$poTrainPerfCV <- shiny::renderPlot(SE$POH$poTrainPerfCV(SE), execOnResize = TRUE)
    SE$output$poTrainPerf   <- shiny::renderPlot(SE$POH$poTrainPerf(SE), execOnResize = TRUE)
    SE$output$poAdjPerfCV   <- shiny::renderPlot(SE$POH$poAdjPerfCV(SE), execOnResize = TRUE)
    SE$output$poAdjPerf     <- shiny::renderPlot(SE$POH$poAdjPerf(SE), execOnResize = TRUE)
    # Table Outputs
    shiny::observe(SE$TBH$tblPredResults(SE))
    shiny::observe(SE$TBH$tblTrainResults(SE))
    shiny::observe(SE$TBH$tblMedoids(SE))
    shiny::observe(SE$TBH$tblClustering(SE))
    # Text Outputs
    btnIDs <- c("btnTrain", "btnSM", "btnPred", "btnAdj")
    vtoIDs <- c("vtoTrainLogs", "vtoSMLogs", "vtoPredLogs", "vtoAdjLogs")
    lapply(seq_along(btnIDs), function(i) {
        btn <- btnIDs[i]
        vto <- vtoIDs[i]
        SE$output[[vto]] <- shiny::renderText({
            shiny::invalidateLater(1000)
            logfile <- file.path(SE$logdir, paste0(btn, ".log"))
            if (file.exists(logfile)) {
                paste(readLines(logfile), collapse = "\n")
            } else {
                NULL
            }
        })
    })
    # catf("Exit: init_outputs")
}

# Inits (Private) #####


# Helpers (Public) #####

#' @export
#' @keywords internal
#'
#' @title Execute an expression with a timeout
#'
#' @param expr
#' The expression to execute
#'
#' @param timeout
#' The timeout in seconds. Default is 2.
#'
#' @return
#' The result of the expression
#'
#' @examples
#' withTimeout(
#'      cat("This works\n"),
#'      timeout = 0.2
#' )
#' try(silent = TRUE, withTimeout(
#'     expr = {Sys.sleep(0.2); cat("This fails\n")},
#'     timeout = 0.1
#' ))
#'
withTimeout <- function(expr, timeout = 2) {
    setTimeLimit(cpu = timeout, elapsed = timeout, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE, after = FALSE)
    expr
}

#' @export
#' @keywords internal
#'
#' @title Execute an expression while redirecting output to a file
#'
#' @param expr
#' The expression to execute
#'
#' @param logfile
#' The file to redirect output to. Default is "tmp.txt".
#'
#' @return The result of the expression
#'
#' @examples
#' logfile <- tempfile(fileext = ".txt")
#' withSink(logfile = logfile, expr = {
#'   cat("Helloworld\n")
#'   message("Goodbye")
#' })
#' readLines(logfile) == c("Helloworld", "Goodbye")
#'
withSink <- function(expr, logfile = tempfile(fileext = ".txt")) {
    zz <- file(logfile, open = "wt")
    on.exit(close(zz), add = TRUE, after = FALSE)
    sink(zz)
    on.exit(sink(), add = TRUE, after = FALSE)
    sink(zz, type = "message")
    on.exit(sink(type = "message"), add = TRUE, after = FALSE)
    expr
}

#' @export
#' @keywords internal
#'
#' @title Try expression with predefined error message
#'
#' @description
#' Executes an expression and prints an error message if it fails
#'
#' @param expr
#' The expression to execute
#'
#' @return
#' The result of the expression
#'
#' @examples
#' f <- function(expr) {
#'   val <- try(expr, silent = TRUE)
#'   err <- if (inherits(val, "try-error")) attr(val, "condition") else NULL
#'   if (!is.null(err)) value <- NULL
#'   list(value = val, error = err)
#' }
#' ret <- f(log("a")) # this error will not show up in the console
#' ret <- f(withStopMessage(log("a"))) # this error will show up in the console
#'
withStopMessage <- function(expr) {
    tryCatch(expr, error = function(e) {
        message("Error in ", deparse(e$call), " : ", e$message)
        stop(e)
    })
}

#' @export
#' @keywords internal
#'
#' @title Initialize log directory
#'
#' @description
#' Initializes the log directory for the session. It creates a new directory if
#' it does not exist.
#'
#' @param SE
#' A list containing session information.
#'
#' @return
#' Updates the logdir element in the SE list with the path to the log directory.
#'
#' @examples
#' SE <- as.environment(list(session = list(token = "asdf")))
#' init_log_dir(SE)
#' dir.exists(SE$logdir)
#'
init_log_dir <- function(SE) {
    catf("Start: init_log_dir")
    token <- SE$session$token
    logdir <- file.path(tempdir(), "FastRet", token)
    if (!dir.exists(logdir)) dir.create(logdir, recursive = TRUE)
    catf("Logdir: %s", logdir)
    SE$logdir <- logdir
}

#' @export
#' @keywords internal
#'
#' @title Add line end
#'
#' @description
#' Checks if a string ends with a newline character. If not, a newline character is appended.
#'
#' @param x
#' A string.
#'
#' @return
#' The input string with a newline character at the end if it was not already present.
#'
#' @examples
#' cat(withLineEnd("Hello"))
#'
withLineEnd <- function(x) {
    if (!grepl("\n$", x)) paste0(x, "\n") else x
}

# Helpers (Private) #####

withShowError <- function(expr, error = NULL) {
    tryCatch(
        expr,
        error = function(e) showError(e$message)
    )
}

showError <- function(msg = NULL, expr = NULL, duration = 10) {
    if (is.null(msg)) {
        msg <- tryCatch(expr, error = function(e) e$message)
    }
    catf("Displaying error message: %s", msg)
    shiny::showNotification(msg, type = "error", duration = duration)
}

#' @noRd
#' @title Create an shiny::ExtendedTask Object
#'
#' @description
#' This function wraps a given function in a [promises::future_promise()] and
#' the result into a [shiny::ExtendedTask()] object. When the
#' shiny::ExtendedTask Object is invoked, the function is executed
#' asynchronously in a seperate process (assuming [future::plan()] has been
#' called with strategy unequal "sequential"). Normal output, messages, warnings
#' and errors from that process get redirected to `logfile`. The status of the
#' task can be checked via the `status()` method. As soon as `status()` returns
#' `"success"`, the result can be retrieved via the `result()` method. If an
#' error has occured, `status()` will return `"error"`. In this case, calling
#' `result()` will reraise the error that occured while executing the task.
#' Querying the status or value of the task requires a shiny::reactive context,
#' e.g. via [shiny::reactive()], [shiny::observe()] or [shiny::reactiveConsole].
#'
#' @param func
#' A function that accepts any number of arguments and returns a value.
#'
#' @return
#' An shiny::ExtendedTask object that wraps the provided function. For further
#' details see: [shiny::ExtendedTask()].
#'
#' @examples
#' shiny::reactiveConsole(enabled = TRUE)
#' on.exit(shiny::reactiveConsole(enabled = FALSE), add = TRUE)
#'
#' f <- function(x) log(x)
#' logfile <- tempfile(fileext = ".log")
#' et <- extendedTask("f", logfile)
#'
#' et$status() == "initial"
#' et$invoke(x = 1)
#' et$status() == "success"
#' et$result() == log(1)
#'
#' et$invoke(x = -1)
#' et$status() == "success"
#' is.na(et$result()) == TRUE
#'
#' et$invoke(x = "a")
#' et$status() == "error"
#' x <- try(et$result(), silent = TRUE)
#' attr(x, "condition")$message == "non-numeric argument to mathematical function"
#'
#' g <- function(x) Sys.sleep(0.02)
#' logfile <- tempfile(fileext = ".log")
#' et <- extendedTask("g", logfile)
#' et$invoke()
#' et$status() == "running"
#' Sys.sleep(0.04)
#' et$status() == "success"
extendedTask <- function(func, logfile = tempfile(fileext = ".log"), timeout = 300) {
    logfile <- logfile
    func <- as.symbol(func)
    langobj <- substitute(
        shiny::ExtendedTask$new(function(...) {
            promises::future_promise(
                seed = TRUE,
                conditions = NULL,
                stdout = NA,
                package = "FastRet",
                expr = {
                    withLineEnd; # make sure this is exported by the future package
                    opts <- options(FastRet.catf.prefix = function() now("%H:%M:%OS2 "))
                    on.exit(options(opts), add = TRUE, after = FALSE)
                    withSink(logfile = logfile, withCallingHandlers(
                        withTimeout(timeout = timeout, func(...)),
                        message = function(m) {
                            cat("Message:", withLineEnd(m$message))
                        },
                        warning = function(w) {
                            cat("Warning:", withLineEnd(w$message))
                            invokeRestart("muffleWarning")
                        },
                        error = function(e)   {
                            cat("Error:", withLineEnd(e$message))
                            # stop(e)
                        }
                    ))
                }
            )
        }) # evaluate func and logfile
    )
    ET <- eval(langobj)
    ET
}

#' @noRd
#' @title Create an shiny::ExtendedTask Handler
#'
#' @description
#' This function creates a handler for an shiny::ExtendedTask object. The
#' handler checks the status of the task and executes the appropriate function
#' based on the status. The status can be "error", "running", "success", or
#' "initial".
#'
#' @param id
#' The ID of the shiny::ExtendedTask object. This ID must be unique and must
#' match the ID of an shiny::ExtendedTask object created via
#' `init_extended_tasks()`.
#'
#' @param onSuccess
#' A function that is executed when the task completes successfully. This
#' function accepts a single argument `SE`, which must point to the environment
#' of the corresponding shiny server function. I.e. inside `server` you should
#' call `SE <- environment()` and pass `SE` to this function.
#'
#' @param onRunning
#' A function that is executed when the task is still running. This function
#' also accepts the session environment `SE` as an argument.
#'
#' @param onError
#' A function that is executed when the task encounters an error. This function
#' also accepts the session environment `SE` as an argument.
#'
#' @param displayError
#' A boolean value that determines whether to display an error message to the
#' user when the task encounters an error. The default value is TRUE.
#'
#' @return
#' A function that checks the status of the shiny::ExtendedTask object and
#' executes the appropriate function based on the status.
#'
#' @examples
#' logfile <- tempfile(fileext = ".log")
#' f <- function(x) log(x)
#' SE <- list(ET = list(task1 = extendedTask(f, logfile)))
#' f <- function(SE) print("Task completed successfully!")
#' extendedTaskHandler(id = "task1", onSuccess = f)
#'
extendedTaskHandler <- function(id,
                                onSuccess = function(SE) {},
                                onRunning = function(SE) {},
                                onError = function(SE) {},
                                displayError = TRUE
                                ) {
    name <- sprintf("SE$ETH$%s", id)
    function(SE) {
        if (SE$ET[[id]]$status() == "error") {
            catf("Task failed.")
            if (isTRUE(displayError)) {
                catf("Displaying error message to user.")
                showError(expr = SE$ET[[id]]$result())
            }
            catf("Executing onError handler.")
            return(onError(SE))
        } else if (SE$ET[[id]]$status() == "running") {
            catf("Task is still running. Executing onRunning handler.")
            return(onRunning(SE))
        } else if (SE$ET[[id]]$status() == "success") {
            catf("Task completed successfully. Executing onSuccess handler.")
            return(onSuccess(SE))
        } else if (SE$ET[[id]]$status() == "initial") {
            catf("Task has not been started yet. Doing nothing.")
            return(NULL)
        } else {
            stop("Unknown status of shiny::ExtendedTask object")
        }
    }
}

renderTbl <- function(expr,
                      rownames = FALSE,
                      pageLength = 5,
                      scrollX = TRUE) {
    opts <- list(
        pageLength = pageLength,
        scrollX = scrollX
    )
    DT::renderDT(expr = expr, rownames = rownames, options = opts)
}

validate_inputdata <- function(df,
                               require = c("RT", "SMILES", "NAME"),
                               min_cds = 1,
                               stop_on_unknown = TRUE) {
    missing_cols <- setdiff(require, colnames(df))
    if (length(missing_cols) > 0) stop(sprintf("missing columns: %s", paste(missing_cols, collapse = ", ")))
    n_cds <- sum(colnames(df) %in% CDFeatures)
    if (n_cds < min_cds) {
        msg <- sprintf("At least %d chemical descriptors are required, but only %d are present", min_cds, n_cds)
        stop(msg)
    }
    unnown_cols <- setdiff(colnames(df), c("RT", "SMILES", "NAME", CDFeatures))
    if (stop_on_unknown && length(unnown_cols) > 0) {
        msg <- sprintf("Unknown columns present: %s", paste(unnown_cols, collapse = ", "))
        stop(msg)
    }
    invisible(df)
}

validate_inputmodel <- function(model) {
    model_nams <- names(model)
    expected_names <- c("model", "df", "cv")
    n_missing <- sum(!expected_names %in% model_nams)
    if (n_missing > 0) {
        if (n_missing < length(expected_names)) {
            missing <- paste(setdiff(expected_names, model_nams), collapse = ", ")
            errmsg1 <- sprintf("Model object is missing required elements: %s.", missing)
        } else {
            errmsg1 <- sprintf("Model object is invalid.")
        }
        errmsg2 <- sprintf("Please upload a model trained with FastRet version %s or greater.", packageVersion("FastRet"))
        errmsg <- paste(errmsg1, errmsg2)
        stop(errmsg)
    }
    invisible(model)
}
