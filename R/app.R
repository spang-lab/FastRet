# Public #####

#' @export
#' @keywords public
#'
#' @title Start the FastRet GUI
#'
#' @description Starts the FastRet GUI
#'
#' @param port
#' The port the application should listen on
#'
#' @param host
#' The address the application should listen on
#'
#' @param reload
#' Whether to reload the application when the source code changes
#'
#' @param nw
#' The number of worker processes started. The first worker always listens for
#' user input from the GUI. The other workers are used for handling long running
#' tasks like model fitting or clustering. If `nw` is 1, the same process is
#' used for both tasks, which means that the GUI will become unresponsive during
#' long running tasks.
#'
#' @param nsw
#' The number of subworkers each worker is allowed to start. The higher this
#' number, the faster individual tasks like model fitting can be processed. A
#' value of 1 means that all subprocesses will run sequentially.
#'
#' @return
#' A shiny app. This function returns a shiny app that can be run to interact
#' with the model.
#'
#' @details
#' If you set `nw = 3` and `nsw = 4`, you should have at least 16 cores, one
#' core for the shiny main process. Three cores for the three worker processes
#' and 12 cores (3 * 4) for the subworkers. For the default case, `nworkers = 2`
#' and `nsw = 1`, you only need 3 cores, as `nsw = 1` means that all
#' subprocesses will run sequentially.
#'
#' @examples
#' if (interactive()) start_gui()
start_gui <- function(port = 8080,
                      host = "0.0.0.0",
                      reload = FALSE,
                      nw = 2,
                      nsw = 1) {
    # Use [start_gui_in_devmode()] for development
    catf("Checking CDK version")
    check_cdk_version()
    oldplan <- future::plan("multisession", workers = nw)
    on.exit(future::plan(oldplan), add = TRUE)
    catf("Starting FastRet GUI")
    app <- fastret_app(port, host, reload, nsw)
    shiny::runApp(app)
}

#' @export
#' @keywords public
#'
#' @title The FastRet GUI
#'
#' @description
#' Creates the FastRet GUI
#'
#' @param port
#' The port the application should listen on
#'
#' @param host
#' The address the application should listen on
#'
#' @param reload
#' Whether to reload the application when the source code changes
#'
#' @param nsw
#' The number of subworkers each worker is allowed to start. The higher this
#' number, the faster individual tasks like model fitting can be processed.
#'
#' @return An object of class `shiny.appobj`.
#'
#' @examples
#' x <- fastret_app()
#' if (interactive()) shiny::runApp(x)
fastret_app <- function(port = 8080,
                        host = "0.0.0.0",
                        reload = FALSE,
                        nsw = 1) {
    shiny::shinyApp(
        ui = function(req) fastret_ui(req),
        server = function(input, output, session) fastret_server(input, output, session, nsw),
        options = list(port = port, host = host, quiet = TRUE, launch.browser = FALSE, reload = reload),
        onStart = function() catf("Listening on http://localhost:%s", port)
    )
}

# Private #####

#' @noRd
#' @title Start the FastRet GUI in development mode
#'
#' @description Starts the FastRet GUI in development mode
#'
#' @param strategy
#' The strategy to use for parallel processing. Can be one of "sequential",
#' "multicore", "multisession"
#'
#' @param startMode
#' The start mode to use. Can be one of "Train new Model", "Predict Retention
#' Times", "Selective Measuring", "Adjust existing Model"
#'
#' @return
#' NULL. Called for side effects.
#'
#' @details
#' By using no subworkers and multicore or sequential, we can ensure that all
#' processes are forked from the current R session and therefore use the
#' functions loaded via devtools. If we use multisession and or subworkers,
#' these processes will use the installed version of FastRet instead. ==> If we
#' work on the UI part, we can use multisession and/or subworkers, because the
#' UI part is handled by the main process, BUT, If we develop train/predict/plot
#' functions, we must use multicore or sequential and NO subworkers! In
#' particular, to use `browser()` in these functions, we must use sequential.
start_gui_in_devmode <- function(strategy = "sequential",
                                 startMode = "Train new Model") {

    catf("Checking args")
    strategies <- c("sequential", "multicore", "multisession")
    startModes <- c("Train new Model", "Predict Retention Times", "Selective Measuring", "Adjust existing Model")
    startMode <- match.arg(startMode, startModes)
    strategy <- match.arg(strategy, strategies)

    catf("Reloading FastRet")
    devtools::load_all() # needs to be called once with updated function

    catf("Setting development options")
    withr::local_options(list(
        shiny.autoreload = TRUE,
        FastRet.UI.startMode = startMode,
        warn = 1
    ))

    catf("Initializing cluster")
    oldplan <- future::plan(strategy)
    on.exit(future::plan(oldplan), add = TRUE, after = FALSE)

    catf("Starting FastRet GUI in development mode")
    pkg_root <- dirname(system.file("DESCRIPTION", package = "FastRet"))
    shiny::with_devmode(TRUE, shiny::runApp(pkg_root), verbose = TRUE)
}

check_cdk_version <- function() {
    if (rcdk::cdk.version() < "2.9") {
        msg <- paste(
            sep = "\n",
            "FastRet requires CDK Version 2.9 or greater, but the installed version is %s.",
            "For details about the installation process see https://github.com/CDK-R/rcdklibs."
        )
        msgf <- sprintf(msg, rcdk::cdk.version())
        stop(msgf, call. = FALSE)
    }
}
