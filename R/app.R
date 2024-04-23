#' @title Start the FastRet GUI
#' @description Starts the FastRet GUI
#' @param port The port the application should listen on
#' @param host The address the application should listen on
#' @param reload Whether to reload the application when the source code changes
#' @param nw The number of worker processes started. The first worker always listens for user input from the GUI. The other workers are used for handling long running tasks like model fitting or clustering. If `nw` is 1, the same process is used for both tasks, which means that the GUI will become unresponsive during long running tasks.
#' @param nsw The number of subworkers each worker is allowed to start. The higher this number, the faster individual tasks like model fitting can be processed. A value of 1 means that all subprocesses will run sequentially.
#' @return A shiny app. This function returns a shiny app that can be run to interact with the model.
#' @details If you set `nw = 3` and `nsw = 4`, you should have at least 16 cores, one core for the shiny main process. Three cores for the three worker processes. And 12 cores (3 * 4) for the subworkers. For the default case, `nworkers = 1` and `nsw = 2`, you should have at least 4 cores.
#' @keywords public
#' @export
start_gui <- function(port = 8080,
                      host = "0.0.0.0",
                      reload = FALSE,
                      nw = 2,
                      nsw = 1) {
    catf("Checking CDK version")
    check_cdk_version()
    oldplan <- future::plan("multisession", workers = nw)
    on.exit(future::plan(oldplan), add = TRUE)
    catf("Starting FastRet GUI")
    app <- fastret_app(port, host, reload, nsw)
    runApp(app)
}

start_gui_in_devmode <- function() {
    patch_shiny()
    patch_pkgload()
    devtools::load_all() # needs to be called once with updated function
    pkg_root <- dirname(system.file("DESCRIPTION", package = "FastRet"))
    opts <- options()
    on.exit(options(opts), add = TRUE) # reset all options set within `pkg_root/app.R`
    oldplan <- future::plan("sequential") # we can set a multicore future inside `pkg_root/app.R`, this way the cores will be forked again after each reload, which is better for development
    shiny::with_devmode(TRUE, shiny::runApp(pkg_root), verbose = TRUE)
    on.exit(future::plan(oldplan), add = TRUE)
}

#' @title The FastRet GUI
#' @description This function creates the FastRet GUI
#' @param port The port the application should listen on
#' @param host The address the application should listen on
#' @param reload Whether to reload the application when the source code changes
#' @param nsw The number of subworkers each worker is allowed to start. The higher this number, the faster individual tasks like model fitting can be processed.
#' @return A shiny app. This function returns a shiny app that can be run to interact with the model.
#' @keywords public
#' @export
fastret_app <- function(port = 8080,
                        host = "0.0.0.0",
                        reload = FALSE,
                        nsw = 0) {
    shinyApp(
        ui = function(req) fastret_ui(req),
        server = function(input, output, session) fastret_server(input, output, session, nsw),
        options = list(port = port, host = host, quiet = TRUE, launch.browser = FALSE, reload = reload),
        onStart = function() catf("Listening on http://localhost:%s", port)
    )
}

check_cdk_version <- function() {
    if (rcdk::cdk.version() != "2.9") {
        msg <- paste(
            sep = "\n",
            "FastRet requires CDK Version 2.9, but the installed version is %s.",
            "For details about the installation process see https://github.com/CDK-R/rcdklibs."
        )
        msgf <- sprintf(msg, rcdk::cdk.version())
        stop(msgf, call. = FALSE)
    }
}
