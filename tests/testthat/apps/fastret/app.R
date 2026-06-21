# Minimal launcher used by the shinytest2 end-to-end tests
# (tests/testthat/test-gui-e2e.R). It builds the exact same Shiny app that
# `FastRet::fastret_app()` builds, but without the hard-coded port so that
# shinytest2 can pick a free one. A sequential future plan keeps the
# ExtendedTask jobs in this process, so the installed FastRet namespace is used
# directly and no worker (sub)processes need to be spawned.
library(shiny)
library(FastRet)
future::plan("sequential")
shiny::shinyApp(
    ui = FastRet:::fastret_ui,
    server = function(input, output, session) {
        FastRet:::fastret_server(input, output, session, nsw = 1)
    }
)
