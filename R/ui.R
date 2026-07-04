# Main (Private) #####

fastret_ui <- function(req = NULL) {
    shiny::tagList(
        shiny::tags$head(
            shiny::tags$style(shiny::HTML(
                "#shiny-notification-panel {
                    top: 0;
                    right: 0;
                    bottom: unset;
                    left: unset;
                    margin-left: auto;
                    margin-right: auto;
                    width: auto;
                    max-width: 450px;
                }
                code {
                    color: black;
                    background-color: #f5f5f5;
                    border: 1px solid #f7f7f7;
                }
                "
            ))
        ),
        shiny::navbarPage(
            title = "FastRet",
            id = "navbar",
            selected = getOption("FastRet.UI.startMode", "Train new Model"),
            ui_tab_train(),
            ui_tab_select(),
            ui_tab_adjust(),
            ui_tab_predict(),
            ui_privacy_policy(),
            ui_contact(),
            ui_about(),
        ),
        # shinythemes::themeSelector(),
        shinyjs::useShinyjs()
    )
}

# Mode tabs (Private) #####

# Each of the four modes is its own navbar tab. The visible label is short; the
# full description is shown as a native browser tooltip on hover (the `title`
# attribute on the label span). The `value` keeps the original mode strings so
# the server-side gates (see `ui_*_results`) only needed `siMode` -> `navbar`.

mode_tab_title <- function(short, long) {
    shiny::span(short, title = long)
}

ui_tab_train <- function() {
    shiny::tabPanel(
        title = mode_tab_title("Train", "Train new Model"),
        value = "Train new Model",
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                ui_sidebar_heading("Train new Model", paste(
                    "<h2>Train new Model</h2>",
                    "<p>Build a retention-time model from your own measurements: the names, SMILES and retention times of compounds measured on your chromatography column. The trained model predicts retention times for new molecules, and can later be adjusted to other columns.</p>"
                )),
                ui_data_upload("ubTrainXlsx", "toTrainXlsxError"),
                ui_train_controls()
            ),
            shiny::mainPanel(
                shiny::uiOutput("ui_train_results")
            )
        )
    )
}

ui_tab_select <- function() {
    shiny::tabPanel(
        title = mode_tab_title("Select", "Selective Measuring"),
        value = "Selective Measuring",
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                ui_sidebar_heading("Selective Measuring", paste(
                    "<h2>Selective Measuring</h2>",
                    "<p>Pick the most informative compounds to measure on a new column. FastRet returns the most representative molecules of your dataset using ridge regression and k-medoids clustering, so you can adjust a model to the new column with as little lab work as possible.</p>"
                )),
                ui_data_upload("ubSmXlsx", "toSmXlsxError"),
                ui_sm_controls()
            ),
            shiny::mainPanel(
                shiny::uiOutput("ui_sm_results")
            )
        )
    )
}

ui_tab_adjust <- function() {
    shiny::tabPanel(
        title = mode_tab_title("Adjust", "Adjust existing Model"),
        value = "Adjust existing Model",
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                ui_sidebar_heading("Adjust existing Model", paste(
                    "<h2>Adjust existing Model</h2>",
                    "<p>Adapt an existing model to a changed setup, such as a different temperature, pH or column age. Upload the model together with a few compounds re-measured on the new column, and FastRet fits a small correction on top of the model.</p>"
                )),
                ui_model_upload("ubAdjFRM", "toAdjFRMError"),
                ui_adjust_controls()
            ),
            shiny::mainPanel(
                shiny::uiOutput("ui_adjust_results")
            )
        )
    )
}

ui_tab_predict <- function() {
    shiny::tabPanel(
        title = mode_tab_title("Predict", "Predict Retention Times"),
        value = "Predict Retention Times",
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                ui_sidebar_heading("Predict Retention Times", paste(
                    "<h2>Predict Retention Times</h2>",
                    "<p>Predict retention times for new molecules using a previously trained or adjusted model. Enter a single SMILES or upload many at once as an Excel file.</p>"
                )),
                ui_model_upload("ubPredFRM", "toPredFRMError"),
                ui_predict_controls()
            ),
            shiny::mainPanel(
                shiny::uiOutput("ui_predict_results")
            )
        )
    )
}

# Info pages (Private) #####

ui_privacy_policy <- function() {
    shiny::tabPanel(
        title = "Privacy Policy",
        shiny::fluidPage(
            shiny::HTML(
                "<div class='mainpanel'>",
                "    <h3>Cookies</h3>",
                "    <div>This website does not use cookies.</div>",
                "    <h3> Server Log</h3>",
                "    <div>The web server keeps a log of all requests, with the following data:</div>",
                "    <ul>",
                "        <li>The request IP adress</li>",
                "        <li>Date and Time of the request</li>",
                "        <li>request type and path</li>",
                "        <li>the User-Agent of the web browser</li>",
                "    </ul>",
                "    <div>This data is only used to diagnose technical problems.</div>",
                "    <h3>Web Analytics / Other Tracking</h3>",
                "    <div>There are no other tracking methods.</div>",
                "    <h3>Privacy Contact</h3>",
                "    <a href='http://www.uni-regensburg.de/universitaet/datenschutzbeauftragte/index.html'>",
                "        Datenschutzbeauftrage der Universit&auml;t",
                "    </a>",
                "</div>"
            )
        )
    )
}

ui_contact <- function() {
    shiny::tabPanel(
        title = "Contact",
        shiny::fluidPage(
            shiny::HTML(
                "<div class='mainpanel'>",
                "    <h3>Contact</h3>",
                "    <div>",
                "        <address>",
                "            Dr. Katja Dettmer-Wilde<br/>",
                "            Institute of Functional Genomics<br/>",
                "            University of Regensburg<br/>",
                "            Am BioPark 9<br/>",
                "            93053 Regensburg, Germany<br/>",
                "            <abbr title='Phone'>P: </abbr>+49 941 943 5051<br/>",
                "            <abbr title='Email'>M: </abbr>katja.dettmer@klinik.uni-regensburg.de",
                "        </address>",
                "    </div>",
                "</div>",
            )
        )
    )
}

ui_about <- function() {
    shiny::tabPanel(
        title = "About",
        shiny::fluidPage(
            shiny::div(
                shiny::h3("Purpose"),
                shiny::tags$p("FastRet is an R package for predicting retention times in liquid chromatography. It can be used through the R console or through a graphical user interface (GUI)."),
                shiny::tags$p("The package's key features include the ability to:"),
                shiny::tags$ul(
                    shiny::tags$li("Train new predictive models specific for your own chromatography column"),
                    shiny::tags$li("Use pre-trained models to predict retention times of molecules"),
                    shiny::tags$li("Adjust pre-trained models to accommodate modifications in chromatography columns")
                ),
                shiny::tags$p("For further details see FastRets ", shiny::tags$a(href = "https://spang-lab.github.io/FastRet/", "documentation site"))
            ),
            shiny::div(
                shiny::h3("Version Info"),
                shiny::pre(paste(capture.output(sessionInfo()), collapse = "\n"))
            )
        )
    )
}

# Sidebar controls (Private) #####

# `inputId`/`errorId` are parameterised because the same upload widget appears in
# two tabs each (data: Train + Select; model: Predict + Adjust). navbarPage
# renders every tab at once, so each instance needs a unique HTML id.
ui_data_upload <- function(inputId, errorId) {
    shiny::tagList(
        with_helptext(
            shiny::fileInput(
                inputId = inputId,
                label = "Data as xlsx file",
                accept = ".xlsx"
            ),
            content = paste(
                "<h2>Training data upload</h2>",
                "<p>Upload your own measurements as an Excel file. FastRet searches every worksheet for the first one that contains the required columns, and any extra columns are ignored.</p>",
                "<h3>Required columns (case sensitive)</h3>",
                "<ul>",
                "<li><code>RT</code>: retention time of each molecule, in any numeric unit (e.g. minutes or seconds). Predictions are reported on the same scale.</li>",
                "<li><code>NAME</code>: a label for each molecule, for example its name.</li>",
                "<li><code>SMILES</code>: isomeric or canonical SMILES, used to compute the chemical descriptors via the Chemistry Development Kit.</li>",
                "</ul>"
            )
        ),
        shiny::div(shiny::textOutput(errorId), style = "color: red;")
    )
}

ui_model_upload <- function(inputId, errorId) {
    shiny::tagList(
        with_helptext(
            shiny::fileInput(
                inputId = inputId,
                label = "Upload a pretrained Model",
                accept = c(".rds", ".RDS")
            ),
            content = paste(
                sep = "\n",
                "<h2>Model upload</h2>",
                "<p>Upload a prediction model that you trained in the <em>Train new Model</em> mode. You can also load and inspect it in R:</p>",
                "<pre><code>model <- readRDS(\"path/to/model.rds\")",
                "coef(model$model)  # model coefficients (Lasso only)",
                "model$df           # the predictor set</code></pre>",
                "<p>See the FastRet online documentation for details.</p>"
            )
        ),
        shiny::div(shiny::textOutput(errorId), style = "color: red;")
    )
}

ui_train_controls <- function() {
    shiny::tagList(
        with_helptext(
            shiny::radioButtons(
                inputId = "rbMethod",
                label = "Method",
                choices = list("XGBoost (recommended)" = 2, "Lasso" = 1),
                selected = 2
            ),
            content = paste(
                "<h2>Method selection</h2>",
                "<p>Choose how the regression model is trained: Lasso or XGBoost.</p>",
                "<h3>Lasso</h3>",
                "<p>Lasso (least absolute shrinkage and selection operator) extends least-squares regression with an L1 penalty. This selects a subset of predictors and helps the model generalise. Implemented with the R package <code>glmnet</code> [2].</p>",
                "<h3>XGBoost</h3>",
                "<p>XGBoost is a boosted-regression-tree method [3]. Unlike a random forest, each tree is fitted on the residuals of its predecessors. Implemented with the R package <code>xgboost</code> [4].</p>",
                "<h3>References</h3>",
                "<p>",
                "[1] Santosa F, Symes WW (1986). Linear inversion of band-limited reflection seismograms. <em>SIAM Journal on Scientific and Statistical Computing</em>, 7(4), 1307-1330.<br>",
                "[2] Friedman J, Hastie T, Tibshirani R (2010). Regularization paths for generalized linear models via coordinate descent. <em>Journal of Statistical Software</em>, 33(1), 1-22.<br>",
                "[3] Friedman JH (2001). Greedy function approximation: a gradient boosting machine. <em>Annals of Statistics</em>, 29(5), 1189-1232.<br>",
                "[4] Chen T, et al. (2021). xgboost: Extreme Gradient Boosting. R package version 1.4.1.1.",
                "</p>"
            )
        ),
        with_helptext(
            shiny::checkboxInput(
                inputId = "cbTuneGrid",
                label = "Tune hyperparameters (grid search)",
                value = FALSE
            ),
            content = paste(
                "<h2>Hyperparameter tuning</h2>",
                "<p>Applies only to XGBoost. When enabled, the model's",
                "hyperparameters are chosen by a small cross-validated grid search",
                "(8 combinations) instead of the fast default settings. This can",
                "improve accuracy at the cost of a few extra seconds of training",
                "time. Has no effect for Lasso.</p>"
            )
        ),
        with_helptext(
            shiny::numericInput(
                inputId = "seed",
                label = "Seed",
                value = 42
            ),
            content = paste(
                "<h2>Seed</h2>",
                "<p>The seed for the random number generator used during model training and cross-validation. It is fixed to a constant default so that training the same data twice yields identical models. Change it only if you want to assess the variability of the results across different random seeds.</p>"
            )
        ),
        bslib::input_task_button(
            id = "btnTrain",
            label = "Train Model",
            type = "default",
            icon = shiny::icon("play")
        ),
        shiny::downloadButton(
            outputId = "dbSaveModel",
            label = "Save Model",
            style = "display: none;"
        ),
        shiny::downloadButton(
            outputId = "dbSavePredictorSet",
            label = "Save Predictor Set",
            style = "display: none;"
        ),
        consoleOutput(
            divID = "divTrainLogs",
            vtoID = "vtoTrainLogs"
        )
    )
}

ui_sm_controls <- function() {
    shiny::tagList(
        with_helptext(
            shiny::numericInput(
                inputId = "niK", # niK for "numeric input K"
                label = "k Cluster",
                value = 25
            ),
            content = paste(
                "<h2>Number of compounds (k)</h2>",
                "<p>How many compounds to select. FastRet combines ridge regression with PAM (k-medoids) clustering and returns the <code>k</code> most representative molecules of your dataset. Measure these on the new column and use them to adjust your model. You can download the representatives and their clusters as an Excel file.</p>"
            )
        ),
        with_helptext(
            shiny::selectInput(
                inputId = "siSmVariant",
                label = "Variant",
                choices = list(
                    "SMmax (recommended)" = "max_ridge_coef",
                    "SM1" = "1",
                    "SM0" = "0",
                    "SMinf" = "inf"
                ),
                selected = "max_ridge_coef"
            ),
            content = paste(
                "<h2>Selective measuring variant</h2>",
                "<p>Controls how strongly the retention time (RT) of the already measured compounds guides the selection. The variants differ only in how RT is scaled before clustering:</p>",
                "<ul>",
                "<li><b>SMmax</b> (recommended): RT is scaled by the largest ridge-regression coefficient, so it carries about the same weight as the most important chemical descriptor.</li>",
                "<li><b>SM1</b>: RT is used as is (standardized, no extra scaling).</li>",
                "<li><b>SM0</b>: RT is excluded; clustering uses the chemical descriptors only.</li>",
                "<li><b>SMinf</b>: chemical descriptors are ignored; clustering uses RT alone.</li>",
                "</ul>"
            )
        ),
        with_helptext(
            shiny::numericInput(
                inputId = "niSmSeed",
                label = "Seed",
                value = 42
            ),
            content = paste(
                "<h2>Seed</h2>",
                "<p>Seed for the random number generator, so the same data and settings always select the same compounds. Change it only to check how sensitive the selection is to the seed.</p>"
            )
        ),
        bslib::input_task_button(
            id = "btnSM",
            label = "Calculate clusters and medoids",
            type = "default"
        ),
        shiny::downloadButton(
            outputId = "dbSaveCluster",
            label = "Download clustering as xlsx",
            style = "display: none;"
        ),
        consoleOutput(
            divID = "divSMLogs",
            vtoID = "vtoSMLogs"
        )
    )
}

ui_predict_controls <- function() {
    shiny::tagList(
        with_helptext(
            shiny::textInput(
                inputId = "tiPredSmiles",
                label = "Input SMILES",
                value = ""
            ),
            content = paste(
                "<h2 id='single-prediction'>Single prediction</h2>",
                "<p>Enter one SMILES string to predict its retention time with the uploaded model. Please use the same SMILES format as in the training data.</p>"
            )
        ),
        shiny::div(shiny::textOutput("toPredSmilesError"), style = "color: red;"),
        with_helptext(
            shiny::fileInput(
                inputId = "ubPredXlsx",
                label = "Upload SMILES as xlsx",
                accept = ".xlsx"
            ),
            content = paste(
                "<h2 id='prediction-data-upload'>Prediction data upload</h2>",
                "<p>Upload an Excel file with the columns <code>NAME</code> and <code>SMILES</code> to predict many molecules at once.</p>"
            )
        ),
        shiny::div(shiny::textOutput("toPredXlsxError"), style = "color: red;"),
        bslib::input_task_button(
            id = "btnPred",
            label = "Predict",
            type = "default"
        ),
        shiny::downloadButton(
            outputId = "dbSavePred",
            label = "Save predictions",
            style = "display: none;"
        ),
        consoleOutput(
            divID = "divPredLogs",
            vtoID = "vtoPredLogs"
        )
    )
}

ui_adjust_controls <- function() {
    shiny::tagList(
        with_helptext(
            shiny::fileInput(
                inputId = "ubAdjXlsx",
                label = "Data for prediction adjustment as xlsx file",
                accept = ".xlsx"
            ),
            content = paste(
                "<h2>Adjustment data</h2>",
                "<p>Upload an Excel file with the columns <code>RT</code>, <code>NAME</code> and <code>SMILES</code> for the compounds re-measured on the new column.</p>"
            )
        ),
        shiny::div(shiny::textOutput("toAdjXlsxError"), style = "color: red;"),
        with_helptext(
            shiny::radioButtons(
                inputId = "rbAdjMethod",
                label = "Adjustment method",
                choices = list("Lasso (recommended)" = "lasso", "Linear model" = "lm", "XGBoost" = "gbtree"),
                selected = "lasso"
            ),
            content = paste(
                "<h2>Adjustment method</h2>",
                "<p>Model adjustment overlays a small correction model on top of an existing FastRet model so that it predicts retention times for an altered chromatographic setup (e.g. changed temperature, pH, or column age). Choose the model family used for this correction:</p>",
                "<ul>",
                "<li><b>Lasso</b> (recommended) and <b>XGBoost</b> adjust using the base retention time together with the molecular descriptors, so they can capture compound-specific shifts.</li>",
                "<li><b>Linear model</b> fits a simple straight-line correction from the base retention time only (no descriptors).</li>",
                "</ul>"
            )
        ),
        bslib::input_task_button(
            id = "btnAdj",
            label = "Adjust Model",
            type = "default"
        ),
        shiny::downloadButton(
            outputId = "dbSaveAdjModel",
            label = "Save adjusted model",
            style = "display: none;"
        ),
        consoleOutput(
            divID = "divAdjLogs",
            vtoID = "vtoAdjLogs"
        )
    )
}

# Results (Private) #####

ui_train_results <- function(SE) {
    shiny::req(SE$RV$tblTrainResults, SE$input$navbar == "Train new Model")
    htmltools::div(
        id = "ui_train_results",
        shiny::tabsetPanel(
            shiny::tabPanel(
                title = "CV Performance",
                shiny::plotOutput("poTrainPerfCV")
            ),
            shiny::tabPanel(
                title = "Training performance",
                shiny::plotOutput("poTrainPerf")
            )
        ),
        DT::DTOutput("tblTrainResults")
    )
}

ui_sm_results <- function(SE) {
    shiny::req(SE$RV$cluster_calc, SE$input$navbar == "Selective Measuring")
    htmltools::div(
        id = "ui_sm_results",
        htmltools::h3("Medoids"),
        DT::DTOutput("tblMedoids"),
        htmltools::h3("Full clustering"),
        DT::DTOutput("tblClustering")
    )
}

ui_predict_results <- function(SE) {
    shiny::req(SE$input$navbar == "Predict Retention Times")
    catf("Showing prediction results")
    htmltools::div(
        id = "ui_predict_results",
        DT::DTOutput("tblPredResults")
    )
}

ui_adjust_results <- function(SE) {
    shiny::req(SE$input$navbar == "Adjust existing Model")
    htmltools::div(
        id = "ui_adjust_results",
        shiny::plotOutput("poAdjPerfCV"),
        shiny::plotOutput("poAdjPerf")
    )
}

# Helpers (Private) #####

with_helptext <- function(..., content) {
    shinyhelper::helper(
        ...,
        icon = "question-circle",
        colour = "#696969",
        type = "inline",
        content = content
    )
}

# Heading shown at the top of a tab's sidebar, with a help button next to it that
# explains what the mode is for.
ui_sidebar_heading <- function(title, content) {
    with_helptext(
        htmltools::h4(title, style = "display: inline-block; margin: 0 0.4em 0.6em 0;"),
        content = content
    )
}

consoleOutput <- function(divID, vtoID) {
    vto <- shiny::verbatimTextOutput(outputId = vtoID)
    vto$attribs$style <- "
        display: block;
        line-height: 1.5em;
        height: 15em;
        overflow-y: auto;
        border: 1px solid #e3e3e3;
        resize: vertical;
        word-wrap: break-word;
        overflow-wrap: break-word;
        word-break: break-all;
        white-space: shiny::pre-wrap;
    "
    shiny::div(
        id = divID,
        style = "margin-top: 15px; margin-bottom: 15px;",
        htmltools::tags$label("Console Log", class = "control-label"),
        vto
    )
}
