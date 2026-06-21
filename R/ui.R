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
                "<h2>Training data upload</h1>",
                "<p>Here you can upload your own data as Excel file. In order for this to work your file must follow a strict format guideline. If any required columns are missing, the FastRet won&#39;t work correctly. FastRet will always load in the first worksheet of an Excel file. Therefore it is suggested that you reduce your file to one sheet beforehand to avoid any errors.</p>",
                "<h3>Required columns</h2>",
                "<p>The file must consist of the following columns (case sensitive):</p>",
                "<ul>",
                "<li>RT: Retention time of your molecules. Can be any numeric input, minutes or seconds. Remember what you put in when you analyze the predictions, since those will be on the same scale as your input data.</li>",
                "<li>NAME: you can put in any characters you like. Preferably the names of your molecules.</li>",
                "<li>SMILES: Isomeric or canonical SMILES. This information is used to calculate chemical descriptors with the chemistry development kit</li>",
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
                "<h2>Model upload</h1>",
                "<p>",
                "Here you need to upload a prediction model generated with this program in the <em>Train new Model</em> mode.",
                "This Model can also be read in and used from within R by calling",
                "<pre style='white-space: pre-wrap;'>",
                "model <- readRDS('path/to/model.rds')",
                "coef(model$model)  # show model coeffcients (LASSO only)",
                "model$df      # show the predictor set",
                "</pre>",
                "For further details see the FastRet online documentation.",
                "</p>"
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
                "<h2>Method Selection</h1>",
                "<p>Here you can choose by which method the regression model should be trained on. You can choose between Lasso or XGBoost. </p>",
                "<h3>Lasso</h2>",
                "<p>Lasso (Least absolut shrinkage and selection operator) is based on the Least Minimum Square approach with the extension of a L1 penalty norm. This leads to a selection of variables as well as a generalization of the trained model.<br>Lasso was implemented with the R-package glmnet [2].</p>",
                "<h3>XGBoost</h2>",
                "<p>XGBoost is a more soffisticated Machine Learning method based on Boosted Regression Trees (BRT) [3]. The main difference to random forest is, that trees are not trained independant from each other but each tree is built with a loss function based on its predecessor. It was implemented with the R-package XGBoost [4].</p>",
                "<h3>References</h2>",
                "<p>",
                "[1] Santosa, Fadil; Symes, William W. (1986). &quot;Linear inversion of band-limited reflection seismograms&quot;. <em>SIAM Journal on Scientific and Statistical Computing</em>. SIAM. <strong>7</strong> (4): 1307<e2><80><93>1330",
                "[2] Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).",
                "Regularization Paths for Generalized Linear Models via",
                "Coordinate Descent. Journal of Statistical Software, 33(1),",
                "1-22.",
                "[3] Jerome H. Friedman. &quot;Greedy function approximation: A gradient boosting machine..&quot; Ann. Statist. 29 (5) 1189 - 1232, October 2001",
                "[4] Tianqi Chen et. Al, (2021). xgboost: Extreme Gradient Boosting. R package",
                "version 1.4.1.1.",
                "</p>"
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
                "<h2>Cluster Calculation</h1>",
                "<p>Here you can choose how many clusters should be calculated. The programm will calculate the best k molecules to be measured for a retention time prediction. It uses a combination of Ridge Regression and PAM (k-medoids) clustering to determine the best representatives of your dataset. Representatives as well as their corresponding clusters can be downloaded afterwards as an excel file. This step should be used once you have a predictive model and/or data set and want to adjust it for a new chromatography column with different gradient/temperature etc.</p>"
            )
        ),
        with_helptext(
            shiny::selectInput(
                inputId = "siSmVariant",
                label = "Variant",
                choices = list(
                    "RT weighted like top descriptor (recommended)" = "max_ridge_coef",
                    "RT left unscaled" = "1",
                    "Exclude RT" = "0",
                    "RT only" = "inf"
                ),
                selected = "max_ridge_coef"
            ),
            content = paste(
                "<h2>Variant</h2>",
                "<p>Controls how strongly the retention time of the already measured compounds influences the selection. The recommended option weights it like the most important chemical descriptor. The remaining options let you exclude the retention time, leave it unscaled, or rely on it alone.</p>"
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
                "<h2 id='single-prediction'>Single Prediction</h1>",
                "<p>Here you can input a single SMILES string to predict the retention time of a molecule. The SMILES string has to be in the same format as the SMILES strings in the training data set. The prediction will be done with the model you uploaded.</p>"
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
                "<h2 id='prediction-data-upload'>Prediction Data Upload</h1>",
                "<p>This file input has to be an excel file with columns NAME and SMILES</p>"
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
                label = "Data for prediction adustment as xlsx file",
                accept = ".xlsx"
            ),
            content = paste(
                "<h2>Adjustment data</h2>",
                "<p>This file input has to be an excel file with columns RT, NAME and SMILES</p>"
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
