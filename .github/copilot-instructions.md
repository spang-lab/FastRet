# GitHub Copilot Instructions for FastRet

> **Important**: This file should be kept up to date whenever any update to the code is made.

## Project Overview

FastRet is an R package for predicting retention times in liquid chromatography. It provides a comprehensive framework for training custom models for specific chromatography columns, predicting retention times using existing models, and adjusting existing models to account for altered experimental conditions.

### Key Features

1. **Train new predictive models** specific for chromatography columns
2. **Use pre-trained models** to predict retention times of molecules
3. **Adjust pre-trained models** to accommodate modifications in chromatography columns
4. **Selective measuring** for optimal subset selection of molecules for re-measurement
5. **Graphical user interface (GUI)** for easy interaction
6. **Command-line interface (CLI)** for programmatic usage

### Technical Details

- **Primary Language**: R (≥ 4.1.0)
- **Dependencies**: Uses chemical descriptors (rCDK), machine learning (glmnet, xgboost), clustering (cluster), and web interface (shiny)
- **Models Supported**: LASSO, Ridge Regression, XGBoost (GBTree)
- **Chemical Descriptors**: 241 features calculated using rCDK version 2.9
- **GUI Framework**: Shiny with background processing using Extended Tasks
- **Caching**: RAM and disk caching for chemical descriptors

## Folder Structure

### Core R Package Structure

- **`R/`** - Main source code directory
  - `app.R` - GUI application setup and configuration
  - `data.R` - Data loading functions and lazy-loaded datasets
  - `getcds.R` - Chemical descriptor calculation with caching
  - `train.R` - Model training, adjustment, and prediction functions
  - `plot.R` - Plotting and visualization functions
  - `prepro.R` - Data preprocessing and feature engineering
  - `sm.R` - Selective measuring (clustering-based subset selection)
  - `server.R` - Shiny server logic
  - `ui.R` - Shiny user interface components
  - `util.R` - Utility functions and imports

### Data and Assets

- **`data/`** - Package datasets (RP.rda)
- **`inst/`** - Installed files
  - `extdata/` - Example datasets (RP.xlsx, RP_adj.xlsx, models)
  - `cachedata/` - Pre-computed chemical descriptors (CDs.rds)

### Documentation and Testing

- **`man/`** - Generated R documentation (.Rd files)
- **`tests/testthat/`** - Test suite covering all major functions
- **`vignettes/`** - User documentation
  - `CLI-Usage.Rmd` - Command-line interface guide
  - `GUI-Usage.Rmd` - Graphical interface guide
  - `Installation.Rmd` - Installation instructions
  - `Contributing.Rmd` - Development guidelines
  - `Package-Internals.Rmd` - Technical documentation

### Development and Deployment

- **`misc/`** - Development utilities
  - `scripts/` - Build and configuration scripts
  - `datasets/` - Dataset management
  - `models/` - Model artifacts
  - `pdfs/` - Generated plots and documentation

## Key Functionality

### Core Workflow Functions

1. **`preprocess_data()`** - Calculates chemical descriptors, removes NAs and near-zero variance features, adds polynomial/interaction terms
2. **`train_frm()`** - Trains LASSO, Ridge, or XGBoost models with cross-validation
3. **`predict.frm()`** - Predicts retention times for new molecules
4. **`adjust_frm()`** - Adjusts existing models for new chromatography conditions
5. **`selective_measuring()`** - Uses PAM clustering to select optimal subset for re-measurement

### Chemical Descriptors System

- **`getCDs()`** - Main function for chemical descriptor calculation with parallel processing
- **`getCDsFor1Molecule()`** - Single molecule descriptor calculation with caching
- **`ram_cache`** - Environment for RAM-based caching of descriptors
- **`CDNames`** - List of chemical descriptor classes used
- **`CDFeatures`** - 241 feature names corresponding to calculated descriptors

### GUI Components

- **`start_gui()`** - Launches the Shiny application
- **`fastret_app()`** - Creates the Shiny app object
- Four main modes: Train Model, Predict RT, Selective Measuring, Adjust Model
- Background processing using futures and Extended Tasks

### Data Management

- **`RP`** - Primary dataset (442 metabolites with RT, SMILES, NAME)
- **`read_rp_xlsx()`**, **`read_rpadj_xlsx()`** - Data loading functions
- **`read_retip_hilic_data()`** - External dataset integration
- Disk and RAM caching system for chemical descriptors

## Coding Guidelines

### Roxygen2 Documentation

- **Always start with tags**, even for title and description: `#' @title ...` instead of `#' ...`
- Every function must be fully documented with `@param`, `@return`, `@examples`, `@keywords`
- Use `@export` for public functions, `@noRd` for internal functions

### Code Style

- **Never use the pipe operator** (`|>` or `%>%`)
- Functions should stay **below 25 lines** when possible
- **Short variable names** are preferred (e.g., `x`, `y` for vectors, `X`, `Y` for 2D objects)
- **Line limit of 80 characters** recommended
- **Maximum indentation of 2 levels** - avoid deep nesting
- Create helper variables instead of deep nesting
- Short, readable function names like `sum`, `sort`, etc. can be nested

### Function Design

- **Avoid deep nesting** of functions or if-else structures
- Prefer creation of helper variables
- Short, readable function names like `sum`, `sort` can be nested
- Use helper functions to break down complex operations
- Validate inputs at function start

### Messaging and Logging

- **Use `catf` instead of `message`, `cat`, or `print`** for info messages during execution
- `catf` supports sprintf-style formatting and timestamps
- Use `verbose` parameters to control output levels (0=silent, 1=progress, 2=detailed)

### Parallel Processing

- Support `nw` (number of workers) parameters for parallelizable operations
- Handle Windows limitations (set `nw=1` automatically)
- Use `mclapply` for simple parallelization, `makeCluster` for complex cases

### Error Handling

- Validate function inputs early with descriptive error messages
- Use `stop()` with clear, actionable error messages
- Handle edge cases (empty data, missing columns, etc.)

### Markdown and RMarkdown Style

- **Empty line after headings** - Always include an empty line after section headings
- **Empty lines around lists** - Add empty lines before and after bullet point and numbered lists
- **Consistent list formatting** - Use `-` for unordered lists, maintain proper indentation
- **Code block formatting** - Use triple backticks with language specification when applicable
- **Table formatting** - Use proper table syntax with header separators

### Testing and Development

- Comprehensive test coverage using `testthat`
- Use `pkg_file()` for package-relative file paths
- Development mode support with `start_gui_in_devmode()`

#### Running Tests

Tests can be executed in several ways:

1. **Run all tests**: `devtools::test()`
2. **Run a specific test file**: `testthat::test_file("tests/testthat/test-getCDs.R")`
3. **Run a specific test file with full path**: `testthat::test_file("c:/Users/tobi/Repos/github/spang/fastret/FastRet/tests/testthat/test-getCDs.R")`

**Important**: Always run `devtools::load_all()` before executing tests to ensure the latest code changes are loaded.

**Critical for Parallel Tests**: Before running tests that involve parallel functionality (multiple workers), the package must be properly installed using `devtools::install()`. This is because spawned subprocesses use `library()` to access package functions and will load the installed version rather than the development version loaded with `devtools::load_all()`. Failure to install will result in "argument 'FUN' is missing" errors when parallel workers cannot find functions that exist in the dev version but not in the installed version.

**Logging Test Output**: To capture test results for analysis, use:
```r
# In R terminal
devtools::load_all()
devtools::test() |> writeLines("test.log")

# Or using system commands (Windows PowerShell)
Rscript -e "devtools::load_all(); devtools::test()" | Tee-Object -FilePath test.log

# Or using system commands (Unix/Linux)
Rscript -e "devtools::load_all(); devtools::test()" | tee test.log
```

#### Test Writing Guidelines

- **Run tests after every code change** - Tests should always be executed after making changes to the codebase
- **Every new function needs tests** - All new functions must have corresponding test cases
- **Speed requirement**: Tests must be as simple as possible and finish within a few milliseconds
- **Performance comments**: If a test cannot be made faster, add an explanation comment explaining why
- **Use subsets**: Use small data subsets (e.g., `RP[1:5, ]` instead of full dataset) for speed
- **Test structure**: Follow the standard `testthat` pattern:
  ```r
  test_that("function_name works correctly", {
      # Setup
      # Test execution  
      # Assertions with expect_*()
  })
  ```

#### Test Performance Patterns

- Cache removal for testing: Use helper functions like `remove_from_cache()` to test both cached and uncached scenarios
- Timing assertions: Include `expect_true(runtime < threshold)` to ensure performance requirements
- Skip expensive tests on CRAN: Use `skip_on_cran()` for tests that require significant resources

#### Common Test Issues and Solutions

**Function Not Found Errors** (`argument "FUN" is missing, with no default`):
- **Problem**: This occurs when `lapply()` or similar functions cannot find the function being called
- **Common causes**: Function not exported, missing imports, or function moved to different file
- **Solution**: Ensure all required functions are properly exported with `@export` or available in namespace
- **Example**: The `getCDsFor1Molecule` function must be exported for use in parallel processing

**Row Names Errors** (`invalid 'row.names' length`):
- **Problem**: Attempting to set row names on a data structure with mismatched dimensions
- **Common causes**: Cache structure changes, data frame dimension mismatches
- **Solution**: Ensure row names match the actual number of rows in the data structure

**Cache Structure Issues**:
- **Problem**: Tests failing after cache system refactoring (matrix vs. list structure)
- **Solution**: Update test helper functions like `remove_from_cache()` to match new cache structure
- **Check**: Verify `ram_cache$CDs` structure matches expectations in test setup

**Missing Dependencies in Tests**:
- **Problem**: Tests fail due to missing variables or functions
- **Solution**: Ensure all dependencies are loaded with `devtools::load_all()` before running tests
- **Check**: Variables like `read_disk`, `write_disk` must be properly defined or defaulted

## Development Patterns

### Model Training Pattern

```r
# 1. Validate input data
validate_inputdata(df, require = c("RT", "SMILES", "NAME"))

# 2. Preprocess (add chemical descriptors)
df <- preprocess_data(df, verbose = verbose, nw = nw)

# 3. Set up cross-validation
folds <- caret::createFolds(y = df$RT, k = nfolds)

# 4. Train models in parallel
models <- mclapply(seq_along(folds), mc.cores = nw, function(i) {
    train_idx <- unname(unlist(folds[-i]))
    test_idx <- folds[[i]]
    fit_model(df[train_idx, ])
})

# 5. Create final model on full data
model <- fit_model(df)

# 6. Return structured frm object
frm <- list(model = model, df = df, cv = cv)
class(frm) <- "frm"
```

### Chemical Descriptor Pattern

```r
# 1. Check RAM cache first
if (smi %in% rownames(ram_cache$CDs)) {
    return(ram_cache$CDs[smi, ])
}

# 2. Calculate descriptors using rCDK
obj <- rcdk::parse.smiles(smi)[[1]]
rcdk::convert.implicit.to.explicit(obj)
rcdk::generate.2d.coordinates(obj)
cds <- rcdk::eval.desc(obj, CDNames, verbose = FALSE)

# 3. Cache results
ram_cache$CDs[smi, ] <- cds
```

### GUI Development Pattern

- Use Extended Tasks for background processing
- Support multiple workers for parallel operations
- Implement progress indicators and console output
- Handle user input validation and error display

## Key Classes and Objects

### S3 Classes

- **`frm`** - FastRet Model object containing model, training data, and CV results
  - Methods: `predict.frm()`, plotting functions

### Important Environments

- **`ram_cache`** - Chemical descriptor caching environment
  - `$CDs` - Data frame of cached descriptors (SMILES as rownames)
  - `$CDRowNr` - List mapping SMILES to row numbers

### Constants

- **`CDNames`** - Chemical descriptor class names (44 descriptors)
- **`CDFeatures`** - Feature names (241 features total)
- **`startModes`** - GUI mode options
- **`strategies`** - Parallel processing strategies

This package represents a mature, well-tested system for retention time prediction with both programmatic and interactive interfaces. When contributing, maintain the existing patterns for caching, parallel processing, error handling, and documentation.
