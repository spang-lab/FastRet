# Package index

## Public API

FastRet’s Public User Interface

- [`adjust_frm()`](https://spang-lab.github.io/FastRet/reference/adjust_frm.md)
  : Adjust an existing FastRet model for use with a new column
- [`clip_predictions()`](https://spang-lab.github.io/FastRet/reference/clip_predictions.md)
  : Clip predictions to observed range
- [`fastret_app()`](https://spang-lab.github.io/FastRet/reference/fastret_app.md)
  : The FastRet GUI
- [`getCDs()`](https://spang-lab.github.io/FastRet/reference/getCDs.md)
  : Get Chemical Descriptors for a list of molecules
- [`plot_frm()`](https://spang-lab.github.io/FastRet/reference/plot_frm.md)
  : Plot predictions for a FastRet model
- [`predict(`*`<frm>`*`)`](https://spang-lab.github.io/FastRet/reference/predict.frm.md)
  : Predict retention times using a FastRet Model
- [`preprocess_data()`](https://spang-lab.github.io/FastRet/reference/preprocess_data.md)
  : Preprocess data
- [`selective_measuring()`](https://spang-lab.github.io/FastRet/reference/selective_measuring.md)
  : Selective Measuring
- [`start_gui()`](https://spang-lab.github.io/FastRet/reference/start_gui.md)
  : Start the FastRet GUI
- [`train_frm()`](https://spang-lab.github.io/FastRet/reference/train_frm.md)
  : Train a new FastRet model (FRM) for retention time prediction

## Internal

Internal helper functions. Although exported, they should not be used
directly. They are only exported to speed up cluster creation during
parallel calculations, as this would require explicit function exporting
of every required private function. Additionally, exporting these
functions allows them to be referenced in the documentation of public
API functions.

- [`CDFeatures`](https://spang-lab.github.io/FastRet/reference/CDFeatures.md)
  : Chemical Descriptor Features
- [`CDNames`](https://spang-lab.github.io/FastRet/reference/CDNames.md)
  : Chemical Descriptors Names
- [`analyzeCDNames()`](https://spang-lab.github.io/FastRet/reference/analyzeCDNames.md)
  : Analyze Chemical Descriptors Names
- [`as_canonical()`](https://spang-lab.github.io/FastRet/reference/as_canonical.md)
  : Canonicalize SMILES
- [`catf()`](https://spang-lab.github.io/FastRet/reference/catf.md) :
  catf function
- [`collect()`](https://spang-lab.github.io/FastRet/reference/collect.md)
  : Collect elements from a list of lists
- [`get_predictors()`](https://spang-lab.github.io/FastRet/reference/get_predictors.md)
  : Extract predictor names from an 'frm' object
- [`init_log_dir()`](https://spang-lab.github.io/FastRet/reference/init_log_dir.md)
  : Initialize log directory
- [`now()`](https://spang-lab.github.io/FastRet/reference/now.md) : now
- [`pkg_file()`](https://spang-lab.github.io/FastRet/reference/pkg_file.md)
  : Get package file
- [`withLineEnd()`](https://spang-lab.github.io/FastRet/reference/withLineEnd.md)
  : Add line end
- [`withSink()`](https://spang-lab.github.io/FastRet/reference/withSink.md)
  : Execute an expression while redirecting output to a file
- [`withStopMessage()`](https://spang-lab.github.io/FastRet/reference/withStopMessage.md)
  : Try expression with predefined error message
- [`withTimeout()`](https://spang-lab.github.io/FastRet/reference/withTimeout.md)
  : Execute an expression with a timeout

## Datasets

Datasets used in FastRet’s examples and vignettes.

- [`RP`](https://spang-lab.github.io/FastRet/reference/RP.md) :
  Retention Times (RT) measured on a Reverse Phase (RP) Column
- [`read_retip_hilic_data()`](https://spang-lab.github.io/FastRet/reference/read_retip_hilic_data.md)
  : Read the HILIC dataset from the Retip package
- [`read_rp_lasso_model_rds()`](https://spang-lab.github.io/FastRet/reference/read_rp_lasso_model_rds.md)
  : LASSO Model trained on RP dataset
- [`read_rp_xlsx()`](https://spang-lab.github.io/FastRet/reference/read_rp_xlsx.md)
  : Read retention times (RT) measured on a reverse phase (RP) column
- [`read_rpadj_xlsx()`](https://spang-lab.github.io/FastRet/reference/read_rpadj_xlsx.md)
  : Hypothetical retention times
