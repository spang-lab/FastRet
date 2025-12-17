# FastRet

FastRet is an R package for predicting retention times in liquid
chromatography. It can be used through the R console or through a
graphical user interface (GUI). The package’s key features include the
ability to

1.  Train new predictive models specific for your own chromatography
    column
2.  Use pre-trained models to predict retention times of molecules
3.  Adjust pre-trained models to accommodate modifications in
    chromatography columns

## Installation

You can install the development version of FastRet from
[GitHub](https://github.com/) by entering the following commands in an R
session:

``` r
if (Sys.which("java")[1] == "") stop("Please install a Java SDK first.")
install.packages("pak")
pak::pkg_install("spang-lab/FastRet")
```

For further details see
[Installation](https://spang-lab.github.io/FastRet/articles/Installation.html).

## Usage

The easiest way to use FastRet is through its GUI. To start the GUI,
[install the package](#installation) and then run the following command
in an interactive R terminal:

``` r
FastRet::start_gui()
```

After running the above code, you should see an output like

    Listening on http://localhost:8080

in your R console. This means that the GUI is now running and you can
access it via the URL `http://localhost:8080` in your browser. If your
terminal supports it, you can also just click on the displayed link.

![start-page.png](https://raw.githubusercontent.com/spang-lab/FastRet/main/vignettes/GUI-Usage/start-page.png)![mode-help.png](https://raw.githubusercontent.com/spang-lab/FastRet/main/vignettes/GUI-Usage/mode-help.png)

By default, the GUI opens in Mode *Train new Model*. To apply or adjust
pretrained models, select mode *Predict Retention Time* or *Adjust
existing Model* instead. For more information about the individual modes
and the various input fields, click on the little question mark symbols
next to the different input fields or have a look at the documentation
page for [GUI
Usage](https://spang-lab.github.io/FastRet/articles/GUI-Usage.html).

## Documentation

FastRet’s documentation is available at
[spang-lab.github.io/FastRet](https://spang-lab.github.io/FastRet/). It
includes pages about

- [GUI
  Usage](https://spang-lab.github.io/FastRet/articles/GUI-Usage.html)
- [CLI
  Usage](https://spang-lab.github.io/FastRet/articles/CLI-Usage.html)
- [Package
  Internals](https://spang-lab.github.io/FastRet/articles/Package-Internals.html)
- [Contribution
  Guidelines](https://spang-lab.github.io/FastRet/articles/Contributing.html)
- [Function
  Reference](https://spang-lab.github.io/FastRet/reference/index.html)
