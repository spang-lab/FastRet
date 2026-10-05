# FastRet

FastRet is an R package for predicting retention times in liquid
chromatography. It can be used through the R console or through a
graphical user interface (GUI). The package is described in [Fadil et
al. (2026)](https://doi.org/10.1021/acs.jcim.6c01344). The package’s key
features include the ability to

1.  Train new predictive models specific for your own chromatography
    column
2.  Use pre-trained models to predict retention times of molecules
3.  Adjust pre-trained models to accommodate modifications in
    chromatography columns
4.  Select a small, representative subset of molecules to re-measure on
    a modified column (Selective Measuring), so that an adjustment model
    can be trained cheaply

## Installation

FastRet requires a Java SDK (for the chemical-descriptor package rcdk).
Once Java is available, you can install the released version of FastRet
from [CRAN](https://cran.r-project.org/package=FastRet) by entering the
following commands in an R session:

``` r

if (Sys.which("java")[1] == "") stop("Please install a Java SDK first.")
install.packages("FastRet")
```

To install the development version from
[GitHub](https://github.com/spang-lab/FastRet) instead, use:

``` r

install.packages("pak")
pak::pkg_install("spang-lab/FastRet")
```

For further details see
[Installation](https://spang-lab.github.io/FastRet/articles/Installation.html).

## Usage

The easiest way to use FastRet is through its GUI. A hosted version of
the GUI is available at <https://fastret.spang-lab.de>, so you can try
FastRet without installing anything. To start the GUI locally, [install
the package](#installation) and then run the following command in an
interactive R terminal:

``` r

FastRet::start_gui()
```

After running the above code, you should see an output like

    Listening on http://localhost:8080

in your R console. This means that the GUI is now running and you can
access it via the URL `http://localhost:8080` in your browser. If your
terminal supports it, you can also just click on the displayed link.

![The FastRet GUI, showing the Train
tab](https://raw.githubusercontent.com/spang-lab/FastRet/main/vignettes/GUI-Usage/train.png)

The GUI is organized into four modes, available as tabs in the
navigation bar: *Train*, *Select*, *Adjust* and *Predict* (the tabs use
short labels; hovering shows the full names *Train new Model*,
*Selective Measuring*, *Adjust existing Model* and *Predict Retention
Times*). By default, the GUI opens on the *Train* tab. For more
information about the individual modes and the various input fields,
click on the little question mark symbols next to the different input
fields or have a look at the documentation page for [GUI
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

## Citation

To cite FastRet in publications, please use:

Fadil F, Schmidt T, Amesoeder C, Heckscher S, Schoen M, Gronwald W,
Oefner PJ, Spang R, Dettmer K (2026). FastRet: Fast and Simple Retention
Time Prediction in Liquid Chromatography. *Journal of Chemical
Information and Modeling*, 66(16), 10412-10425.
[doi:10.1021/acs.jcim.6c01344](https://doi.org/10.1021/acs.jcim.6c01344)

A BibTeX entry is available via `citation("FastRet")` in R.
