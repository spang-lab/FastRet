# Read retention times (RT) measured on a reverse phase (RP) column

Reads retention times from a reverse phase liquid chromatography
experiment, performed at 35\\^\circ\\C and a flow rate of 0.3 mL/min.
The data is also available as a dataframe in the package; to access it
directly, use [RP](https://spang-lab.github.io/FastRet/reference/RP.md).

## Usage

``` r
read_rp_xlsx()
```

## Source

Measured by the Institute of Functional Genomics at the University of
Regensburg.

## Value

A dataframe of 442 metabolites with columns `RT`, `SMILES` and `NAME`.

## See also

RP

## Examples

``` r
x <- read_rp_xlsx()
all.equal(x, RP)
#> [1] TRUE
```
