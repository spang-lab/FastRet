# Retention Times (RT) measured on a Reverse Phase (RP) Column

Retention time data from a reverse phase liquid chromatography measured
with a temperature of 35\\^\circ\\C and a flowrate of 0.3ml/min. The
same data is available as an xlsx file in the package. To read it into R
use
[`read_rp_xlsx()`](https://spang-lab.github.io/FastRet/reference/read_rp_xlsx.md).
@format A dataframe of 442 metabolites with the following columns:

- RT:

  Retention time

- SMILES:

  SMILES notation of the metabolite

- NAME:

  Name of the metabolite

## Usage

``` r
RP
```

## Format

An object of class `data.frame` with 442 rows and 3 columns.

## Source

Measured by the Institute of Functional Genomics at the University of
Regensburg.

## See also

read_rp_xlsx
