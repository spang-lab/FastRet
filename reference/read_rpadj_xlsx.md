# Hypothetical retention times

Subset of the data from
[`read_rp_xlsx()`](https://spang-lab.github.io/FastRet/reference/read_rp_xlsx.md)
with some slight modifications to simulate changes in temperature and/or
flowrate.

## Usage

``` r
read_rpadj_xlsx()
```

## Value

A dataframe with 25 rows (metabolites) and 3 columns: RT, SMILES and
NAME.

## Examples

``` r
# \donttest{
x <- read_rpadj_xlsx()
# }
```
