# Retention times re-measured under a steeper gradient

25 metabolites from the
[RP](https://spang-lab.github.io/FastRet/reference/RP.md) dataset,
re-measured on the same reverse-phase column under a **steeper
gradient** (real data, dataset `RP_Steep` in the FastRet paper). Used to
demonstrate model adjustment via
[`adjust_frm()`](https://spang-lab.github.io/FastRet/reference/adjust_frm.md).

## Usage

``` r
read_rpadj_xlsx()
```

## Value

A dataframe with 25 rows (metabolites) and 4 columns: RT, NAME, SMILES
and INCHIKEY.

## Examples

``` r
# \donttest{
x <- read_rpadj_xlsx()
# }
```
