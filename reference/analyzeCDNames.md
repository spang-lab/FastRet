# Analyze Chemical Descriptors Names

Analyze the chemical descriptor names and return a dataframe with their
names and a boolean column indicating if all values are NA.

## Usage

``` r
analyzeCDNames(df, descriptors = rcdk::get.desc.names(type = "all"))
```

## Arguments

- df:

  dataframe with two mandatory columns: "NAME" and "SMILES"

- descriptors:

  Vector of chemical descriptor names

## Value

A dataframe with two columns `descriptor` and `all_na`. Column
`descriptor` contains the names of the chemical descriptors. Column
`all_na` contains a boolean value indicating if all values obtained for
the corresponding descriptor are NA.

## Details

This function is used to analyze the chemical descriptor names and to
identify which descriptors produce only NAs in the test datasets. The
function is used to generate the CDNames object.

## Examples

``` r
X <- analyzeCDNames(df = head(RP, 2), descriptors = CDNames[1:2])
#> 2025-12-17 23:27:07.04 Descriptor 1/2
#> 2025-12-17 23:27:07.87 Descriptor 2/2
```
