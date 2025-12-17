# Get Chemical Descriptors for a list of molecules

Calculate Chemical Descriptors (CDs) for a list of molecules. Molecules
can appear multiple times in the list.

## Usage

``` r
getCDs(df, verbose = 1, nw = 1, keepdf = TRUE)
```

## Arguments

- df:

  dataframe with two mandatory columns: "NAME" and "SMILES"

- verbose:

  0: no output, 1: progress, 2: more progress and warnings

- nw:

  number of workers for parallel processing

- keepdf:

  If TRUE, `cbind(df, CDs)` is returned. Else `CDs`.

## Value

A dataframe with all input columns (if `keepdf` is TRUE) and chemical
descriptors as remaining columns.

## Examples

``` r
cds <- getCDs(head(RP, 3), verbose = 1, nw = 1)
#> 2025-12-17 23:24:47.77 Calculating chemical descriptors for 3 molecules
#> 2025-12-17 23:24:47.77 Finished calculating chemical descriptors in 0.00s
```
