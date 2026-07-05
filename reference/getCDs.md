# Get Chemical Descriptors for a list of molecules

Calculate Chemical Descriptors (CDs) for a list of molecules. Molecules
can appear multiple times in the list.

## Usage

``` r
getCDs(df, verbose = 1, nw = 1, keepdf = TRUE, cache = TRUE)
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

- cache:

  If TRUE (default), descriptors are cached on disk (an SQLite database
  under tools::R_user_dir`("FastRet", "cache")`): already-known
  descriptors are read from the cache and newly computed ones are
  written back, so repeated calls on the same molecules are fast and the
  cache is shared across processes. If FALSE, the cache is bypassed
  entirely and every descriptor is (re)computed via rCDK – useful for
  benchmarking the uncached runtime.

## Value

A dataframe with all input columns (if `keepdf` is TRUE) and chemical
descriptors as remaining columns.

## Examples

``` r
cds <- getCDs(head(RP, 3), verbose = 1, nw = 1)
#> 2026-07-05 17:57:58.94 Obtaining chemical descriptors for 3 molecules
#> 2026-07-05 17:57:58.96 Finished obtaining chemical descriptors in 0.02s
```
