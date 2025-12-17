# Canonicalize SMILES

Convert SMILES to canonical form.

## Usage

``` r
as_canonical(smiles)
```

## Arguments

- smiles:

  Character vector of SMILES.

## Value

Character vector of canonical SMILES.

## Examples

``` r
as_canonical(c("CCO", "C(C)O"))
#> [1] "OCC" "OCC"
```
