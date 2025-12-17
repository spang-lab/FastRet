# Collect elements from a list of lists

Takes a list of lists where each inner list has the same names. It
returns a list where each element corresponds to a name of the inner
list that is extracted from each inner list. Especially useful for
collecting results from lapply.

## Usage

``` r
collect(xx)
```

## Arguments

- xx:

  A list of lists where each inner list has the same names.

## Value

A list where each element corresponds to a name of the inner list that
is extracted from each inner list.

## Examples

``` r
xx <- lapply(1:3, function(i) list(a = i, b = i^2, c = i^3))
ret <- collect(xx)
```
