# now

Returns the current system time formatted according to the provided
format string.

## Usage

``` r
now(format = "%Y-%m-%d %H:%M:%OS2")
```

## Arguments

- format:

  A string representing the desired time format. Default is "%Y-%m-%d
  %H:%M:%OS2".

## Value

A string representing the current system time in the specified format.

## Examples

``` r
now()            # e.g. "2024-06-12 16:09:32.41"
#> [1] "2025-12-17 23:27:09.51"
now("%H:%M:%S")  # e.g. "16:09:32"
#> [1] "23:27:09"
```
