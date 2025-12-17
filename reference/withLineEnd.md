# Add line end

Checks if a string ends with a newline character. If not, a newline
character is appended.

## Usage

``` r
withLineEnd(x)
```

## Arguments

- x:

  A string.

## Value

The input string with a newline character at the end if it was not
already present.

## Examples

``` r
cat(withLineEnd("Hello"))
#> Hello
```
