# Execute an expression with a timeout

Execute an expression with a timeout

## Usage

``` r
withTimeout(expr, timeout = 2)
```

## Arguments

- expr:

  The expression to execute

- timeout:

  The timeout in seconds. Default is 2.

## Value

The result of the expression

## Examples

``` r
withTimeout(
     cat("This works\n"),
     timeout = 0.2
)
#> This works
try(silent = TRUE, withTimeout(
    expr = {Sys.sleep(0.2); cat("This fails\n")},
    timeout = 0.1
))
#> This fails
```
