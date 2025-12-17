# Try expression with predefined error message

Executes an expression and prints an error message if it fails

## Usage

``` r
withStopMessage(expr)
```

## Arguments

- expr:

  The expression to execute

## Value

The result of the expression

## Examples

``` r
f <- function(expr) {
  val <- try(expr, silent = TRUE)
  err <- if (inherits(val, "try-error")) attr(val, "condition") else NULL
  if (!is.null(err)) value <- NULL
  list(value = val, error = err)
}
ret <- f(log("a")) # this error will not show up in the console
ret <- f(withStopMessage(log("a"))) # this error will show up in the console
#> Error in log("a") : non-numeric argument to mathematical function
```
