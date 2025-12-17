# Execute an expression while redirecting output to a file

Execute an expression while redirecting output to a file

## Usage

``` r
withSink(expr, logfile = tempfile(fileext = ".txt"))
```

## Arguments

- expr:

  The expression to execute

- logfile:

  The file to redirect output to. Default is "tmp.txt".

## Value

The result of the expression

## Examples

``` r
logfile <- tempfile(fileext = ".txt")
withSink(logfile = logfile, expr = {
  cat("Helloworld\n")
  message("Goodbye")
})
#> Goodbye
readLines(logfile) == c("Helloworld", "Goodbye")
#> [1]  TRUE FALSE
```
