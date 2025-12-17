# The FastRet GUI

Creates the FastRet GUI

## Usage

``` r
fastret_app(port = 8080, host = "0.0.0.0", reload = FALSE, nsw = 1)
```

## Arguments

- port:

  The port the application should listen on

- host:

  The address the application should listen on

- reload:

  Whether to reload the application when the source code changes

- nsw:

  The number of subworkers each worker is allowed to start. The higher
  this number, the faster individual tasks like model fitting can be
  processed.

## Value

An object of class `shiny.appobj`.

## Examples

``` r
x <- fastret_app()
if (interactive()) shiny::runApp(x)
```
