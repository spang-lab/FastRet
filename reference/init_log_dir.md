# Initialize log directory

Initializes the log directory for the session. It creates a new
directory if it does not exist.

## Usage

``` r
init_log_dir(SE)
```

## Arguments

- SE:

  A list containing session information.

## Value

Updates the logdir element in the SE list with the path to the log
directory.

## Examples

``` r
SE <- as.environment(list(session = list(token = "asdf")))
init_log_dir(SE)
#> 2026-06-11 20:22:15.99 Start: init_log_dir
#> 2026-06-11 20:22:15.99 Logdir: /tmp/Rtmpu0yfme/FastRet/asdf
dir.exists(SE$logdir)
#> [1] TRUE
```
