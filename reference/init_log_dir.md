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
#> 2025-12-17 23:27:09.36 Start: init_log_dir
#> 2025-12-17 23:27:09.36 Logdir: /tmp/RtmpeHumg1/FastRet/asdf
dir.exists(SE$logdir)
#> [1] TRUE
```
