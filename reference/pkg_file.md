# Get package file

Returns the path to a file within the FastRet package.

## Usage

``` r
pkg_file(path, mustWork = FALSE)
```

## Arguments

- path:

  The path to the file within the package.

- mustWork:

  If TRUE, an error is thrown if the file does not exist.

## Value

The path to the file.

## Examples

``` r
path <- pkg_file("extdata/RP.xlsx")
```
