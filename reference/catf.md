# catf function

Prints a formatted string with optional prefix and end strings.

## Usage

``` r
catf(
  ...,
  prefix = .Options$FastRet.catf.prefix,
  end = .Options$FastRet.catf.end
)
```

## Arguments

- ...:

  Arguments to be passed to sprintf for string formatting.

- prefix:

  A function returning a string to be used as the prefix. Default is a
  timestamp.

- end:

  A string to be used as the end of the message. Default is a newline
  character.

## Value

No return value. This function is called for its side effect of printing
a message.

## Examples

``` r
catf("Hello, %s!", "world")
#> 2025-12-17 23:27:08.18 Hello, world!
catf("Goodbye", prefix = NULL, end = "!\n")
#> 2025-12-17 23:27:08.18 Goodbye!
```
