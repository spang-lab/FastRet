# Clip predictions to observed range

Clips predicted retention times by fitting a log-normal distribution to
the observed training RTs and bounding predictions to the central 99.99%
interval. All observed RTs must be positive to estimate the
distribution. If the estimated lower bound would be negative, it is
replaced by 1% of the observed minimum RT instead.

## Usage

``` r
clip_predictions(yhat, y)
```

## Arguments

- yhat:

  Numeric vector of predicted retention times.

- y:

  Numeric vector of observed retention times used to derive bounds.

## Value

Numeric vector of clipped (bounded) predictions.

## Examples

``` r
# Draw only a few samples (10) and clip based on these. The allowed range will
# be much bigger than the observed range.

set.seed(42)
y <- rlnorm(n = 1000, meanlog = 2, sdlog = 0.1)
yhat <- y
yhat[1] <- -100 # way too low to be realistic
yhat[2] <- 1000 # way too high to be realistic
yhat <- clip_predictions(yhat, y)
range(y)  # [ 6.18,  8.93]
#> [1]  5.274195 10.480647
yhat[1:2] # [ 4.96, 10.61] # Limited by theoretical bounds
#> [1]  5.274195 10.480647


# Draw more samples (1000) and clip based on these. The allowed range will
# be almost identical to the observed range.

set.seed(42)
y <- rnorm(n = 100, mean = 100, sd = 5)
yhat <- y
yhat[1] <- -100
yhat[2] <- 1000
yhat <- clip_predictions(yhat, y)
range(y)  # 83.14, 117.47
#> [1]  85.03455 111.43323
yhat[1:2] # 83.14, 117.72
#> [1]  84.06834 119.01241
```
