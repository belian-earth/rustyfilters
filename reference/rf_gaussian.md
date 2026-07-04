# Gaussian blur

Smooths with a separable Gaussian kernel: two 1-D passes, so cost grows
with the window edge length rather than its area. Under
`na_policy = "omit"` this is normalized convolution: missing cells are
excluded and the remaining weights are rescaled, which also fills
isolated `NA` holes with their neighbourhood estimate.

## Usage

``` r
rf_gaussian(x, ...)

# S3 method for class 'matrix'
rf_gaussian(
  x,
  sigma = 1,
  window = NULL,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_gaussian(
  x,
  sigma = 1,
  window = NULL,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_gaussian(x, ...)

# S3 method for class 'SpatRaster'
rf_gaussian(x, ...)
```

## Arguments

- x:

  A numeric matrix or 3-D array (filtered layer by layer). Methods for
  terra `SpatRaster` objects are provided when terra is installed.

- ...:

  Passed on to methods.

- sigma:

  Single positive number: the Gaussian standard deviation in cells.

- window:

  Kernel size in cells: a single odd positive integer or a pair
  `c(rows, cols)`. The default `NULL` uses `2 * ceiling(3 * sigma) + 1`
  in both dimensions, which captures effectively all of the kernel mass.

- edge:

  How to treat windows that overhang the matrix edge: `"shrink"`
  (default) truncates the window to the cells that exist; `"reflect"`
  mirrors the matrix across the boundary; `"nearest"` repeats the
  closest edge cell; `"constant"` pads with `edge_value`.

- edge_value:

  Single number used to pad when `edge = "constant"`. `NA` is allowed
  and behaves as missing data under `na_policy = "omit"`.

- na_policy:

  How to treat missing values inside a window. The default `"omit"`
  excludes them from the statistics (like `na.rm = TRUE`); a window with
  no valid cells yields `NA`. `"propagate"` is the fast path: no
  per-cell NA handling is compiled into the inner loop, and any `NA` in
  a window makes the result `NA`. Use it when the input has no missing
  values (or when spreading `NA` is acceptable) for maximum speed.

## Value

An object of the same class and dimensions as `x` (dimnames are
preserved), containing the filtered values as doubles.

## See also

[`rf_mean()`](https://belian-earth.github.io/rustyfilters/reference/rf_mean.md),
[`rf_median()`](https://belian-earth.github.io/rustyfilters/reference/rf_median.md),
[`rf_focal()`](https://belian-earth.github.io/rustyfilters/reference/rf_focal.md),
[`rf_set_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)

## Examples

``` r
m <- matrix(as.numeric(1:25), 5)
rf_gaussian(m, sigma = 1)
#>          [,1]      [,2]     [,3]     [,4]     [,5]
#> [1,] 4.116513  7.163616 11.51942 15.87522 18.92232
#> [2,] 4.725934  7.773037 12.12884 16.48464 19.53175
#> [3,] 5.597094  8.644198 13.00000 17.35580 20.40291
#> [4,] 6.468255  9.515358 13.87116 18.22696 21.27407
#> [5,] 7.077675 10.124779 14.48058 18.83638 21.88349
rf_gaussian(m, sigma = 2, window = 7L, edge = "reflect")
#>          [,1]      [,2]     [,3]     [,4]     [,5]
#> [1,] 6.608405  8.642498 11.93473 15.22697 17.26106
#> [2,] 7.015223  9.049317 12.34155 15.63379 17.66788
#> [3,] 7.673670  9.707764 13.00000 16.29224 18.32633
#> [4,] 8.332118 10.366211 13.65845 16.95068 18.98478
#> [5,] 8.738936 10.773030 14.06527 17.35750 19.39160
```
