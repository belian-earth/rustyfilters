# Moving-window mean (boxcar) filter

Smooths by replacing each cell with the mean of its window. Runs on the
separable sliding-sum engine, so cost per cell is independent of the
window size.

## Usage

``` r
rf_mean(x, ...)

# S3 method for class 'matrix'
rf_mean(
  x,
  window = 3L,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_mean(
  x,
  window = 3L,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_mean(x, ...)

# S3 method for class 'SpatRaster'
rf_mean(x, ...)
```

## Arguments

- x:

  A numeric matrix or 3-D array (filtered layer by layer). Methods for
  terra `SpatRaster` objects are provided when terra is installed.

- ...:

  Passed on to methods.

- window:

  Window size in cells: a single odd positive integer, or a pair
  `c(rows, cols)` of odd positive integers.

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

[`rf_median()`](https://belian-earth.github.io/rustyfilters/reference/rf_median.md),
[`rf_gaussian()`](https://belian-earth.github.io/rustyfilters/reference/rf_gaussian.md),
[`rf_focal()`](https://belian-earth.github.io/rustyfilters/reference/rf_focal.md)
for other window statistics,
[`rf_set_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)

## Examples

``` r
m <- matrix(as.numeric(1:25), 5)
rf_mean(m)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  4.0  6.5 11.5 16.5 19.0
#> [2,]  4.5  7.0 12.0 17.0 19.5
#> [3,]  5.5  8.0 13.0 18.0 20.5
#> [4,]  6.5  9.0 14.0 19.0 21.5
#> [5,]  7.0  9.5 14.5 19.5 22.0
rf_mean(m, window = c(3L, 5L), edge = "reflect")
#>          [,1]      [,2]     [,3]     [,4]     [,5]
#> [1,] 5.333333  7.333333 11.33333 15.33333 17.33333
#> [2,] 6.000000  8.000000 12.00000 16.00000 18.00000
#> [3,] 7.000000  9.000000 13.00000 17.00000 19.00000
#> [4,] 8.000000 10.000000 14.00000 18.00000 20.00000
#> [5,] 8.666667 10.666667 14.66667 18.66667 20.66667
```
