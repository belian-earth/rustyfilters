# Moving-window median filter

Replaces each cell with the median of its window: a robust smoother that
preserves edges and removes salt-and-pepper noise. Windows with an even
number of valid cells (possible at edges or around missing values)
average the two middle values, matching
[`stats::median()`](https://rdrr.io/r/stats/median.html).

## Usage

``` r
rf_median(x, ...)

# S3 method for class 'matrix'
rf_median(
  x,
  window = 3L,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_median(
  x,
  window = 3L,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_median(x, ...)

# S3 method for class 'SpatRaster'
rf_median(x, ...)
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

[`rf_mean()`](https://belian-earth.github.io/rustyfilters/reference/rf_mean.md),
[`rf_gaussian()`](https://belian-earth.github.io/rustyfilters/reference/rf_gaussian.md),
[`rf_focal()`](https://belian-earth.github.io/rustyfilters/reference/rf_focal.md),
[`rf_set_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)

## Examples

``` r
m <- matrix(as.numeric(1:25), 5)
m[3, 3] <- 1000
rf_median(m)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  4.0  6.5 11.5 16.5 19.0
#> [2,]  4.5  7.0 12.0 18.0 19.5
#> [3,]  5.5  8.0 14.0 19.0 20.5
#> [4,]  6.5  9.0 15.0 20.0 21.5
#> [5,]  7.0  9.5 14.5 19.5 22.0
```
