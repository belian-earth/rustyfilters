# Moving-window mean (boxcar) filter

Smooths by replacing each cell with the mean of its window. Runs on the
separable sliding-sum engine, so cost per cell is independent of the
window size.

## Usage

``` r
# S3 method for class 'Rcpp_GDALRaster'
rf_mean(x, window = 3L, ...)

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

  A numeric matrix or 3-D array (filtered layer by layer). Methods are
  also provided for terra `SpatRaster` objects (when terra is installed)
  and for open gdalraster `GDALRaster` datasets (when gdalraster is
  installed). `GDALRaster` methods return a new `GDALRaster` object open
  in update mode on a Float64 dataset with the source's geometry. Small
  datasets are filtered in memory and land on an in-memory `/vsimem`
  GTiff by default; datasets whose decoded size exceeds
  `options(rustyfilters.block_memory)` (default 2 GiB) stream through
  full-width row bands with a halo sized to the filter's window, writing
  to a GeoTIFF tempfile instead. Interior band seams are exact (the halo
  supplies the true neighbouring data; `edge` fires only at real raster
  edges). `GDALRaster` methods accept three extra arguments: `filename`
  (output path, replacing the tempfile/`/vsimem` default), `by_block`
  (`TRUE`/`FALSE` to force or forbid streaming) and `block_rows` (rows
  per band, sized from the memory budget by default).

- window:

  Window size in cells: a single odd positive integer, or a pair
  `c(rows, cols)` of odd positive integers.

- ...:

  Passed on to methods.

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
op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
rf_plot(volcano, main = "volcano")
rf_plot(rf_mean(volcano, window = 7L), main = "7 x 7 mean")

par(op)
```
