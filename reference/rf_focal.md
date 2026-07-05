# Focal statistics over a moving window

Computes a summary statistic of the cells in a moving window centred on
each cell. All statistics share the same multi-threaded engine; see
[`rf_mean()`](https://belian-earth.github.io/rustyfilters/reference/rf_mean.md),
[`rf_median()`](https://belian-earth.github.io/rustyfilters/reference/rf_median.md)
and
[`rf_gaussian()`](https://belian-earth.github.io/rustyfilters/reference/rf_gaussian.md)
for the dedicated smoothing filters and
[`rf_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee.md)
and friends for adaptive speckle filters.

## Usage

``` r
rf_focal(x, ...)

# S3 method for class 'matrix'
rf_focal(
  x,
  window = 3L,
  stat = c("mean", "median", "min", "max", "range", "sd", "sum", "mode"),
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_focal(
  x,
  window = 3L,
  stat = c("mean", "median", "min", "max", "range", "sd", "sum", "mode"),
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_focal(x, ...)

# S3 method for class 'Rcpp_GDALRaster'
rf_focal(x, window = 3L, ...)

# S3 method for class 'SpatRaster'
rf_focal(x, ...)
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

- ...:

  Passed on to methods.

- window:

  Window size in cells: a single odd positive integer, or a pair
  `c(rows, cols)` of odd positive integers.

- stat:

  Statistic to compute for each window. One of `"mean"`, `"median"`,
  `"min"`, `"max"`, `"range"` (max minus min), `"sd"` (sample standard
  deviation), `"sum"` or `"mode"`.

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

## Details

`"mode"` returns the most frequent value using exact floating-point
equality, with ties resolved to the lowest value. It is intended for
categorical data encoded as numbers and is not meaningful for continuous
data.

## See also

[`rf_mean()`](https://belian-earth.github.io/rustyfilters/reference/rf_mean.md),
[`rf_median()`](https://belian-earth.github.io/rustyfilters/reference/rf_median.md),
[`rf_gaussian()`](https://belian-earth.github.io/rustyfilters/reference/rf_gaussian.md),
[`rf_set_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)

## Examples

``` r
op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
rf_plot(rf_focal(volcano, window = 5L, stat = "sd"), main = "focal sd")
rf_plot(rf_focal(volcano, window = 5L, stat = "range"), main = "focal range")

par(op)
```
