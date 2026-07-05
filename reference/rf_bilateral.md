# Bilateral filter

Smooths with weights that are the product of a spatial Gaussian
(distance from the centre cell, `sigma_d`) and a range Gaussian
(difference from the centre value, `sigma_r`), so averaging happens
within regions of similar value while sharp transitions survive.

## Usage

``` r
rf_bilateral(x, ...)

# S3 method for class 'matrix'
rf_bilateral(
  x,
  sigma_d = 1.5,
  sigma_r = NULL,
  window = NULL,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_bilateral(
  x,
  sigma_d = 1.5,
  sigma_r = NULL,
  window = NULL,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_bilateral(x, ...)

# S3 method for class 'Rcpp_GDALRaster'
rf_bilateral(x, sigma_d = 1.5, sigma_r = NULL, window = NULL, ...)

# S3 method for class 'SpatRaster'
rf_bilateral(x, ...)
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

- sigma_d:

  Single positive number: the spatial standard deviation in cells.

- sigma_r:

  Single positive number: the range standard deviation in value units.
  Values differing from the centre by much more than `sigma_r` are
  effectively excluded. The default `NULL` uses the standard deviation
  of the valid cells of `x`, a serviceable starting point that you
  should expect to tune.

- window:

  Window size: a single odd integer or a `c(rows, cols)` pair. The
  default `NULL` uses `2 * ceiling(2 * sigma_d) + 1`.

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

Cells whose centre value is `NA` stay `NA` (the range weight needs the
centre). Cost grows with the window area; this is the slowest smoother
in the package.

## See also

[`rf_guided()`](https://belian-earth.github.io/rustyfilters/reference/rf_guided.md),
[`rf_gaussian()`](https://belian-earth.github.io/rustyfilters/reference/rf_gaussian.md),
[`rf_median()`](https://belian-earth.github.io/rustyfilters/reference/rf_median.md),
[`rf_set_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)

## Examples

``` r
noisy <- volcano + matrix(rnorm(length(volcano), sd = 8), nrow(volcano))
op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
rf_plot(noisy, main = "noisy volcano")
rf_plot(rf_bilateral(noisy, sigma_d = 2, sigma_r = 20), main = "bilateral")

par(op)
```
