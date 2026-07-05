# Frost speckle filter

The Frost (1982) filter convolves each window with a damped exponential
kernel `exp(-B * d)` where `d` is the distance from the centre and
`B = damping * ci^2` adapts to the local variation coefficient:
homogeneous windows average broadly, heterogeneous windows concentrate
weight on the centre pixel.

## Usage

``` r
# S3 method for class 'Rcpp_GDALRaster'
rf_frost(x, window = 7L, ...)

rf_frost(x, ...)

# S3 method for class 'matrix'
rf_frost(
  x,
  window = 7L,
  damping = 2,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_frost(
  x,
  window = 7L,
  damping = 2,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_frost(x, ...)

# S3 method for class 'SpatRaster'
rf_frost(x, ...)
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

- damping:

  Single positive number scaling the kernel decay. Larger values
  preserve more edges and smooth less; `2` is a common default.

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

## References

Frost, V. S., Stiles, J. A., Shanmugan, K. S., & Holtzman, J. C. (1982).
A model for radar images and its application to adaptive digital
filtering of multiplicative noise. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 4(2), 157-166.

## See also

[`rf_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee.md),
[`rf_enhanced_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_enhanced_lee.md),
[`rf_gamma_map()`](https://belian-earth.github.io/rustyfilters/reference/rf_gamma_map.md)

## Examples

``` r
op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
rf_plot(s1_sthelens, main = "Sentinel-1 VV")
rf_plot(rf_frost(s1_sthelens, window = 7L, damping = 2), main = "Frost")

par(op)
```
