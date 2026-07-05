# Kuan speckle filter

The Kuan (1985) filter: like
[`rf_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee.md)
but derived without linearising the multiplicative noise model, giving
weight `W = (1 - cu^2/ci^2) / (1 + cu^2)` (clamped to `[0, 1]`).

## Usage

``` r
# S3 method for class 'Rcpp_GDALRaster'
rf_kuan(x, window = 7L, ...)

rf_kuan(x, ...)

# S3 method for class 'matrix'
rf_kuan(
  x,
  window = 7L,
  looks = 1,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_kuan(
  x,
  window = 7L,
  looks = 1,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_kuan(x, ...)

# S3 method for class 'SpatRaster'
rf_kuan(x, ...)
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

- looks:

  Single positive number: the effective number of looks (ENL) of the
  intensity image. Controls the assumed speckle strength
  `cu = 1 / sqrt(looks)`; single-look SAR intensity is `looks = 1`.

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

Kuan, D. T., Sawchuk, A. A., Strand, T. C., & Chavel, P. (1985).
Adaptive noise smoothing filter for images with signal-dependent noise.
*IEEE Transactions on Pattern Analysis and Machine Intelligence*, 7(2),
165-177.

## See also

[`rf_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee.md),
[`rf_gamma_map()`](https://belian-earth.github.io/rustyfilters/reference/rf_gamma_map.md)

## Examples

``` r
op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
rf_plot(s1_sthelens, main = "Sentinel-1 VV")
rf_plot(rf_kuan(s1_sthelens, window = 7L), main = "Kuan")

par(op)
```
