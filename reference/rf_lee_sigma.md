# Lee sigma speckle filter

The classic two-sigma Lee filter (1983): each pixel is replaced by the
mean of the window pixels lying within `x * (1 +/- k * cu)`, the range a
pure-speckle observation of the centre value would span. If fewer than
`min_count` pixels qualify, the full window mean is used instead
(suppressing isolated dark/bright noise).

## Usage

``` r
# S3 method for class 'Rcpp_GDALRaster'
rf_lee_sigma(x, window = 7L, ...)

rf_lee_sigma(x, ...)

# S3 method for class 'matrix'
rf_lee_sigma(
  x,
  window = 7L,
  looks = 1,
  k = 2,
  min_count = 3L,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_lee_sigma(
  x,
  window = 7L,
  looks = 1,
  k = 2,
  min_count = 3L,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_lee_sigma(x, ...)

# S3 method for class 'SpatRaster'
rf_lee_sigma(x, ...)
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

- k:

  Single positive number: the sigma multiplier defining the inclusion
  range. `2` covers roughly 95.5% of pure speckle.

- min_count:

  Single positive integer: minimum number of in-range pixels below which
  the filter falls back to the full window mean.

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

Lee, J.-S. (1983). Digital image smoothing and the sigma filter.
*Computer Vision, Graphics, and Image Processing*, 24(2), 255-269.

## See also

[`rf_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee.md),
[`rf_enhanced_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_enhanced_lee.md)

## Examples

``` r
op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
rf_plot(s1_sthelens, main = "Sentinel-1 VV")
rf_plot(rf_lee_sigma(s1_sthelens, window = 7L), main = "Lee sigma")

par(op)
```
