# Improved Lee sigma speckle filter

The improved sigma filter of Lee et al. (2009), as popularised by ESA
SNAP. It addresses the classic sigma filter's dark bias by (1)
estimating an a priori mean with a small MMSE filter and using published
sigma range bounds around it, (2) preserving point targets, detected as
clusters of pixels above the scene's 98th percentile, and (3) filtering
the in-range pixels with an MMSE weight based on a revised speckle
variation coefficient.

## Usage

``` r
# S3 method for class 'Rcpp_GDALRaster'
rf_lee_sigma_improved(
  x,
  window = 7L,
  looks = 1,
  sigma = 0.9,
  target_window = 3L,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

rf_lee_sigma_improved(x, ...)

# S3 method for class 'matrix'
rf_lee_sigma_improved(
  x,
  window = 7L,
  looks = 1,
  sigma = 0.9,
  target_window = 3L,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_lee_sigma_improved(
  x,
  window = 7L,
  looks = 1,
  sigma = 0.9,
  target_window = 3L,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_lee_sigma_improved(x, ...)

# S3 method for class 'SpatRaster'
rf_lee_sigma_improved(x, ...)
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

- looks:

  Effective number of looks. Must be 1, 2, 3 or 4; the published sigma
  range tables only cover these.

- sigma:

  Sigma confidence level: one of 0.5, 0.6, 0.7, 0.8 or 0.9 (the default,
  covering the widest speckle range).

- target_window:

  Single odd integer (typically 3 or 5): the window used for the a
  priori mean estimate and for point-target detection.

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

- ...:

  Passed on to methods.

## Value

An object of the same class and dimensions as `x` (dimnames are
preserved), containing the filtered values as doubles.

## Details

Point-target marking is order-independent here: a pixel is preserved
when it exceeds the layer's 98th percentile and lies within
`target_window` of a cluster of more than five such pixels. SNAP marks
clusters during a sequential scan, which can differ at cluster fringes.
The 98th percentile is computed per layer over valid cells.

## References

Lee, J.-S., Wen, J.-H., Ainsworth, T. L., Chen, K.-S., & Chen, A. J.
(2009). Improved sigma filter for speckle filtering of SAR imagery.
*IEEE Transactions on Geoscience and Remote Sensing*, 47(1), 202-213.

## See also

[`rf_lee_sigma()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee_sigma.md)
for the classic 1983 filter,
[`rf_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee.md),
[`rf_enhanced_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_enhanced_lee.md)

## Examples

``` r
op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
rf_plot(s1_sthelens, main = "Sentinel-1 VV")
rf_plot(
  rf_lee_sigma_improved(s1_sthelens, window = 7L, sigma = 0.9),
  main = "improved Lee sigma"
)

par(op)
```
