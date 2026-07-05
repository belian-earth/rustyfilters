# Guided filter

The guided filter (He et al. 2013) fits a local linear model between a
guide image and the input in every window, then averages the
coefficients: an edge-preserving smoother whose cost is independent of
the window size (it runs entirely on the package's O(1) box-mean
engine). With the default self-guidance (`guide = NULL`) it behaves like
a fast bilateral alternative; with a separate `guide` it transfers the
guide's structure onto `x` (e.g. detail-preserving smoothing of a noisy
band guided by a clean one).

## Usage

``` r
rf_guided(x, ...)

# S3 method for class 'matrix'
rf_guided(
  x,
  guide = NULL,
  window = 5L,
  eps = NULL,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_guided(
  x,
  guide = NULL,
  window = 5L,
  eps = NULL,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_guided(x, ...)

# S3 method for class 'Rcpp_GDALRaster'
rf_guided(x, guide = NULL, window = 5L, ...)

# S3 method for class 'SpatRaster'
rf_guided(x, guide = NULL, ...)
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

- guide:

  A matrix or 3-D array with the same dimensions as `x`, or `NULL` (the
  default) to use `x` itself.

- window:

  Window size: a single odd integer or a `c(rows, cols)` pair.

- eps:

  Single positive number: the regularisation added to the window
  variance of the guide. Windows whose variance is well below `eps` are
  smoothed; windows well above it (edges) are preserved.
  Scale-dependent: the default `NULL` uses `(0.1 * sd(guide))^2`.

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

Missing cells in either `x` or the guide are excluded from every window
fit; cells whose centre is `NA` in the guide stay `NA` in the output.

## References

He, K., Sun, J., & Tang, X. (2013). Guided image filtering. *IEEE
Transactions on Pattern Analysis and Machine Intelligence*, 35(6),
1397-1409.

## See also

[`rf_bilateral()`](https://belian-earth.github.io/rustyfilters/reference/rf_bilateral.md),
[`rf_mean()`](https://belian-earth.github.io/rustyfilters/reference/rf_mean.md),
[`rf_set_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)

## Examples

``` r
noisy <- volcano + matrix(rnorm(length(volcano), sd = 8), nrow(volcano))
op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
rf_plot(noisy, main = "noisy volcano")
rf_plot(rf_guided(noisy, window = 5L, eps = 64), main = "guided")

par(op)
```
