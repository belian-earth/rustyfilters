# Filter with an arbitrary kernel

Applies a user-supplied kernel over a moving window: each output cell is
the sum of the window values times the matching kernel weights. This
gives sharpening, embossing, gradients or any custom linear filter;
[`rf_sobel()`](https://belian-earth.github.io/rustyfilters/reference/rf_sobel.md)
and
[`rf_laplacian()`](https://belian-earth.github.io/rustyfilters/reference/rf_laplacian.md)
are thin wrappers with fixed kernels.

## Usage

``` r
rf_convolve(x, ...)

# S3 method for class 'matrix'
rf_convolve(
  x,
  kernel,
  normalize = FALSE,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_convolve(
  x,
  kernel,
  normalize = FALSE,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_convolve(x, ...)

# S3 method for class 'Rcpp_GDALRaster'
rf_convolve(x, kernel, ...)

# S3 method for class 'SpatRaster'
rf_convolve(x, ...)
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

- kernel:

  A numeric matrix with odd dimensions and no missing values. Applied
  as-is (cross-correlation): the kernel is not flipped, which only
  matters for asymmetric kernels.

- normalize:

  If `TRUE`, divide each result by the sum of the kernel weights
  actually used, so partial windows (at edges under `edge = "shrink"`,
  or around missing values) keep the input's scale. Sensible for
  smoothing kernels; leave `FALSE` for derivative kernels whose weights
  sum to zero.

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

Under `na_policy = "omit"`, missing cells simply drop out of the
weighted sum. For zero-sum kernels this biases results near missing
values; there is no principled general correction, so consider
`na_policy = "propagate"` for derivative kernels on gappy data.

## See also

[`rf_sobel()`](https://belian-earth.github.io/rustyfilters/reference/rf_sobel.md),
[`rf_laplacian()`](https://belian-earth.github.io/rustyfilters/reference/rf_laplacian.md),
[`rf_gaussian()`](https://belian-earth.github.io/rustyfilters/reference/rf_gaussian.md),
[`rf_set_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)

## Examples

``` r
sharpen <- matrix(c(0, -1, 0, -1, 5, -1, 0, -1, 0), 3)
op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
rf_plot(volcano, main = "volcano")
rf_plot(rf_convolve(volcano, sharpen), main = "sharpened")

par(op)
```
