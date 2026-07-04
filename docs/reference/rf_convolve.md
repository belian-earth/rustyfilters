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
rf_convolve(x, ...)

# S3 method for class 'SpatRaster'
rf_convolve(x, ...)
```

## Arguments

- x:

  A numeric matrix or 3-D array (filtered layer by layer). Methods are
  also provided for terra `SpatRaster` objects (when terra is installed)
  and for open gdalraster `GDALRaster` datasets (when gdalraster is
  installed). `GDALRaster` methods read the dataset into memory, filter
  it, and return a new `GDALRaster` object open in update mode on a
  Float64 dataset with the source's geometry: an in-memory `/vsimem`
  GTiff by default, or pass `filename` to write to disk.

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
m <- matrix(as.numeric(1:25), 5)
sharpen <- matrix(c(0, -1, 0, -1, 5, -1, 0, -1, 0), 3)
rf_convolve(m, sharpen)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   -3   11   21   31   67
#> [2,]   -1    7   12   17   49
#> [3,]    1    8   13   18   51
#> [4,]    3    9   14   19   53
#> [5,]   11   21   31   41   81
rf_convolve(m, matrix(1, 3, 3) / 9, normalize = TRUE)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  4.0  6.5 11.5 16.5 19.0
#> [2,]  4.5  7.0 12.0 17.0 19.5
#> [3,]  5.5  8.0 13.0 18.0 20.5
#> [4,]  6.5  9.0 14.0 19.0 21.5
#> [5,]  7.0  9.5 14.5 19.5 22.0
```
