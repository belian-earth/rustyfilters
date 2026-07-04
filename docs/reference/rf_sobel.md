# Sobel edge detection

Convolves with the 3 x 3 Sobel kernels to estimate the local gradient:
`"x"` is the rate of change across columns, `"y"` across rows, and
`"magnitude"` (the default) is `sqrt(gx^2 + gy^2)`.

## Usage

``` r
rf_sobel(x, ...)

# S3 method for class 'matrix'
rf_sobel(
  x,
  direction = c("magnitude", "x", "y"),
  edge = c("nearest", "shrink", "reflect", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_sobel(
  x,
  direction = c("magnitude", "x", "y"),
  edge = c("nearest", "shrink", "reflect", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_sobel(x, ...)

# S3 method for class 'Rcpp_GDALRaster'
rf_sobel(x, ...)

# S3 method for class 'SpatRaster'
rf_sobel(x, ...)
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

- direction:

  `"magnitude"`, `"x"` or `"y"`.

- edge:

  Edge policy, as in
  [`rf_convolve()`](https://belian-earth.github.io/rustyfilters/reference/rf_convolve.md).
  Defaults to `"nearest"` for gradients: replicating the border cell
  gives a near-zero gradient there, whereas a shrinking window would
  fabricate strong edges.

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

[`rf_convolve()`](https://belian-earth.github.io/rustyfilters/reference/rf_convolve.md),
[`rf_laplacian()`](https://belian-earth.github.io/rustyfilters/reference/rf_laplacian.md)

## Examples

``` r
m <- matrix(rep(c(1, 1, 5, 5), each = 5), 5)
rf_sobel(m)
#>      [,1] [,2] [,3] [,4]
#> [1,]    0   16   16    0
#> [2,]    0   16   16    0
#> [3,]    0   16   16    0
#> [4,]    0   16   16    0
#> [5,]    0   16   16    0
rf_sobel(m, direction = "x")
#>      [,1] [,2] [,3] [,4]
#> [1,]    0   16   16    0
#> [2,]    0   16   16    0
#> [3,]    0   16   16    0
#> [4,]    0   16   16    0
#> [5,]    0   16   16    0
```
