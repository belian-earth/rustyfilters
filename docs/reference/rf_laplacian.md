# Laplacian filter

Convolves with the discrete Laplacian, highlighting local extrema and
rapid changes (zero over flat and linearly varying regions).

## Usage

``` r
rf_laplacian(x, ...)

# S3 method for class 'matrix'
rf_laplacian(
  x,
  neighbours = 4L,
  edge = c("nearest", "shrink", "reflect", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_laplacian(
  x,
  neighbours = 4L,
  edge = c("nearest", "shrink", "reflect", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_laplacian(x, ...)

# S3 method for class 'Rcpp_GDALRaster'
rf_laplacian(x, ...)

# S3 method for class 'SpatRaster'
rf_laplacian(x, ...)
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

- neighbours:

  `4L` for the cross-shaped kernel (`c(0, 1, 0, 1, -4, 1, 0, 1, 0)`),
  `8L` to include diagonals.

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
[`rf_sobel()`](https://belian-earth.github.io/rustyfilters/reference/rf_sobel.md)

## Examples

``` r
m <- matrix(as.numeric(1:25), 5)
rf_laplacian(m) # linear ramp: zero in the interior
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    6    1    1    1   -4
#> [2,]    5    0    0    0   -5
#> [3,]    5    0    0    0   -5
#> [4,]    5    0    0    0   -5
#> [5,]    4   -1   -1   -1   -6
```
