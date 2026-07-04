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
rf_bilateral(x, ...)

# S3 method for class 'SpatRaster'
rf_bilateral(x, ...)
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
set.seed(1)
m <- matrix(rep(c(0, 10), each = 50), 10) + rnorm(100)
rf_bilateral(m, sigma_d = 1.5, sigma_r = 3)
#>               [,1]         [,2]         [,3]        [,4]        [,5]      [,6]
#>  [1,]  0.143146839  0.433277981  0.456661594  0.48432052  0.33117897 10.295527
#>  [2,]  0.091220595  0.178886144  0.266057685  0.22151094  0.22038602 10.031554
#>  [3,] -0.135325785 -0.103181079  0.001907338  0.07393837  0.16563127 10.060620
#>  [4,]  0.137387958 -0.396204830 -0.365987489 -0.13390398 -0.01592159  9.821866
#>  [5,]  0.022827722  0.007365092 -0.129788085 -0.33477032 -0.26430912 10.206756
#>  [6,] -0.007906971 -0.040124131 -0.138907906 -0.22283801 -0.24704180 10.321837
#>  [7,]  0.180719365  0.058077761 -0.043296829 -0.10035536 -0.01882685  9.950505
#>  [8,]  0.285379566  0.215151386  0.012208003  0.08098876  0.21064326  9.863862
#>  [9,]  0.331861436  0.294488014  0.196290635  0.31540896  0.30453343 10.203144
#> [10,]  0.316862382  0.347542676  0.354019174  0.42084108  0.51903322 10.207370
#>            [,7]      [,8]      [,9]     [,10]
#>  [1,] 10.521672 10.256929 10.118779 10.121320
#>  [2,] 10.145819 10.098389 10.184529 10.381397
#>  [3,] 10.118693 10.149121 10.298332 10.443616
#>  [4,]  9.984859  9.925367  9.983365 10.388405
#>  [5,]  9.880397  9.855923 10.156033 10.404318
#>  [6,] 10.007047 10.015495 10.080259 10.173837
#>  [7,]  9.750834  9.939104 10.059597  9.847196
#>  [8,] 10.194989 10.007926  9.911676  9.788703
#>  [9,] 10.142247 10.059224  9.944596  9.699547
#> [10,] 10.431933 10.053900  9.954558  9.742056
```
