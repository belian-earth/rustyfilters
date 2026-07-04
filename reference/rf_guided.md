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
rf_guided(x, guide = NULL, ...)

# S3 method for class 'SpatRaster'
rf_guided(x, guide = NULL, ...)
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
set.seed(1)
m <- matrix(rep(c(0, 10), each = 50), 10) + rnorm(100)
rf_guided(m, window = 5L)
#>             [,1]        [,2]        [,3]         [,4]        [,5]      [,6]
#>  [1,] -0.4446212  1.27072828  0.81915110  1.258741425 -0.09268275 10.336216
#>  [2,]  0.1581293  0.34100630  0.69351313 -0.055952437 -0.18093902  9.385098
#>  [3,] -0.6339756 -0.49767087  0.08143414  0.376709599  0.70735758 10.282630
#>  [4,]  1.2446933 -1.86280580 -1.72494781 -0.029540381  0.57286826  8.886872
#>  [5,]  0.2399302  0.91670647  0.53435449 -1.211203523 -0.59935720 11.313911
#>  [6,] -0.5947828 -0.03156511 -0.03566379 -0.342573282 -0.60997570 11.830191
#>  [7,]  0.3546261  0.00950134 -0.10050241 -0.307979104  0.39321509  9.603210
#>  [8,]  0.5402274  0.74854161 -1.12690722 -0.003332007  0.77163404  8.964137
#>  [9,]  0.4318143  0.65172931 -0.31715571  0.987036663 -0.03314232 10.486859
#> [10,] -0.1162973  0.49221889  0.38392104  0.706207238  0.87862095  9.821295
#>            [,7]      [,8]      [,9]     [,10]
#>  [1,] 12.162229 10.426100  9.574170  9.648303
#>  [2,]  9.964366  9.418367  9.926463 10.972635
#>  [3,] 10.616740 10.537895 10.995960 10.937527
#>  [4,] 10.010225  9.198988  8.752920 10.569083
#>  [5,]  9.310022  8.922491 10.496821 11.242985
#>  [6,] 10.149220 10.238928 10.270652 10.431914
#>  [7,]  8.341597  9.605124 10.843838  9.033063
#>  [8,] 11.292477  9.973314  9.736392  9.563920
#>  [9,] 10.112988 10.032811 10.262659  9.109205
#> [10,] 11.925016  9.489169 10.173975  9.640103
rf_guided(m, window = 5L, eps = 100) # large eps: plain smoothing
#>               [,1]      [,2]      [,3]     [,4]     [,5]     [,6]     [,7]
#>  [1,]  0.105351301 0.6065634 1.1841527 2.249748 3.441085 6.903453 8.534911
#>  [2,]  0.052944913 0.5016732 1.1151475 2.022579 3.371906 6.737767 8.230899
#>  [3,] -0.003095555 0.4130992 1.0161606 2.031545 3.467631 6.847381 8.281176
#>  [4,] -0.050509007 0.2619536 0.7713559 1.889494 3.364442 6.562651 8.117100
#>  [5,] -0.048235247 0.4263132 0.9745021 1.715474 3.163539 6.924368 8.003294
#>  [6,] -0.049637788 0.3832513 0.9295146 1.840340 3.179042 7.010997 8.108079
#>  [7,]  0.029941146 0.4500084 0.9809800 1.903756 3.385856 6.671052 7.870406
#>  [8,]  0.124152748 0.5761430 0.9684115 2.029957 3.520297 6.613830 8.275422
#>  [9,]  0.180734462 0.6293627 1.1008557 2.220324 3.457589 6.875588 8.164415
#> [10,]  0.204137897 0.6478031 1.1982541 2.232098 3.645155 6.793225 8.408890
#>           [,8]     [,9]     [,10]
#>  [1,] 9.247047 9.753216 10.233378
#>  [2,] 9.142860 9.762101 10.243368
#>  [3,] 9.216490 9.797603 10.220200
#>  [4,] 9.017092 9.597677 10.124401
#>  [5,] 8.969579 9.666082 10.100947
#>  [6,] 9.072037 9.634008 10.053059
#>  [7,] 8.983899 9.615356  9.965143
#>  [8,] 9.021042 9.546371  9.940272
#>  [9,] 9.048812 9.588906  9.940903
#> [10,] 9.006235 9.577931  9.929005
```
