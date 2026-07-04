# Lee speckle filter

The Lee (1980) minimum mean-square-error filter: each pixel becomes
`m + W * (x - m)` with weight `W = v / (v + (m * cu)^2)`, so homogeneous
areas take the window mean while high-contrast features stay close to
the observed value.

## Usage

``` r
# S3 method for class 'Rcpp_GDALRaster'
rf_lee(x, ...)

rf_lee(x, ...)

# S3 method for class 'matrix'
rf_lee(
  x,
  window = 7L,
  looks = 1,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_lee(
  x,
  window = 7L,
  looks = 1,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_lee(x, ...)

# S3 method for class 'SpatRaster'
rf_lee(x, ...)
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

- window:

  Window size in cells: a single odd positive integer, or a pair
  `c(rows, cols)` of odd positive integers.

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

## Details

Cells whose centre value is `NA` remain `NA`. Speckle filters expect
intensity (linear power) data, not dB.

## References

Lee, J.-S. (1980). Digital image enhancement and noise filtering by use
of local statistics. *IEEE Transactions on Pattern Analysis and Machine
Intelligence*, 2(2), 165-168.

## See also

[`rf_enhanced_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_enhanced_lee.md),
[`rf_lee_sigma()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee_sigma.md),
[`rf_kuan()`](https://belian-earth.github.io/rustyfilters/reference/rf_kuan.md),
[`rf_frost()`](https://belian-earth.github.io/rustyfilters/reference/rf_frost.md),
[`rf_gamma_map()`](https://belian-earth.github.io/rustyfilters/reference/rf_gamma_map.md),
[`rf_set_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)

## Examples

``` r
set.seed(1)
truth <- matrix(rep(c(1, 4), each = 200), 20, 20)
speckled <- truth * rexp(400)
rf_lee(speckled, window = 7L, looks = 1)
#>            [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#>  [1,] 0.8254263 1.3339167 0.8991960 0.5762854 0.7559167 1.2044998 0.6995902
#>  [2,] 1.0895728 0.8006348 0.9926381 1.0154354 1.0035295 0.8619246 0.9142819
#>  [3,] 0.5215517 0.6097852 1.1434416 0.5971624 0.8921337 0.7910391 1.3112362
#>  [4,] 0.5898735 0.8102551 1.1373611 1.2889670 0.5832238 0.8304846 1.0490610
#>  [5,] 0.7767612 0.6441977 0.8594833 2.8237619 1.6418123 1.2284343 0.7203543
#>  [6,] 2.0312064 0.6155725 0.7387520 0.7935470 0.9838012 1.7802188 0.7603482
#>  [7,] 1.2040862 0.9371920 1.2426455 1.8606519 1.3561404 0.8656322 2.0144827
#>  [8,] 0.9191058 2.3383192 1.1799027 1.1966117 1.2209796 1.6471567 1.3045747
#>  [9,] 1.1179283 1.2001134 1.0115167 1.0830602 1.2460165 1.0341523 0.8565955
#> [10,] 0.9556425 1.2326298 1.5132687 1.1190706 1.5579403 1.1659480 1.2230777
#> [11,] 1.3849773 1.3624107 0.9461692 1.4254361 1.1346938 1.3334932 0.8563818
#> [12,] 1.0362603 0.7289508 1.5440176 1.5762276 0.8364235 1.3926383 1.1455585
#> [13,] 1.2223518 0.7934960 1.9802172 1.8312370 1.0387643 0.7625736 0.8407044
#> [14,] 2.6170895 1.1759940 0.8126226 0.6356924 0.6188067 0.6680445 0.7442176
#> [15,] 1.0801318 0.5910915 0.7577630 0.6145873 0.7626410 1.1612058 0.5645948
#> [16,] 1.0282944 0.9649555 0.9066689 0.3908145 0.5198399 0.7918605 0.8515711
#> [17,] 1.3813977 0.5333055 0.4551102 0.5253954 0.5692359 0.6161855 0.4831838
#> [18,] 0.7325320 0.7056485 0.5189477 1.0517799 0.4527911 0.4504600 0.4445028
#> [19,] 0.6418525 0.7297778 0.8521990 0.7365715 0.3996614 0.4800134 1.3409034
#> [20,] 0.7374739 0.5083262 0.7073876 1.6831391 0.4868060 0.5009294 0.7371782
#>            [,8]      [,9]     [,10]    [,11]     [,12]    [,13]     [,14]
#>  [1,] 0.9018979 1.5369430 1.6149860 1.705047  3.909392 5.147485  5.555309
#>  [2,] 1.0908063 1.3010802 1.2794305 3.469697  2.368527 4.219339  7.676861
#>  [3,] 1.5808587 1.5634661 1.1425374 2.621423  3.365061 4.441263  2.180150
#>  [4,] 1.2540020 1.6028342 1.5448769 2.949348  3.720137 6.608257  3.388991
#>  [5,] 0.8108741 1.6394532 1.7774634 1.584134  4.740764 2.954709  3.469886
#>  [6,] 1.0751089 0.9616217 1.2868803 2.601542  1.913667 1.735381  2.149932
#>  [7,] 1.0299085 2.2605277 1.4364371 5.757561  5.019611 2.011691  3.799724
#>  [8,] 0.9179626 1.9383760 1.1709919 2.033270  2.633220 7.107994  2.390445
#>  [9,] 1.6111618 0.9580479 1.0217878 2.146261  2.357602 2.887276  6.461367
#> [10,] 1.8198703 2.4339584 2.4540914 3.618090  4.075666 2.716552  4.077627
#> [11,] 0.9692718 1.2264733 1.5648134 1.462016  2.541874 9.316163  2.770604
#> [12,] 0.8299630 1.8099492 1.2722415 4.397006 11.834226 3.104220  5.088646
#> [13,] 1.5052679 1.0199605 0.8364470 3.794480  1.774884 3.603689  2.577210
#> [14,] 1.9754748 2.0417582 1.5439178 1.386460  2.466614 1.369991  1.593268
#> [15,] 1.6659180 1.4488667 0.9869748 3.991347  2.726107 2.015304  6.178107
#> [16,] 0.9063158 1.4943869 1.2659960 1.311416  2.166094 3.224441 16.024441
#> [17,] 1.5967210 1.7194583 1.2223786 2.690965  3.437266 4.805081 10.604136
#> [18,] 0.8439719 1.0034277 1.0522077 6.976863  2.240748 4.415883  2.394947
#> [19,] 0.3365112 0.8093863 1.2034136 5.265216  2.922305 1.789366  4.862858
#> [20,] 0.3219682 1.0398453 0.8569460 1.313759  3.713123 2.365485  3.981152
#>          [,15]    [,16]    [,17]    [,18]    [,19]    [,20]
#>  [1,] 3.897689 3.273708 4.285779 2.880584 3.999557 6.353870
#>  [2,] 7.559955 8.465115 3.192286 5.509295 2.737308 5.018002
#>  [3,] 2.912337 5.920349 2.558429 2.804465 3.594903 2.409643
#>  [4,] 2.843298 3.941525 2.634617 3.828622 2.703895 2.815097
#>  [5,] 3.940523 5.375627 3.449827 4.332103 3.778145 2.212409
#>  [6,] 2.571497 2.325278 1.868592 4.563882 3.585596 7.415677
#>  [7,] 2.178238 4.203443 3.133693 2.763093 3.138550 2.000922
#>  [8,] 2.723087 2.271655 6.552451 2.095860 3.164181 1.736236
#>  [9,] 2.587868 5.603196 2.806348 1.481645 5.861382 2.033965
#> [10,] 5.694943 6.295785 1.550652 3.520426 2.078744 2.294561
#> [11,] 2.193850 1.893608 2.015579 1.361353 1.926449 1.274603
#> [12,] 1.894057 2.742245 4.168659 1.740233 2.016111 2.699904
#> [13,] 4.108229 1.616690 3.099882 2.406157 5.630050 1.664240
#> [14,] 2.264809 2.197400 2.625869 4.131039 2.719077 1.596706
#> [15,] 3.347375 7.554197 4.310597 4.844876 4.596731 4.328201
#> [16,] 2.714283 9.179409 2.282733 4.081415 2.999476 4.203739
#> [17,] 4.300395 3.378790 2.636713 2.763602 3.167133 5.879195
#> [18,] 5.178891 3.679545 5.224450 7.052358 4.473423 3.335320
#> [19,] 4.821004 2.761963 7.171890 2.425859 4.480103 4.813084
#> [20,] 3.188894 4.802245 2.523299 4.671834 2.412165 2.284633
```
