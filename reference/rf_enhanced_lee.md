# Enhanced Lee speckle filter

The enhanced Lee filter (Lopes et al. 1990) splits each window into
three regimes by its variation coefficient `ci`: homogeneous windows
(`ci <= cu`) take the window mean, heterogeneous windows blend mean and
centre with an exponential damping weight, and point targets
(`ci >= cmax = sqrt(1 + 2/looks)`) are preserved untouched.

## Usage

``` r
rf_enhanced_lee(x, ...)

# S3 method for class 'matrix'
rf_enhanced_lee(
  x,
  window = 7L,
  looks = 1,
  damping = 1,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_enhanced_lee(
  x,
  window = 7L,
  looks = 1,
  damping = 1,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_enhanced_lee(x, ...)

# S3 method for class 'SpatRaster'
rf_enhanced_lee(x, ...)
```

## Arguments

- x:

  A numeric matrix or 3-D array (filtered layer by layer). Methods for
  terra `SpatRaster` objects are provided when terra is installed.

- ...:

  Passed on to methods.

- window:

  Window size in cells: a single odd positive integer, or a pair
  `c(rows, cols)` of odd positive integers.

- looks:

  Single positive number: the effective number of looks (ENL) of the
  intensity image. Controls the assumed speckle strength
  `cu = 1 / sqrt(looks)`; single-look SAR intensity is `looks = 1`.

- damping:

  Single positive number controlling how quickly the blend moves towards
  the observed value as heterogeneity grows. Larger values preserve more
  detail; `1` is the conventional default.

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

Lopes, A., Touzi, R., & Nezry, E. (1990). Adaptive speckle filters and
scene heterogeneity. *IEEE Transactions on Geoscience and Remote
Sensing*, 28(6), 992-1000.

## See also

[`rf_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee.md),
[`rf_gamma_map()`](https://belian-earth.github.io/rustyfilters/reference/rf_gamma_map.md),
[`rf_frost()`](https://belian-earth.github.io/rustyfilters/reference/rf_frost.md)

## Examples

``` r
set.seed(1)
speckled <- matrix(rexp(400), 20, 20)
rf_enhanced_lee(speckled, window = 7L, looks = 1, damping = 1)
#>            [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#>  [1,] 0.8599573 0.8015771 0.8031746 0.8195226 0.9184840 0.9783807 0.8933462
#>  [2,] 1.0038048 0.9550046 0.9579469 0.9265457 0.9995365 1.0866342 1.0403570
#>  [3,] 0.8564401 0.9108603 0.9978559 0.9469250 0.9584951 1.0350582 1.0022292
#>  [4,] 1.0387718 1.0328991 1.0365385 1.0339432 1.0300705 1.1382968 1.0993726
#>  [5,] 1.1138560 1.1239948 1.1219752 1.1302138 1.1282241 1.2165455 1.1747412
#>  [6,] 1.1921660 1.1192763 1.1192216 1.1145488 1.1309728 1.1955573 1.1633072
#>  [7,] 1.1805387 1.2133108 1.2050642 1.1775398 1.1910713 1.2611564 1.2548055
#>  [8,] 1.2341914 1.2841140 1.2924744 1.2377898 1.2120123 1.2543359 1.2691215
#>  [9,] 1.2113680 1.2133867 1.2475082 1.2244149 1.1924136 1.2426661 1.2138137
#> [10,] 1.3543133 1.3345840 1.2852617 1.2608707 1.2818339 1.3348544 1.2379129
#> [11,] 1.3813145 1.3183887 1.2700810 1.1890727 1.1912340 1.2190849 1.1466535
#> [12,] 1.2245808 1.1772142 1.1314763 1.0446611 1.0673161 1.1550845 1.0873055
#> [13,] 1.2113664 1.1357312 1.1005479 1.0347725 1.0357261 1.1504050 1.0847400
#> [14,] 1.1667527 1.0527896 1.0236368 0.9477388 0.9078025 1.0231075 0.9634901
#> [15,] 1.1032392 0.9784219 0.9189696 0.8315151 0.8398192 0.9640946 0.8925520
#> [16,] 1.0217537 0.9083867 0.8410741 0.6977294 0.7753439 0.8573428 0.8063806
#> [17,] 0.9027411 0.7538838 0.6831732 0.6875463 0.6851451 0.7717176 0.7692921
#> [18,] 0.7787054 0.6924190 0.6646908 0.6680620 0.6617819 0.7424470 0.7236045
#> [19,] 0.8223736 0.7139953 0.6528878 0.6744733 0.6420964 0.7077106 0.6920222
#> [20,] 0.8350031 0.7249361 0.6465357 0.6654011 0.6197817 0.7005440 0.6895813
#>            [,8]      [,9]     [,10]     [,11]     [,12]     [,13]     [,14]
#>  [1,] 0.8776376 0.9195465 1.0309285 1.1020230 1.0841217 1.1443039 1.1338134
#>  [2,] 0.8946832 0.9077877 0.9778178 1.0551681 1.0623356 1.1189994 1.0963282
#>  [3,] 0.8866042 0.8859433 0.8813687 0.9442671 0.9398629 0.9874846 0.9582973
#>  [4,] 0.9898067 0.9950189 0.9886367 0.9966982 0.9813018 0.9939239 0.9676109
#>  [5,] 1.0505380 1.0289826 1.0395111 0.9944731 0.9767558 0.9728153 0.9823340
#>  [6,] 1.0322917 1.0064578 1.0064308 0.9557697 0.8686498 0.8430147 0.8574029
#>  [7,] 1.1183090 1.0710970 1.0475782 1.0277754 0.9755530 0.9348305 0.9063315
#>  [8,] 1.1150827 1.0369882 1.0094658 1.0007878 0.9616399 0.9122278 0.8805136
#>  [9,] 1.1356990 1.1566527 1.1200270 1.1043985 1.0467003 0.9640272 0.9976558
#> [10,] 1.1164310 1.1310693 1.1653649 1.1516575 1.0883953 0.9927405 1.0022404
#> [11,] 1.0323087 1.0459528 1.0801816 1.1046864 1.0005487 1.0843132 0.8526462
#> [12,] 1.0020781 1.0262536 1.0265955 1.1099579 0.9789594 0.9368986 0.9184107
#> [13,] 1.0102184 1.0557458 1.0559623 1.2060789 1.0147206 1.0102263 0.7840446
#> [14,] 0.8920802 0.9658906 0.9912819 1.2197476 1.0153069 0.6812003 0.5985630
#> [15,] 0.8670060 0.9541324 0.9636328 1.1883681 1.0316533 0.8224408 1.2698252
#> [16,] 0.8146965 0.8278320 0.8566908 1.0100077 0.9270141 0.9141381 2.0975275
#> [17,] 0.7619090 0.8055938 0.8166064 0.9780943 0.9624561 1.0649265 1.6137440
#> [18,] 0.7188664 0.7642500 0.7800472 1.1696619 0.9794715 1.0936443 1.0627729
#> [19,] 0.6664656 0.7236520 0.7768958 1.1479601 1.0038201 0.9137040 1.1517731
#> [20,] 0.6524951 0.7186932 0.7841219 0.8985596 0.9569594 0.9298899 0.9961682
#>           [,15]     [,16]     [,17]     [,18]     [,19]     [,20]
#>  [1,] 1.1658892 1.1221718 1.0812518 1.0255817 1.0060977 0.8058343
#>  [2,] 1.1529452 1.1002902 1.0583139 1.0214343 1.0071637 0.8087487
#>  [3,] 1.0267323 0.9979056 1.0336292 1.0250397 1.0356107 0.9091940
#>  [4,] 0.9912029 0.9418048 0.9744684 0.9577462 0.9824742 0.8519432
#>  [5,] 0.9963232 0.9516668 0.9086866 0.9164117 0.9396431 0.8052788
#>  [6,] 0.8482155 0.8392318 0.7735632 0.7805097 0.8543760 0.7592208
#>  [7,] 0.9065221 0.8858548 0.8446163 0.8397568 0.8806271 0.7837155
#>  [8,] 0.8757341 0.8673883 0.7964467 0.7876985 0.8191443 0.7260558
#>  [9,] 0.8405286 0.8474435 0.7471901 0.6452021 0.8393416 0.6591811
#> [10,] 0.9809291 0.8563352 0.7134170 0.6913189 0.6537223 0.6091405
#> [11,] 0.7708922 0.8202714 0.6670550 0.5960858 0.6117545 0.5765631
#> [12,] 0.8491391 0.8735333 0.8056605 0.7553861 0.7732169 0.6896577
#> [13,] 1.0223494 0.6889099 0.8347612 0.7734928 0.8194084 0.6682670
#> [14,] 0.7145054 0.7236141 0.7578271 0.7520952 0.7725858 0.7090628
#> [15,] 0.9742201 1.2708296 1.0701412 0.9047316 0.9436544 0.8943941
#> [16,] 0.8911971 1.2936952 1.1073515 0.9981683 1.0225926 0.9986372
#> [17,] 1.0640444 1.0314787 1.1190460 0.9853701 1.0217646 0.9461283
#> [18,] 1.1792789 1.1781366 1.2623775 1.0684836 1.1039147 1.0046836
#> [19,] 1.1694815 1.0250518 1.3065345 0.9915845 1.0023187 0.9318554
#> [20,] 1.0048420 1.0053763 1.0605579 0.9740414 0.9552689 1.0028764
```
