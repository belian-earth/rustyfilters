# Frost speckle filter

The Frost (1982) filter convolves each window with a damped exponential
kernel `exp(-B * d)` where `d` is the distance from the centre and
`B = damping * ci^2` adapts to the local variation coefficient:
homogeneous windows average broadly, heterogeneous windows concentrate
weight on the centre pixel.

## Usage

``` r
rf_frost(x, ...)

# S3 method for class 'matrix'
rf_frost(
  x,
  window = 7L,
  damping = 2,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_frost(
  x,
  window = 7L,
  damping = 2,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_frost(x, ...)

# S3 method for class 'SpatRaster'
rf_frost(x, ...)
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

- damping:

  Single positive number scaling the kernel decay. Larger values
  preserve more edges and smooth less; `2` is a common default.

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

Frost, V. S., Stiles, J. A., Shanmugan, K. S., & Holtzman, J. C. (1982).
A model for radar images and its application to adaptive digital
filtering of multiplicative noise. *IEEE Transactions on Pattern
Analysis and Machine Intelligence*, 4(2), 157-166.

## See also

[`rf_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee.md),
[`rf_enhanced_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_enhanced_lee.md),
[`rf_gamma_map()`](https://belian-earth.github.io/rustyfilters/reference/rf_gamma_map.md)

## Examples

``` r
set.seed(1)
speckled <- matrix(rexp(400), 20, 20)
rf_frost(speckled, window = 7L, damping = 2)
#>            [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#>  [1,] 0.9880259 1.2540708 0.9809786 0.7236168 0.7884373 0.9438676 0.8376948
#>  [2,] 1.0339086 0.8192408 0.9825914 0.9234711 0.8624601 0.8188529 0.9061336
#>  [3,] 0.2740817 0.4878429 1.0525315 0.6597502 0.7909160 0.7845516 1.2221052
#>  [4,] 0.2864767 0.6101091 1.1158252 1.4114014 0.7975413 0.8438963 0.9699484
#>  [5,] 0.6007021 0.4876356 0.9717877 2.5115749 1.6969940 1.1818869 0.7982455
#>  [6,] 2.0856072 0.6060940 0.7987871 1.1621809 1.3147432 1.5073272 0.9586281
#>  [7,] 1.2811402 1.0878945 1.2653665 1.7171074 1.4379483 1.2427823 1.5460661
#>  [8,] 1.0642984 1.8192921 1.2919631 1.2863611 1.3242699 1.4288866 1.2778052
#>  [9,] 1.0914367 1.2412814 1.1828053 1.1826455 1.2497127 1.1621572 1.0971832
#> [10,] 0.9525519 1.1744157 1.2966894 1.2410237 1.3273638 1.1696740 1.1376871
#> [11,] 1.1160844 1.1851062 1.1879812 1.4233505 1.2379557 1.1792752 1.0054318
#> [12,] 1.0074012 0.9918754 1.5644228 1.6872089 1.1241875 1.1421515 1.0238004
#> [13,] 1.3771289 1.1061432 1.8241563 1.7533360 1.0412884 0.8247830 0.9069414
#> [14,] 2.6059324 1.3227756 0.9647125 0.7370035 0.6044644 0.6254360 0.8571584
#> [15,] 1.2723718 0.6644144 0.6439800 0.4580936 0.6174918 0.8556088 0.7502719
#> [16,] 1.0662969 0.9000457 0.7845836 0.2308162 0.3976539 0.6660470 0.8187508
#> [17,] 1.3891011 0.5356510 0.3701538 0.4297904 0.4449407 0.4646003 0.5156119
#> [18,] 0.7517660 0.6978047 0.6204622 0.8612714 0.4159981 0.3625227 0.4770625
#> [19,] 0.5917537 0.7021354 0.8993028 0.8725033 0.4179717 0.4328351 1.1116611
#> [20,] 0.5789509 0.5158839 0.9073045 1.8415816 0.5440214 0.4252365 0.7211763
#>            [,8]      [,9]     [,10]     [,11]     [,12]     [,13]     [,14]
#>  [1,] 0.8916532 0.9482358 0.8869797 0.8078962 0.9880426 1.3753169 1.6454808
#>  [2,] 1.0143560 0.9623562 0.8657161 0.8853971 0.9188924 1.2303554 1.8141817
#>  [3,] 1.2236690 1.0740785 0.8684343 0.8556187 0.9707712 1.1948744 0.8821342
#>  [4,] 1.0825589 1.1077404 0.9584848 0.9030160 1.0301336 1.3594060 0.8729377
#>  [5,] 0.8786558 1.0415964 0.9922697 0.8064547 1.0265976 0.8002634 0.7294053
#>  [6,] 0.9530759 0.8999442 0.8912026 0.8887542 0.7521691 0.5659651 0.4825963
#>  [7,] 1.1189343 1.3099558 1.0249223 1.1464772 1.0980294 0.7460008 0.7357983
#>  [8,] 1.1114427 1.3022860 0.8675001 0.8191336 0.9125042 1.4280050 0.7544065
#>  [9,] 1.2345502 1.0552391 0.7754274 0.7308644 0.7301581 0.9116626 1.6018871
#> [10,] 1.2798879 1.5553718 1.4374582 1.0281634 1.0558867 0.9494431 1.1296040
#> [11,] 1.0193426 1.1632132 1.1650272 0.8925147 1.1571495 2.3908254 0.8436548
#> [12,] 0.9523847 1.2219499 0.9792394 1.2568956 2.4059375 1.1969414 1.2179675
#> [13,] 1.1564287 1.0592579 0.7676481 0.9512470 0.6801936 0.8237763 0.4846224
#> [14,] 1.3611647 1.4897918 1.0190924 0.4784936 0.4987309 0.1238212 0.1382790
#> [15,] 1.2848106 1.2727845 0.8098734 0.8204439 0.5881128 0.4418961 1.7742335
#> [16,] 1.0198535 1.1955123 0.8883687 0.4370911 0.5676136 1.0749932 4.7937563
#> [17,] 1.2171276 1.2306537 0.9126087 0.7579561 0.8981122 1.4555344 3.2571112
#> [18,] 0.7494434 0.8102218 0.8230107 1.6389047 0.6834600 1.0542639 0.7496948
#> [19,] 0.4595623 0.5977462 0.7410967 1.3267393 0.7037596 0.3143043 1.1010928
#> [20,] 0.2873266 0.5740572 0.4589761 0.4588997 0.7778375 0.5960911 0.8628005
#>           [,15]     [,16]     [,17]     [,18]     [,19]     [,20]
#>  [1,] 1.3577040 1.0910949 1.0314376 0.7157003 1.0718360 1.9370019
#>  [2,] 1.7957897 1.8749216 0.9854354 1.1192094 0.7928721 1.2984168
#>  [3,] 1.0742814 1.4700482 0.7225724 0.6460002 0.7073400 0.6077559
#>  [4,] 0.7921256 0.9970890 0.6711854 0.8014804 0.6042367 0.5939297
#>  [5,] 0.8348614 1.1673795 0.8383761 1.0465206 0.9051172 0.6872931
#>  [6,] 0.5607378 0.6142186 0.5870348 1.0923801 1.0160142 1.9601954
#>  [7,] 0.5525891 0.8455343 0.8290430 0.7493991 0.7781489 0.6022463
#>  [8,] 0.6555947 0.7364303 1.4704661 0.6168898 0.7397608 0.3434247
#>  [9,] 0.7805644 1.4895254 0.8364545 0.4221737 1.4124019 0.4748717
#> [10,] 1.5087729 1.5163438 0.5050319 0.8003476 0.5230080 0.5116513
#> [11,] 0.4920223 0.4675202 0.4450957 0.2430465 0.3411028 0.1863002
#> [12,] 0.3784503 0.5229601 0.7819637 0.3790086 0.4201308 0.5221422
#> [13,] 0.8733769 0.1666771 0.6604576 0.6310211 1.2783102 0.4710397
#> [14,] 0.3509677 0.3961045 0.5617775 1.1021153 0.7781479 0.4899553
#> [15,] 0.9010501 2.0402420 1.1188946 1.1782242 1.1138217 1.0408595
#> [16,] 1.0533992 2.3866903 0.6217694 0.9554349 0.8530843 1.0922257
#> [17,] 1.3339759 0.8869342 0.5870531 0.8528170 0.9028187 1.3755508
#> [18,] 1.2563731 0.8611058 1.2275995 1.3655977 1.1158489 1.0025817
#> [19,] 1.1228962 0.6285067 1.6974829 0.8976816 1.0655047 1.1476705
#> [20,] 0.8078913 1.0835159 0.7300031 1.0148241 0.5831073 0.4085167
```
