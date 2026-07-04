# Improved Lee sigma speckle filter

The improved sigma filter of Lee et al. (2009), as popularised by ESA
SNAP. It addresses the classic sigma filter's dark bias by (1)
estimating an a priori mean with a small MMSE filter and using published
sigma range bounds around it, (2) preserving point targets, detected as
clusters of pixels above the scene's 98th percentile, and (3) filtering
the in-range pixels with an MMSE weight based on a revised speckle
variation coefficient.

## Usage

``` r
# S3 method for class 'Rcpp_GDALRaster'
rf_lee_sigma_improved(x, ...)

rf_lee_sigma_improved(x, ...)

# S3 method for class 'matrix'
rf_lee_sigma_improved(
  x,
  window = 7L,
  looks = 1,
  sigma = 0.9,
  target_window = 3L,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# S3 method for class 'array'
rf_lee_sigma_improved(
  x,
  window = 7L,
  looks = 1,
  sigma = 0.9,
  target_window = 3L,
  edge = c("shrink", "reflect", "nearest", "constant"),
  edge_value = 0,
  na_policy = c("omit", "propagate"),
  ...
)

# Default S3 method
rf_lee_sigma_improved(x, ...)

# S3 method for class 'SpatRaster'
rf_lee_sigma_improved(x, ...)
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

  Effective number of looks. Must be 1, 2, 3 or 4; the published sigma
  range tables only cover these.

- sigma:

  Sigma confidence level: one of 0.5, 0.6, 0.7, 0.8 or 0.9 (the default,
  covering the widest speckle range).

- target_window:

  Single odd integer (typically 3 or 5): the window used for the a
  priori mean estimate and for point-target detection.

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

Point-target marking is order-independent here: a pixel is preserved
when it exceeds the layer's 98th percentile and lies within
`target_window` of a cluster of more than five such pixels. SNAP marks
clusters during a sequential scan, which can differ at cluster fringes.
The 98th percentile is computed per layer over valid cells.

## References

Lee, J.-S., Wen, J.-H., Ainsworth, T. L., Chen, K.-S., & Chen, A. J.
(2009). Improved sigma filter for speckle filtering of SAR imagery.
*IEEE Transactions on Geoscience and Remote Sensing*, 47(1), 202-213.

## See also

[`rf_lee_sigma()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee_sigma.md)
for the classic 1983 filter,
[`rf_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee.md),
[`rf_enhanced_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_enhanced_lee.md)

## Examples

``` r
set.seed(1)
speckled <- matrix(rexp(400), 20, 20)
rf_lee_sigma_improved(speckled, window = 7L, looks = 1, sigma = 0.9)
#>            [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#>  [1,] 0.9113095 0.8823671 0.8355186 0.8476805 0.9503070 1.0124222 0.9242382
#>  [2,] 0.7819000 0.8414312 0.8516547 0.8598993 0.9397329 1.0349961 0.9843804
#>  [3,] 0.6325447 0.8034424 0.9385242 0.9163240 0.8850933 0.9887211 0.9533667
#>  [4,] 0.4757308 0.7014977 1.0129850 1.1697321 1.0320913 1.1057439 1.0642812
#>  [5,] 0.7137988 0.8039585 1.0055925 1.3818162 1.2339606 1.1890958 1.1445651
#>  [6,] 1.2446846 0.9762653 1.0782415 1.1943698 1.2518531 1.2698193 1.1323854
#>  [7,] 1.2607831 1.3200345 1.2919981 1.1962858 1.3187810 1.3396966 1.3329315
#>  [8,] 1.2765096 1.3201338 1.3529603 1.3416628 1.3135940 1.3053425 1.2444083
#>  [9,] 1.3008307 1.2839947 1.3074659 1.3001851 1.2660967 1.2666569 1.2372034
#> [10,] 1.4030927 1.3727403 1.3476642 1.3395029 1.3618332 1.3618898 1.2629287
#> [11,] 1.1266201 1.2631031 1.3317245 1.2630224 1.2653247 1.2682319 1.1697682
#> [12,] 1.1471946 1.2107420 1.2133851 1.1542613 1.1560701 1.1783748 1.1091837
#> [13,] 1.3011129 1.2018566 1.1819625 1.1215017 1.1005240 1.0397124 1.1290012
#> [14,] 1.3583610 1.1310282 1.0116405 0.9308620 0.8688221 0.8203201 1.0025916
#> [15,] 1.1714289 0.8991769 0.7637985 0.5835020 0.6264797 0.7514717 0.9473137
#> [16,] 0.9295777 0.7485785 0.6000906 0.5515795 0.6266990 0.7326394 0.8746637
#> [17,] 0.8000293 0.7013706 0.5929977 0.5685819 0.5138987 0.4796222 0.6904785
#> [18,] 0.8102989 0.6414739 0.6211050 0.6170575 0.5527548 0.4751970 0.6728230
#> [19,] 0.7204330 0.6287766 0.7094348 0.7168963 0.5780660 0.5378255 0.8019138
#> [20,] 0.7067204 0.6178673 0.6865214 1.0272456 0.5719102 0.6042785 0.6544517
#>            [,8]      [,9]     [,10]     [,11]     [,12]     [,13]     [,14]
#>  [1,] 0.9079478 0.9195465 1.0309285 1.0077546 1.1226070 1.2276872 1.2163897
#>  [2,] 0.9737118 0.9598999 1.0059815 0.9771283 1.0922567 1.1831617 1.1591165
#>  [3,] 0.9508267 0.9278436 0.9236752 0.9330931 1.0077788 1.0590637 1.0771563
#>  [4,] 1.0509868 1.0353202 1.0715860 1.0370844 1.0870808 1.1011373 0.9700858
#>  [5,] 1.0464677 1.0499979 1.0822462 1.0347646 1.0382643 0.9734686 0.8401189
#>  [6,] 1.0745273 1.0270038 1.0906177 1.0143857 0.9231078 0.8958009 0.5649677
#>  [7,] 1.1875330 1.0929897 1.1354226 1.0481410 1.0363153 0.9345781 0.9231896
#>  [8,] 1.1601637 1.0581703 1.0550627 0.9690856 1.0191184 1.0047502 0.8537452
#>  [9,] 1.1574613 1.1788516 1.1128497 1.0008244 0.9917211 0.9532840 1.0453405
#> [10,] 1.1612098 1.1764709 1.1157687 1.0298689 1.0579127 1.0248603 1.1065197
#> [11,] 1.0735078 1.0839356 1.1257344 1.1422390 1.0717371 1.3563342 0.8696742
#> [12,] 1.0221807 0.9715586 0.9703926 1.0825155 1.2160744 0.9784350 0.9217900
#> [13,] 1.0304906 1.0769665 1.0032857 1.0705413 0.9494421 0.9234635 0.6768728
#> [14,] 0.9098912 1.0339662 0.9358530 0.9217059 0.7955241 0.6380508 0.6577370
#> [15,] 0.9370685 1.0378431 0.9070274 0.9431505 0.7271194 0.8442093 1.1779061
#> [16,] 0.8643249 0.8956386 0.9102823 0.8998280 0.8598420 1.0352617 2.1547442
#> [17,] 0.8641259 0.8724404 0.9049682 0.9661470 0.8649604 1.1854829 1.6861553
#> [18,] 0.7888975 0.8390584 0.8570676 0.9094673 0.8656431 1.0084234 1.1055597
#> [19,] 0.6243678 0.5756799 0.8709014 0.8750160 0.8643434 0.7928470 0.8830217
#> [20,] 0.5646218 0.4593008 0.8412934 0.9150036 0.8778902 0.8804079 0.9546758
#>           [,15]     [,16]     [,17]     [,18]     [,19]     [,20]
#>  [1,] 1.2474010 1.1250249 1.1119458 0.8421151 1.0040597 1.0052644
#>  [2,] 1.1855312 1.1313276 1.0479463 0.8510043 0.9532200 0.8254314
#>  [3,] 1.1264389 1.0714215 0.9638884 0.7645266 0.9424029 0.7363176
#>  [4,] 0.9957748 0.9419858 0.9457769 0.7508953 0.7845855 0.6518470
#>  [5,] 0.8564396 0.8043405 0.8771269 0.8532259 0.9411779 0.7259205
#>  [6,] 0.9013408 0.8502496 0.7656579 0.8848955 0.8826529 1.1429139
#>  [7,] 0.6292269 0.9413133 0.8760935 0.8971577 0.9265710 0.7320103
#>  [8,] 0.9208903 0.9676478 0.9634848 0.8456007 0.6742883 0.6455763
#>  [9,] 0.9713383 1.0378609 0.8446023 0.6444761 0.8840715 0.6068241
#> [10,] 1.1669174 1.0570429 0.6938321 0.7399312 0.5080930 0.5749530
#> [11,] 0.8575267 0.7691457 0.6452497 0.4268026 0.4231023 0.4398833
#> [12,] 0.6672092 0.6690271 0.5898842 0.6861327 0.5261527 0.7396563
#> [13,] 0.6951016 0.7230172 0.7410682 0.8499279 0.8684025 0.6672543
#> [14,] 0.7113755 0.7469317 0.8992090 0.9154625 0.9546935 0.8548255
#> [15,] 1.0748089 1.0674017 1.1152281 0.9976033 1.0953504 0.9984352
#> [16,] 1.0544819 1.1165598 1.1759850 0.9875493 1.1874004 1.1583651
#> [17,] 1.2064468 1.1025323 1.1765398 1.0245657 1.1855908 1.1433721
#> [18,] 1.1512021 1.1726953 1.1826909 1.1960901 1.2664198 1.1960940
#> [19,] 1.0035052 1.1193281 1.1958545 1.1368080 1.1841639 1.0897356
#> [20,] 1.0012849 1.0019082 1.0670156 1.0586293 1.0565682 1.0029915
```
