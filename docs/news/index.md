# Changelog

## rustyfilters 0.0.0.9000

- All filters gained methods for open gdalraster `GDALRaster` datasets
  (gdalraster in Suggests): results are returned as a new `GDALRaster`
  object on an in-memory `/vsimem` GTiff, or written to `filename`.

- Added convolution and gradient filters
  ([`rf_convolve()`](https://belian-earth.github.io/rustyfilters/reference/rf_convolve.md),
  [`rf_sobel()`](https://belian-earth.github.io/rustyfilters/reference/rf_sobel.md),
  [`rf_laplacian()`](https://belian-earth.github.io/rustyfilters/reference/rf_laplacian.md)),
  edge-preserving smoothers
  ([`rf_bilateral()`](https://belian-earth.github.io/rustyfilters/reference/rf_bilateral.md),
  [`rf_guided()`](https://belian-earth.github.io/rustyfilters/reference/rf_guided.md))
  and the improved Lee sigma speckle filter of Lee et al. (2009)
  ([`rf_lee_sigma_improved()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee_sigma_improved.md)).

- Results are now bitwise identical across thread counts for every
  filter: all filters execute one shared parallel code path on a
  dedicated rayon pool, whatever
  [`rf_set_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)
  is set to.

- Initial release: SAR speckle filters
  ([`rf_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee.md),
  [`rf_enhanced_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_enhanced_lee.md),
  [`rf_lee_sigma()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee_sigma.md),
  [`rf_frost()`](https://belian-earth.github.io/rustyfilters/reference/rf_frost.md),
  [`rf_kuan()`](https://belian-earth.github.io/rustyfilters/reference/rf_kuan.md),
  [`rf_gamma_map()`](https://belian-earth.github.io/rustyfilters/reference/rf_gamma_map.md)),
  smoothing filters
  ([`rf_mean()`](https://belian-earth.github.io/rustyfilters/reference/rf_mean.md),
  [`rf_gaussian()`](https://belian-earth.github.io/rustyfilters/reference/rf_gaussian.md),
  [`rf_median()`](https://belian-earth.github.io/rustyfilters/reference/rf_median.md))
  and focal statistics
  ([`rf_focal()`](https://belian-earth.github.io/rustyfilters/reference/rf_focal.md))
  for matrices, 3-D arrays and terra SpatRaster objects, with thread
  control via
  [`rf_set_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md).
