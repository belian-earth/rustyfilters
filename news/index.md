# Changelog

## rustyfilters 0.0.0.9000

- New `s1_sthelens` dataset: a real Sentinel-1 RTC backscatter patch
  (VV, linear power, 300 x 300 at 10 m) over the crater of Mount St
  Helens, for realistic speckle-filter examples. Contains modified
  Copernicus Sentinel data (2024).

- New
  [`rf_plot()`](https://belian-earth.github.io/rustyfilters/reference/rf_plot.md)
  helper: correctly oriented
  [`image()`](https://rdrr.io/r/graphics/image.html) plots with a
  percentile stretch and an Inferno palette; all documentation examples
  now plot instead of printing matrices.

- `GDALRaster` methods stream rasters larger than
  `options(rustyfilters.block_memory)` (default 2 GiB) through
  full-width row bands with a filter-sized halo, writing to a GeoTIFF
  tempfile (or `filename`). Interior band seams are exact; `by_block`
  and `block_rows` control the behaviour.
  [`rf_lee_sigma_improved()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee_sigma_improved.md)
  computes its global 98th percentile exactly when streaming, so tiled
  results match the whole-image run bit for bit.

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
