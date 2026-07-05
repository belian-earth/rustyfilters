# rustyfilters 0.0.0.9000

* New `s1_sthelens` dataset: a real Sentinel-1 RTC backscatter patch (VV,
  linear power, 300 x 300 at 10 m) over the crater of Mount St Helens, for
  realistic speckle-filter examples. Contains modified Copernicus Sentinel
  data (2024).

* New `rf_plot()` helper: correctly oriented `image()` plots with a
  percentile stretch and an Inferno palette; all documentation examples now
  plot instead of printing matrices.

* `GDALRaster` methods stream rasters larger than
  `options(rustyfilters.block_memory)` (default 2 GiB) through full-width
  row bands with a filter-sized halo, writing to a GeoTIFF tempfile (or
  `filename`). Interior band seams are exact; `by_block` and `block_rows`
  control the behaviour. `rf_lee_sigma_improved()` computes its global 98th
  percentile exactly when streaming, so tiled results match the whole-image
  run bit for bit.

* All filters gained methods for open gdalraster `GDALRaster` datasets
  (gdalraster in Suggests): results are returned as a new `GDALRaster`
  object on an in-memory `/vsimem` GTiff, or written to `filename`.

* Added convolution and gradient filters (`rf_convolve()`, `rf_sobel()`,
  `rf_laplacian()`), edge-preserving smoothers (`rf_bilateral()`,
  `rf_guided()`) and the improved Lee sigma speckle filter of Lee et
  al. (2009) (`rf_lee_sigma_improved()`).

* Results are now bitwise identical across thread counts for every filter:
  all filters execute one shared parallel code path on a dedicated rayon
  pool, whatever `rf_set_threads()` is set to.

* Initial release: SAR speckle filters (`rf_lee()`, `rf_enhanced_lee()`,
  `rf_lee_sigma()`, `rf_frost()`, `rf_kuan()`, `rf_gamma_map()`), smoothing
  filters (`rf_mean()`, `rf_gaussian()`, `rf_median()`) and focal statistics
  (`rf_focal()`) for matrices, 3-D arrays and terra SpatRaster objects, with
  thread control via `rf_set_threads()`.
