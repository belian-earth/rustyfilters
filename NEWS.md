# rustyfilters 0.0.0.9000

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
