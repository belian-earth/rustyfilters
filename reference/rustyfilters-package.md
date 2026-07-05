# rustyfilters: Fast Focal and Speckle Filters for Matrices and Arrays

Minimal, high-performance moving-window filters for numeric matrices,
3-D arrays and, optionally, terra SpatRaster objects, powered by Rust
via extendr with multi-threading from rayon.

## Speckle filters

- [`rf_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee.md)
  — Lee (1980) minimum mean-square-error filter

- [`rf_enhanced_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_enhanced_lee.md)
  — enhanced Lee (Lopes et al. 1990)

- [`rf_lee_sigma()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee_sigma.md)
  — Lee sigma (1983) two-sigma filter

- [`rf_lee_sigma_improved()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee_sigma_improved.md)
  — improved Lee sigma (Lee et al. 2009)

- [`rf_frost()`](https://belian-earth.github.io/rustyfilters/reference/rf_frost.md)
  — Frost (1982) damped exponential filter

- [`rf_kuan()`](https://belian-earth.github.io/rustyfilters/reference/rf_kuan.md)
  — Kuan (1985) filter

- [`rf_gamma_map()`](https://belian-earth.github.io/rustyfilters/reference/rf_gamma_map.md)
  — Gamma maximum a posteriori (Lopes et al. 1990)

## Smoothing

- [`rf_mean()`](https://belian-earth.github.io/rustyfilters/reference/rf_mean.md)
  — moving-window (boxcar) mean

- [`rf_gaussian()`](https://belian-earth.github.io/rustyfilters/reference/rf_gaussian.md)
  — separable Gaussian blur

- [`rf_median()`](https://belian-earth.github.io/rustyfilters/reference/rf_median.md)
  — moving-window median

- [`rf_bilateral()`](https://belian-earth.github.io/rustyfilters/reference/rf_bilateral.md)
  — bilateral edge-preserving smoother

- [`rf_guided()`](https://belian-earth.github.io/rustyfilters/reference/rf_guided.md)
  — guided filter (He et al. 2013)

## Focal statistics

- [`rf_focal()`](https://belian-earth.github.io/rustyfilters/reference/rf_focal.md)
  — min, max, range, standard deviation, sum, mode and friends over a
  moving window

## Convolution and edges

- [`rf_convolve()`](https://belian-earth.github.io/rustyfilters/reference/rf_convolve.md)
  — arbitrary kernel filtering

- [`rf_sobel()`](https://belian-earth.github.io/rustyfilters/reference/rf_sobel.md)
  — Sobel gradient

- [`rf_laplacian()`](https://belian-earth.github.io/rustyfilters/reference/rf_laplacian.md)
  — discrete Laplacian

## Data and visualisation

- [s1_sthelens](https://belian-earth.github.io/rustyfilters/reference/s1_sthelens.md)
  — a real Sentinel-1 backscatter patch to experiment on

- [`rf_plot()`](https://belian-earth.github.io/rustyfilters/reference/rf_plot.md)
  — quick image plots with a percentile stretch

## Configuration

- [`rf_set_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)
  /
  [`rf_get_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)
  — multi-threading control

## See also

Useful links:

- <https://github.com/belian-earth/rustyfilters>

- <https://belian-earth.github.io/rustyfilters/>

- Report bugs at <https://github.com/belian-earth/rustyfilters/issues>

## Author

**Maintainer**: Hugh Graham <hugh@belian.earth>

Authors:

- Hugh Graham <hugh@belian.earth>
