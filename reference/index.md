# Package index

## Speckle filters

Adaptive local-statistics filters for SAR intensity data.

- [`rf_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee.md)
  : Lee speckle filter
- [`rf_enhanced_lee()`](https://belian-earth.github.io/rustyfilters/reference/rf_enhanced_lee.md)
  : Enhanced Lee speckle filter
- [`rf_lee_sigma()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee_sigma.md)
  : Lee sigma speckle filter
- [`rf_lee_sigma_improved()`](https://belian-earth.github.io/rustyfilters/reference/rf_lee_sigma_improved.md)
  : Improved Lee sigma speckle filter
- [`rf_frost()`](https://belian-earth.github.io/rustyfilters/reference/rf_frost.md)
  : Frost speckle filter
- [`rf_kuan()`](https://belian-earth.github.io/rustyfilters/reference/rf_kuan.md)
  : Kuan speckle filter
- [`rf_gamma_map()`](https://belian-earth.github.io/rustyfilters/reference/rf_gamma_map.md)
  : Gamma-MAP speckle filter

## Smoothing

Classic and edge-preserving moving-window smoothers.

- [`rf_mean()`](https://belian-earth.github.io/rustyfilters/reference/rf_mean.md)
  : Moving-window mean (boxcar) filter
- [`rf_gaussian()`](https://belian-earth.github.io/rustyfilters/reference/rf_gaussian.md)
  : Gaussian blur
- [`rf_median()`](https://belian-earth.github.io/rustyfilters/reference/rf_median.md)
  : Moving-window median filter
- [`rf_bilateral()`](https://belian-earth.github.io/rustyfilters/reference/rf_bilateral.md)
  : Bilateral filter
- [`rf_guided()`](https://belian-earth.github.io/rustyfilters/reference/rf_guided.md)
  : Guided filter

## Focal statistics

Window summaries beyond the mean.

- [`rf_focal()`](https://belian-earth.github.io/rustyfilters/reference/rf_focal.md)
  : Focal statistics over a moving window

## Convolution & edges

Arbitrary kernels and gradient filters.

- [`rf_convolve()`](https://belian-earth.github.io/rustyfilters/reference/rf_convolve.md)
  : Filter with an arbitrary kernel
- [`rf_sobel()`](https://belian-earth.github.io/rustyfilters/reference/rf_sobel.md)
  : Sobel edge detection
- [`rf_laplacian()`](https://belian-earth.github.io/rustyfilters/reference/rf_laplacian.md)
  : Laplacian filter

## Data & visualisation

Example data and quick image plotting.

- [`s1_sthelens`](https://belian-earth.github.io/rustyfilters/reference/s1_sthelens.md)
  : Sentinel-1 backscatter over Mount St Helens
- [`rf_plot()`](https://belian-earth.github.io/rustyfilters/reference/rf_plot.md)
  : Plot a matrix or array layer as an image

## Configuration

Multi-threading control.

- [`rf_set_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)
  [`rf_get_threads()`](https://belian-earth.github.io/rustyfilters/reference/rf_threads.md)
  : Control the number of threads used by rustyfilters
