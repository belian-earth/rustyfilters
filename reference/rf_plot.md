# Plot a matrix or array layer as an image

A small convenience around
[`graphics::image()`](https://rdrr.io/r/graphics/image.html) used
throughout the package examples: correct raster orientation (row 1 at
the top), square pixels, no axes, a percentile stretch that keeps skewed
data (like SAR backscatter) readable, and an Inferno palette from
[`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html).

## Usage

``` r
rf_plot(x, layer = 1L, palette = "Inferno", stretch = c(0.02, 0.98), ...)
```

## Arguments

- x:

  A numeric matrix or 3-D array.

- layer:

  Which layer of a 3-D array to plot.

- palette:

  A palette name understood by
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html).

- stretch:

  Two probabilities: values outside these quantiles are clipped before
  plotting. Use `c(0, 1)` for no stretch.

- ...:

  Passed on to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html) (e.g.
  `main`).

## Value

`x`, invisibly.

## Examples

``` r
rf_plot(volcano, main = "Maungawhau")

rf_plot(s1_sthelens, main = "Mount St Helens, Sentinel-1 VV")
```
