
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rustyfilters

<!-- badges: start -->

<!-- badges: end -->

Minimal, blazing-fast moving-window filters for R matrices, 3-D arrays
and terra SpatRaster objects, powered by Rust and rayon. Includes SAR
speckle filters (Lee, enhanced Lee, Lee sigma, Frost, Kuan, Gamma-MAP),
smoothing filters (mean, Gaussian, median) and focal statistics (min,
max, range, sd, sum, mode).

Under construction; a full README with examples and benchmarks will
follow once the API lands.

## Installation

``` r
# install.packages("pak")
pak::pak("belian-earth/rustyfilters")
```
