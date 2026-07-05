#' @title rustyfilters: Fast Focal and Speckle Filters for Matrices and Arrays
#'
#' @description
#' Minimal, high-performance moving-window filters for numeric matrices, 3-D
#' arrays and, optionally, terra SpatRaster objects, powered by Rust via
#' extendr with multi-threading from rayon.
#'
#' @section Speckle filters:
#' - [rf_lee()] --- Lee (1980) minimum mean-square-error filter
#' - [rf_enhanced_lee()] --- enhanced Lee (Lopes et al. 1990)
#' - [rf_lee_sigma()] --- Lee sigma (1983) two-sigma filter
#' - [rf_lee_sigma_improved()] --- improved Lee sigma (Lee et al. 2009)
#' - [rf_frost()] --- Frost (1982) damped exponential filter
#' - [rf_kuan()] --- Kuan (1985) filter
#' - [rf_gamma_map()] --- Gamma maximum a posteriori (Lopes et al. 1990)
#'
#' @section Smoothing:
#' - [rf_mean()] --- moving-window (boxcar) mean
#' - [rf_gaussian()] --- separable Gaussian blur
#' - [rf_median()] --- moving-window median
#' - [rf_bilateral()] --- bilateral edge-preserving smoother
#' - [rf_guided()] --- guided filter (He et al. 2013)
#'
#' @section Focal statistics:
#' - [rf_focal()] --- min, max, range, standard deviation, sum, mode and
#'   friends over a moving window
#'
#' @section Convolution and edges:
#' - [rf_convolve()] --- arbitrary kernel filtering
#' - [rf_sobel()] --- Sobel gradient
#' - [rf_laplacian()] --- discrete Laplacian
#'
#' @section Data and visualisation:
#' - [s1_sthelens] --- a real Sentinel-1 backscatter patch to experiment on
#' - [rf_plot()] --- quick image plots with a percentile stretch
#'
#' @section Configuration:
#' - [rf_set_threads()] / [rf_get_threads()] --- multi-threading control
#'
#' @importFrom rlang %||%
#' @keywords internal
"_PACKAGE"
