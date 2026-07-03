# SpatRaster methods for every filter. terra is Suggests-only: the methods
# are registered statically but only touch terra inside filter_spatraster(),
# which guards with check_installed(). The raster is materialised in memory
# as a (rows, cols, layers) array, filtered layer by layer, and rebuilt with
# the source geometry.

filter_spatraster <- function(x, fun, ..., call = rlang::caller_env()) {
  rlang::check_installed("terra", reason = "to filter SpatRaster objects.")
  a <- terra::as.array(x)
  out <- fun(a, ...)
  r <- terra::rast(out, crs = terra::crs(x), extent = terra::ext(x))
  names(r) <- names(x)
  r
}

#' @rdname rf_lee
#' @export
rf_lee.SpatRaster <- function(x, ...) {
  filter_spatraster(x, rf_lee, ...)
}

#' @rdname rf_enhanced_lee
#' @export
rf_enhanced_lee.SpatRaster <- function(x, ...) {
  filter_spatraster(x, rf_enhanced_lee, ...)
}

#' @rdname rf_lee_sigma
#' @export
rf_lee_sigma.SpatRaster <- function(x, ...) {
  filter_spatraster(x, rf_lee_sigma, ...)
}

#' @rdname rf_frost
#' @export
rf_frost.SpatRaster <- function(x, ...) {
  filter_spatraster(x, rf_frost, ...)
}

#' @rdname rf_kuan
#' @export
rf_kuan.SpatRaster <- function(x, ...) {
  filter_spatraster(x, rf_kuan, ...)
}

#' @rdname rf_gamma_map
#' @export
rf_gamma_map.SpatRaster <- function(x, ...) {
  filter_spatraster(x, rf_gamma_map, ...)
}

#' @rdname rf_mean
#' @export
rf_mean.SpatRaster <- function(x, ...) {
  filter_spatraster(x, rf_mean, ...)
}

#' @rdname rf_gaussian
#' @export
rf_gaussian.SpatRaster <- function(x, ...) {
  filter_spatraster(x, rf_gaussian, ...)
}

#' @rdname rf_median
#' @export
rf_median.SpatRaster <- function(x, ...) {
  filter_spatraster(x, rf_median, ...)
}

#' @rdname rf_focal
#' @export
rf_focal.SpatRaster <- function(x, ...) {
  filter_spatraster(x, rf_focal, ...)
}
