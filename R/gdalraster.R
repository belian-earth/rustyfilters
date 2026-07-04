# GDALRaster methods for every filter. gdalraster is Suggests-only: the
# methods are registered statically but only touch gdalraster inside
# filter_gdalraster(), which guards with check_installed(). The dataset is
# read into memory as a (rows, cols, bands) array (gdalraster masks nodata
# to NA on read), filtered layer by layer, and written to a new Float64
# dataset with the source's geometry: by default an in-memory /vsimem GTiff,
# or `filename` if given. The returned GDALRaster object is open in update
# mode; NA cells are stored as NaN with the band nodata set accordingly.

filter_gdalraster <- function(x, fun, ..., filename = NULL,
                              call = rlang::caller_env()) {
  rlang::check_installed("gdalraster", reason = "to filter GDALRaster objects.")
  src_file <- x$getFilename()
  if (!nzchar(src_file)) {
    cli::cli_abort(
      "{.arg x} must be backed by a file or {.path /vsimem} dataset.",
      call = call
    )
  }
  dm <- x$dim() # xsize (cols), ysize (rows), nbands
  v <- as.numeric(gdalraster::read_ds(x))
  a <- aperm(array(v, dim = dm), c(2L, 1L, 3L)) # rows, cols, bands
  out <- fun(a, ...)
  filename <- filename %||%
    paste0("/vsimem/", basename(tempfile("rustyfilters_")), ".tif")
  gdalraster::rasterFromRaster(src_file, filename,
    fmt = "GTiff", dtName = "Float64"
  )
  dst <- methods::new(gdalraster::GDALRaster, filename, read_only = FALSE)
  for (b in seq_len(dm[[3L]])) {
    dst$setNoDataValue(b, NaN)
    dst$write(
      band = b, xoff = 0L, yoff = 0L,
      xsize = dm[[1L]], ysize = dm[[2L]],
      rasterData = as.numeric(t(out[, , b]))
    )
  }
  dst
}

#' @rdname rf_lee
#' @export
rf_lee.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_lee, ...)
}

#' @rdname rf_enhanced_lee
#' @export
rf_enhanced_lee.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_enhanced_lee, ...)
}

#' @rdname rf_lee_sigma
#' @export
rf_lee_sigma.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_lee_sigma, ...)
}

#' @rdname rf_lee_sigma_improved
#' @export
rf_lee_sigma_improved.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_lee_sigma_improved, ...)
}

#' @rdname rf_frost
#' @export
rf_frost.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_frost, ...)
}

#' @rdname rf_kuan
#' @export
rf_kuan.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_kuan, ...)
}

#' @rdname rf_gamma_map
#' @export
rf_gamma_map.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_gamma_map, ...)
}

#' @rdname rf_mean
#' @export
rf_mean.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_mean, ...)
}

#' @rdname rf_gaussian
#' @export
rf_gaussian.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_gaussian, ...)
}

#' @rdname rf_median
#' @export
rf_median.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_median, ...)
}

#' @rdname rf_focal
#' @export
rf_focal.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_focal, ...)
}

#' @rdname rf_convolve
#' @export
rf_convolve.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_convolve, ...)
}

#' @rdname rf_sobel
#' @export
rf_sobel.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_sobel, ...)
}

#' @rdname rf_laplacian
#' @export
rf_laplacian.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_laplacian, ...)
}

#' @rdname rf_bilateral
#' @export
rf_bilateral.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_bilateral, ...)
}

#' @rdname rf_guided
#' @export
rf_guided.Rcpp_GDALRaster <- function(x, guide = NULL, ...) {
  if (inherits(guide, "Rcpp_GDALRaster")) {
    gd <- guide$dim()
    guide <- aperm(
      array(as.numeric(gdalraster::read_ds(guide)), dim = gd),
      c(2L, 1L, 3L)
    )
  }
  filter_gdalraster(x, rf_guided, guide = guide, ...)
}
