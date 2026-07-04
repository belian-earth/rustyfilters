# GDALRaster methods for every filter. gdalraster is Suggests-only: the
# methods are registered statically but only touch gdalraster inside
# filter_gdalraster(), which guards with check_installed().
#
# Small datasets are read whole, filtered, and written to a new Float64
# dataset with the source's geometry (an in-memory /vsimem GTiff by default,
# or `filename`). Datasets whose decoded size exceeds
# `options(rustyfilters.block_memory)` (default 2 GiB) stream through the
# halo'd row-band engine in R/blocks.R instead, writing to `filename` or a
# GeoTIFF tempfile; pass `by_block = TRUE/FALSE` to force either path and
# `block_rows` to override the band height. The returned GDALRaster object
# is open in update mode; NA cells are stored as NaN with the band nodata
# set accordingly.

filter_gdalraster <- function(x, fun, ..., radius, filename = NULL,
                              by_block = NULL, block_rows = NULL,
                              block_prep = NULL, global_prep = NULL,
                              call = rlang::caller_env()) {
  rlang::check_installed("gdalraster", reason = "to filter GDALRaster objects.")
  src_file <- x$getFilename()
  if (!nzchar(src_file)) {
    cli::cli_abort(
      "{.arg x} must be backed by a file or {.path /vsimem} dataset.",
      call = call
    )
  }
  if (!is.null(by_block) && !rlang::is_bool(by_block)) {
    cli::cli_abort("{.arg by_block} must be `TRUE`, `FALSE` or `NULL`.",
      call = call
    )
  }
  dm <- x$dim()
  budget <- getOption("rustyfilters.block_memory", 2^31)
  stream <- by_block %||% (prod(dm) * 8 > budget)
  if (!stream) {
    v <- as.numeric(gdalraster::read_ds(x))
    a <- aperm(array(v, dim = dm), c(2L, 1L, 3L))
    # The eager path is one full-height block.
    extra <- if (is.null(block_prep)) list() else block_prep(0L, dm[[2L]])
    out <- rlang::exec(fun, a, ..., !!!extra)
    filename <- filename %||%
      paste0("/vsimem/", basename(tempfile("rustyfilters_")), ".tif")
    dst <- gdal_create_like(src_file, filename, dm[[3L]])
    gdal_write_rows(dst, out, 0L)
    return(dst)
  }
  # Streaming: like terra, default to a GeoTIFF tempfile on disk (a /vsimem
  # destination would put the full result back into memory).
  filename <- filename %||% tempfile("rustyfilters_", fileext = ".tif")
  block_rows <- block_rows %||% gdal_block_rows(dm, budget)
  block_rows <- vctrs::vec_cast(block_rows, integer(), call = call)
  if (length(block_rows) != 1L || is.na(block_rows) || block_rows < 1L) {
    cli::cli_abort("{.arg block_rows} must be a single positive integer.",
      call = call
    )
  }
  extra <- if (is.null(global_prep)) list() else global_prep(block_rows)
  rlang::exec(
    gdal_block_apply,
    x, fun, ..., !!!extra,
    radius = radius, filename = filename,
    block_rows = block_rows, block_prep = block_prep
  )
}

# Halo radius in rows for a window argument (scalar or c(rows, cols)).
halo_rows <- function(window) {
  w <- rep(as.integer(window), length.out = 2L)[[1L]]
  (w - 1L) %/% 2L
}

#' @rdname rf_lee
#' @export
rf_lee.Rcpp_GDALRaster <- function(x, window = 7L, ...) {
  filter_gdalraster(x, rf_lee, window = window, ..., radius = halo_rows(window))
}

#' @rdname rf_enhanced_lee
#' @export
rf_enhanced_lee.Rcpp_GDALRaster <- function(x, window = 7L, ...) {
  filter_gdalraster(x, rf_enhanced_lee,
    window = window, ...,
    radius = halo_rows(window)
  )
}

#' @rdname rf_lee_sigma
#' @export
rf_lee_sigma.Rcpp_GDALRaster <- function(x, window = 7L, ...) {
  filter_gdalraster(x, rf_lee_sigma,
    window = window, ...,
    radius = halo_rows(window)
  )
}

#' @rdname rf_lee_sigma_improved
#' @export
rf_lee_sigma_improved.Rcpp_GDALRaster <- function(x, window = 7L, looks = 1,
                                                  sigma = 0.9,
                                                  target_window = 3L,
                                                  edge = c("shrink", "reflect", "nearest", "constant"),
                                                  edge_value = 0,
                                                  na_policy = c("omit", "propagate"),
                                                  ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  # The 98th percentile is global: when streaming, compute it exactly up
  # front (histogram refinement over the blocks) so tiled results match the
  # whole-image run; the eager path lets Rust compute it per layer.
  global_prep <- function(block_rows) {
    z98 <- vapply(
      seq_len(x$dim()[[3L]]),
      function(b) gdal_band_z98(x, b, block_rows),
      numeric(1)
    )
    list(z98 = z98)
  }
  filter_gdalraster(
    x, lee_sigma_improved_impl,
    window = window, looks = looks, sigma = sigma,
    target_window = target_window, edge = edge, edge_value = edge_value,
    na_policy = na_policy, ...,
    radius = max(halo_rows(window), as.integer(target_window) - 1L),
    global_prep = global_prep
  )
}

#' @rdname rf_frost
#' @export
rf_frost.Rcpp_GDALRaster <- function(x, window = 7L, ...) {
  filter_gdalraster(x, rf_frost,
    window = window, ...,
    radius = halo_rows(window)
  )
}

#' @rdname rf_kuan
#' @export
rf_kuan.Rcpp_GDALRaster <- function(x, window = 7L, ...) {
  filter_gdalraster(x, rf_kuan,
    window = window, ...,
    radius = halo_rows(window)
  )
}

#' @rdname rf_gamma_map
#' @export
rf_gamma_map.Rcpp_GDALRaster <- function(x, window = 7L, ...) {
  filter_gdalraster(x, rf_gamma_map,
    window = window, ...,
    radius = halo_rows(window)
  )
}

#' @rdname rf_mean
#' @export
rf_mean.Rcpp_GDALRaster <- function(x, window = 3L, ...) {
  filter_gdalraster(x, rf_mean,
    window = window, ...,
    radius = halo_rows(window)
  )
}

#' @rdname rf_gaussian
#' @export
rf_gaussian.Rcpp_GDALRaster <- function(x, sigma = 1, window = NULL, ...) {
  window <- window %||% (2L * as.integer(ceiling(3 * sigma)) + 1L)
  filter_gdalraster(x, rf_gaussian,
    sigma = sigma, window = window, ...,
    radius = halo_rows(window)
  )
}

#' @rdname rf_median
#' @export
rf_median.Rcpp_GDALRaster <- function(x, window = 3L, ...) {
  filter_gdalraster(x, rf_median,
    window = window, ...,
    radius = halo_rows(window)
  )
}

#' @rdname rf_focal
#' @export
rf_focal.Rcpp_GDALRaster <- function(x, window = 3L, ...) {
  filter_gdalraster(x, rf_focal,
    window = window, ...,
    radius = halo_rows(window)
  )
}

#' @rdname rf_convolve
#' @export
rf_convolve.Rcpp_GDALRaster <- function(x, kernel, ...) {
  kernel <- check_kernel(kernel)
  filter_gdalraster(x, rf_convolve,
    kernel = kernel, ...,
    radius = (nrow(kernel) - 1L) %/% 2L
  )
}

#' @rdname rf_sobel
#' @export
rf_sobel.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_sobel, ..., radius = 1L)
}

#' @rdname rf_laplacian
#' @export
rf_laplacian.Rcpp_GDALRaster <- function(x, ...) {
  filter_gdalraster(x, rf_laplacian, ..., radius = 1L)
}

#' @rdname rf_bilateral
#' @export
rf_bilateral.Rcpp_GDALRaster <- function(x, sigma_d = 1.5, sigma_r = NULL,
                                         window = NULL, ...) {
  window <- window %||% (2L * as.integer(ceiling(2 * sigma_d)) + 1L)
  filter_gdalraster(x, rf_bilateral,
    sigma_d = sigma_d, sigma_r = sigma_r, window = window, ...,
    radius = halo_rows(window)
  )
}

#' @rdname rf_guided
#' @export
rf_guided.Rcpp_GDALRaster <- function(x, guide = NULL, window = 5L, ...) {
  # The guide must be windowed alongside `x`, whether it arrives as another
  # GDALRaster or as an in-memory array, so it is always delivered per block
  # through block_prep (the eager path asks for one full-height block).
  block_prep <- NULL
  if (inherits(guide, "Rcpp_GDALRaster")) {
    gds <- guide
    if (!identical(gds$dim(), x$dim())) {
      cli::cli_abort("{.arg guide} must have the same dimensions as {.arg x}.")
    }
    block_prep <- function(y0, nrows) {
      list(guide = gdal_read_rows(gds, y0, nrows))
    }
  } else if (!is.null(guide)) {
    garr <- check_input(guide)
    if (is.matrix(garr)) garr <- array(garr, dim = c(dim(garr), 1L))
    dm <- as.integer(x$dim())
    if (!identical(dim(garr), c(dm[[2L]], dm[[1L]], dm[[3L]]))) {
      cli::cli_abort("{.arg guide} must have the same dimensions as {.arg x}.")
    }
    block_prep <- function(y0, nrows) {
      list(guide = garr[y0 + seq_len(nrows), , , drop = FALSE])
    }
  }
  # Two chained box means: the dependency radius is a full window minus one.
  filter_gdalraster(
    x, rf_guided,
    window = window, ...,
    radius = 2L * halo_rows(window),
    block_prep = block_prep
  )
}
