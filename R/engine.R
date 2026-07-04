# Shared validate -> Rust -> restore-attributes funnel used by every filter
# method. `args` carries the filter-specific parameters in the order the
# `*_rs` function expects them (between the window and edge arguments).
rf_dispatch <- function(x, rs_fun, args, window, edge, edge_value, na_policy,
                        call = rlang::caller_env()) {
  x <- check_input(x, call = call)
  w <- check_window(window, call = call)
  edge_value <- check_edge_value(edge_value, call = call)
  d <- dim(x)
  nl <- if (length(d) == 3L) d[[3L]] else 1L
  out <- rlang::exec(
    rs_fun,
    x, d[[1L]], d[[2L]], nl, w[[1L]], w[[2L]], !!!args,
    edge, edge_value, na_policy == "omit"
  )
  dim(out) <- d
  dimnames(out) <- dimnames(x)
  out
}

# Shared @param blocks -----------------------------------------------------

#' Common filter parameters
#'
#' @param x A numeric matrix or 3-D array (filtered layer by layer). Methods
#'   are also provided for terra `SpatRaster` objects (when terra is
#'   installed) and for open gdalraster `GDALRaster` datasets (when
#'   gdalraster is installed). `GDALRaster` methods return a new
#'   `GDALRaster` object open in update mode on a Float64 dataset with the
#'   source's geometry. Small datasets are filtered in memory and land on an
#'   in-memory `/vsimem` GTiff by default; datasets whose decoded size
#'   exceeds `options(rustyfilters.block_memory)` (default 2 GiB) stream
#'   through full-width row bands with a halo sized to the filter's window,
#'   writing to a GeoTIFF tempfile instead. Interior band seams are exact
#'   (the halo supplies the true neighbouring data; `edge` fires only at
#'   real raster edges). `GDALRaster` methods accept three extra arguments:
#'   `filename` (output path, replacing the tempfile/`/vsimem` default),
#'   `by_block` (`TRUE`/`FALSE` to force or forbid streaming) and
#'   `block_rows` (rows per band, sized from the memory budget by default).
#' @param window Window size in cells: a single odd positive integer, or a
#'   pair `c(rows, cols)` of odd positive integers.
#' @param edge How to treat windows that overhang the matrix edge:
#'   `"shrink"` (default) truncates the window to the cells that exist;
#'   `"reflect"` mirrors the matrix across the boundary; `"nearest"` repeats
#'   the closest edge cell; `"constant"` pads with `edge_value`.
#' @param edge_value Single number used to pad when `edge = "constant"`.
#'   `NA` is allowed and behaves as missing data under `na_policy = "omit"`.
#' @param na_policy How to treat missing values inside a window. The default
#'   `"omit"` excludes them from the statistics (like `na.rm = TRUE`); a
#'   window with no valid cells yields `NA`. `"propagate"` is the fast path:
#'   no per-cell NA handling is compiled into the inner loop, and any `NA`
#'   in a window makes the result `NA`. Use it when the input has no missing
#'   values (or when spreading `NA` is acceptable) for maximum speed.
#' @param ... Passed on to methods.
#' @returns An object of the same class and dimensions as `x` (dimnames are
#'   preserved), containing the filtered values as doubles.
#' @keywords internal
#' @name rf_params
NULL
