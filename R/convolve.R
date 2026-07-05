#' Filter with an arbitrary kernel
#'
#' Applies a user-supplied kernel over a moving window: each output cell is
#' the sum of the window values times the matching kernel weights. This gives
#' sharpening, embossing, gradients or any custom linear filter;
#' [rf_sobel()] and [rf_laplacian()] are thin wrappers with fixed kernels.
#'
#' @inheritParams rf_params
#' @inherit rf_params return
#' @param kernel A numeric matrix with odd dimensions and no missing values.
#'   Applied as-is (cross-correlation): the kernel is not flipped, which only
#'   matters for asymmetric kernels.
#' @param normalize If `TRUE`, divide each result by the sum of the kernel
#'   weights actually used, so partial windows (at edges under
#'   `edge = "shrink"`, or around missing values) keep the input's scale.
#'   Sensible for smoothing kernels; leave `FALSE` for derivative kernels
#'   whose weights sum to zero.
#' @details
#' Under `na_policy = "omit"`, missing cells simply drop out of the weighted
#' sum. For zero-sum kernels this biases results near missing values; there
#' is no principled general correction, so consider `na_policy =
#' "propagate"` for derivative kernels on gappy data.
#' @seealso [rf_sobel()], [rf_laplacian()], [rf_gaussian()],
#'   [rf_set_threads()]
#' @examples
#' sharpen <- matrix(c(0, -1, 0, -1, 5, -1, 0, -1, 0), 3)
#' op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
#' rf_plot(volcano, main = "volcano")
#' rf_plot(rf_convolve(volcano, sharpen), main = "sharpened")
#' par(op)
#' @export
rf_convolve <- function(x, ...) {
  UseMethod("rf_convolve")
}

#' @rdname rf_convolve
#' @export
rf_convolve.matrix <- function(x, kernel, normalize = FALSE,
                               edge = c("shrink", "reflect", "nearest", "constant"),
                               edge_value = 0,
                               na_policy = c("omit", "propagate"),
                               ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  kernel <- check_kernel(kernel)
  if (!rlang::is_bool(normalize)) {
    cli::cli_abort("{.arg normalize} must be `TRUE` or `FALSE`.")
  }
  rf_dispatch(
    x, rf_convolve_rs, list(kernel, normalize),
    dim(kernel), edge, edge_value, na_policy
  )
}

#' @rdname rf_convolve
#' @export
rf_convolve.array <- rf_convolve.matrix

#' @rdname rf_convolve
#' @export
rf_convolve.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}

#' Sobel edge detection
#'
#' Convolves with the 3 x 3 Sobel kernels to estimate the local gradient:
#' `"x"` is the rate of change across columns, `"y"` across rows, and
#' `"magnitude"` (the default) is `sqrt(gx^2 + gy^2)`.
#'
#' @inheritParams rf_params
#' @inherit rf_params return
#' @param direction `"magnitude"`, `"x"` or `"y"`.
#' @param edge Edge policy, as in [rf_convolve()]. Defaults to `"nearest"`
#'   for gradients: replicating the border cell gives a near-zero gradient
#'   there, whereas a shrinking window would fabricate strong edges.
#' @seealso [rf_convolve()], [rf_laplacian()]
#' @examples
#' op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
#' rf_plot(volcano, main = "volcano")
#' rf_plot(rf_sobel(volcano), main = "Sobel gradient magnitude")
#' par(op)
#' @export
rf_sobel <- function(x, ...) {
  UseMethod("rf_sobel")
}

#' @rdname rf_sobel
#' @export
rf_sobel.matrix <- function(x, direction = c("magnitude", "x", "y"),
                            edge = c("nearest", "shrink", "reflect", "constant"),
                            edge_value = 0,
                            na_policy = c("omit", "propagate"),
                            ...) {
  direction <- rlang::arg_match(direction)
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  gx_kernel <- matrix(c(-1, -2, -1, 0, 0, 0, 1, 2, 1), 3)
  conv <- function(k) {
    rf_dispatch(x, rf_convolve_rs, list(k, FALSE), 3L, edge, edge_value, na_policy)
  }
  switch(direction,
    x = conv(gx_kernel),
    y = conv(t(gx_kernel)),
    magnitude = sqrt(conv(gx_kernel)^2 + conv(t(gx_kernel))^2)
  )
}

#' @rdname rf_sobel
#' @export
rf_sobel.array <- rf_sobel.matrix

#' @rdname rf_sobel
#' @export
rf_sobel.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}

#' Laplacian filter
#'
#' Convolves with the discrete Laplacian, highlighting local extrema and
#' rapid changes (zero over flat and linearly varying regions).
#'
#' @inheritParams rf_sobel
#' @inherit rf_params return
#' @param neighbours `4L` for the cross-shaped kernel
#'   (`c(0, 1, 0, 1, -4, 1, 0, 1, 0)`), `8L` to include diagonals.
#' @seealso [rf_convolve()], [rf_sobel()]
#' @examples
#' op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
#' rf_plot(volcano, main = "volcano")
#' rf_plot(rf_laplacian(volcano), main = "Laplacian")
#' par(op)
#' @export
rf_laplacian <- function(x, ...) {
  UseMethod("rf_laplacian")
}

#' @rdname rf_laplacian
#' @export
rf_laplacian.matrix <- function(x, neighbours = 4L,
                                edge = c("nearest", "shrink", "reflect", "constant"),
                                edge_value = 0,
                                na_policy = c("omit", "propagate"),
                                ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  neighbours <- vctrs::vec_cast(neighbours, integer())
  if (!length(neighbours) == 1L || !neighbours %in% c(4L, 8L)) {
    cli::cli_abort("{.arg neighbours} must be 4 or 8.")
  }
  k <- if (neighbours == 4L) {
    matrix(c(0, 1, 0, 1, -4, 1, 0, 1, 0), 3)
  } else {
    matrix(c(1, 1, 1, 1, -8, 1, 1, 1, 1), 3)
  }
  rf_dispatch(x, rf_convolve_rs, list(k, FALSE), 3L, edge, edge_value, na_policy)
}

#' @rdname rf_laplacian
#' @export
rf_laplacian.array <- rf_laplacian.matrix

#' @rdname rf_laplacian
#' @export
rf_laplacian.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}
