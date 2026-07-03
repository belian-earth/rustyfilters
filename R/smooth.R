#' Moving-window mean (boxcar) filter
#'
#' Smooths by replacing each cell with the mean of its window. Runs on the
#' separable sliding-sum engine, so cost per cell is independent of the
#' window size.
#'
#' @inheritParams rf_params
#' @inherit rf_params return
#' @seealso [rf_median()], [rf_gaussian()], [rf_focal()] for other window
#'   statistics, [rf_set_threads()]
#' @examples
#' m <- matrix(as.numeric(1:25), 5)
#' rf_mean(m)
#' rf_mean(m, window = c(3L, 5L), edge = "reflect")
#' @export
rf_mean <- function(x, ...) {
  UseMethod("rf_mean")
}

#' @rdname rf_mean
#' @export
rf_mean.matrix <- function(x, window = 3L,
                           edge = c("shrink", "reflect", "nearest", "constant"),
                           edge_value = 0,
                           na_policy = c("omit", "propagate"),
                           ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  rf_dispatch(x, rf_focal_rs, list("mean"), window, edge, edge_value, na_policy)
}

#' @rdname rf_mean
#' @export
rf_mean.array <- rf_mean.matrix

#' @rdname rf_mean
#' @export
rf_mean.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}

#' Gaussian blur
#'
#' Smooths with a separable Gaussian kernel: two 1-D passes, so cost grows
#' with the window edge length rather than its area. Under
#' `na_policy = "omit"` this is normalized convolution: missing cells are
#' excluded and the remaining weights are rescaled, which also fills isolated
#' `NA` holes with their neighbourhood estimate.
#'
#' @inheritParams rf_params
#' @inherit rf_params return
#' @param sigma Single positive number: the Gaussian standard deviation in
#'   cells.
#' @param window Kernel size in cells: a single odd positive integer or a
#'   pair `c(rows, cols)`. The default `NULL` uses
#'   `2 * ceiling(3 * sigma) + 1` in both dimensions, which captures
#'   effectively all of the kernel mass.
#' @seealso [rf_mean()], [rf_median()], [rf_focal()], [rf_set_threads()]
#' @examples
#' m <- matrix(as.numeric(1:25), 5)
#' rf_gaussian(m, sigma = 1)
#' rf_gaussian(m, sigma = 2, window = 7L, edge = "reflect")
#' @export
rf_gaussian <- function(x, ...) {
  UseMethod("rf_gaussian")
}

#' @rdname rf_gaussian
#' @export
rf_gaussian.matrix <- function(x, sigma = 1, window = NULL,
                               edge = c("shrink", "reflect", "nearest", "constant"),
                               edge_value = 0,
                               na_policy = c("omit", "propagate"),
                               ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  sigma <- check_positive_scalar(sigma)
  window <- window %||% (2L * as.integer(ceiling(3 * sigma)) + 1L)
  rf_dispatch(x, rf_gaussian_rs, list(sigma), window, edge, edge_value, na_policy)
}

#' @rdname rf_gaussian
#' @export
rf_gaussian.array <- rf_gaussian.matrix

#' @rdname rf_gaussian
#' @export
rf_gaussian.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}

#' Moving-window median filter
#'
#' Replaces each cell with the median of its window: a robust smoother that
#' preserves edges and removes salt-and-pepper noise. Windows with an even
#' number of valid cells (possible at edges or around missing values) average
#' the two middle values, matching [stats::median()].
#'
#' @inheritParams rf_params
#' @inherit rf_params return
#' @seealso [rf_mean()], [rf_gaussian()], [rf_focal()], [rf_set_threads()]
#' @examples
#' m <- matrix(as.numeric(1:25), 5)
#' m[3, 3] <- 1000
#' rf_median(m)
#' @export
rf_median <- function(x, ...) {
  UseMethod("rf_median")
}

#' @rdname rf_median
#' @export
rf_median.matrix <- function(x, window = 3L,
                             edge = c("shrink", "reflect", "nearest", "constant"),
                             edge_value = 0,
                             na_policy = c("omit", "propagate"),
                             ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  rf_dispatch(x, rf_focal_rs, list("median"), window, edge, edge_value, na_policy)
}

#' @rdname rf_median
#' @export
rf_median.array <- rf_median.matrix

#' @rdname rf_median
#' @export
rf_median.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}
