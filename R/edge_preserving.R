# Edge-preserving smoothers: bilateral and guided filters.

#' Bilateral filter
#'
#' Smooths with weights that are the product of a spatial Gaussian (distance
#' from the centre cell, `sigma_d`) and a range Gaussian (difference from the
#' centre value, `sigma_r`), so averaging happens within regions of similar
#' value while sharp transitions survive.
#'
#' @inheritParams rf_params
#' @inherit rf_params return
#' @param sigma_d Single positive number: the spatial standard deviation in
#'   cells.
#' @param sigma_r Single positive number: the range standard deviation in
#'   value units. Values differing from the centre by much more than
#'   `sigma_r` are effectively excluded. The default `NULL` uses the standard
#'   deviation of the valid cells of `x`, a serviceable starting point that
#'   you should expect to tune.
#' @param window Window size: a single odd integer or a `c(rows, cols)`
#'   pair. The default `NULL` uses `2 * ceiling(2 * sigma_d) + 1`.
#' @details
#' Cells whose centre value is `NA` stay `NA` (the range weight needs the
#' centre). Cost grows with the window area; this is the slowest smoother in
#' the package.
#' @seealso [rf_guided()], [rf_gaussian()], [rf_median()],
#'   [rf_set_threads()]
#' @examples
#' noisy <- volcano + matrix(rnorm(length(volcano), sd = 8), nrow(volcano))
#' op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
#' rf_plot(noisy, main = "noisy volcano")
#' rf_plot(rf_bilateral(noisy, sigma_d = 2, sigma_r = 20), main = "bilateral")
#' par(op)
#' @export
rf_bilateral <- function(x, ...) {
  UseMethod("rf_bilateral")
}

#' @rdname rf_bilateral
#' @export
rf_bilateral.matrix <- function(x, sigma_d = 1.5, sigma_r = NULL, window = NULL,
                                edge = c("shrink", "reflect", "nearest", "constant"),
                                edge_value = 0,
                                na_policy = c("omit", "propagate"),
                                ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  sigma_d <- check_positive_scalar(sigma_d)
  if (is.null(sigma_r)) {
    sigma_r <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(sigma_r) || sigma_r <= 0) {
      cli::cli_abort(
        "Cannot infer {.arg sigma_r} from {.arg x}; supply it explicitly."
      )
    }
  }
  sigma_r <- check_positive_scalar(sigma_r)
  window <- window %||% (2L * as.integer(ceiling(2 * sigma_d)) + 1L)
  rf_dispatch(
    x, rf_bilateral_rs, list(sigma_d, sigma_r),
    window, edge, edge_value, na_policy
  )
}

#' @rdname rf_bilateral
#' @export
rf_bilateral.array <- rf_bilateral.matrix

#' @rdname rf_bilateral
#' @export
rf_bilateral.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}

#' Guided filter
#'
#' The guided filter (He et al. 2013) fits a local linear model between a
#' guide image and the input in every window, then averages the coefficients:
#' an edge-preserving smoother whose cost is independent of the window size
#' (it runs entirely on the package's O(1) box-mean engine). With the default
#' self-guidance (`guide = NULL`) it behaves like a fast bilateral
#' alternative; with a separate `guide` it transfers the guide's structure
#' onto `x` (e.g. detail-preserving smoothing of a noisy band guided by a
#' clean one).
#'
#' @inheritParams rf_params
#' @inherit rf_params return
#' @param guide A matrix or 3-D array with the same dimensions as `x`, or
#'   `NULL` (the default) to use `x` itself.
#' @param eps Single positive number: the regularisation added to the window
#'   variance of the guide. Windows whose variance is well below `eps` are
#'   smoothed; windows well above it (edges) are preserved. Scale-dependent:
#'   the default `NULL` uses `(0.1 * sd(guide))^2`.
#' @param window Window size: a single odd integer or a `c(rows, cols)`
#'   pair.
#' @details
#' Missing cells in either `x` or the guide are excluded from every window
#' fit; cells whose centre is `NA` in the guide stay `NA` in the output.
#' @references He, K., Sun, J., & Tang, X. (2013). Guided image filtering.
#'   *IEEE Transactions on Pattern Analysis and Machine Intelligence*,
#'   35(6), 1397-1409.
#' @seealso [rf_bilateral()], [rf_mean()], [rf_set_threads()]
#' @examples
#' noisy <- volcano + matrix(rnorm(length(volcano), sd = 8), nrow(volcano))
#' op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
#' rf_plot(noisy, main = "noisy volcano")
#' rf_plot(rf_guided(noisy, window = 5L, eps = 64), main = "guided")
#' par(op)
#' @export
rf_guided <- function(x, ...) {
  UseMethod("rf_guided")
}

#' @rdname rf_guided
#' @export
rf_guided.matrix <- function(x, guide = NULL, window = 5L, eps = NULL,
                             edge = c("shrink", "reflect", "nearest", "constant"),
                             edge_value = 0,
                             na_policy = c("omit", "propagate"),
                             ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  guide <- guide %||% x
  guide <- check_input(guide)
  if (!identical(dim(guide), dim(check_input(x)))) {
    cli::cli_abort(
      "{.arg guide} must have the same dimensions as {.arg x}."
    )
  }
  if (is.null(eps)) {
    eps <- (0.1 * stats::sd(guide, na.rm = TRUE))^2
    if (!is.finite(eps) || eps <= 0) {
      cli::cli_abort(
        "Cannot infer {.arg eps} from the guide; supply it explicitly."
      )
    }
  }
  eps <- check_positive_scalar(eps)
  rf_dispatch(
    x, rf_guided_rs, list(guide, eps),
    window, edge, edge_value, na_policy
  )
}

#' @rdname rf_guided
#' @export
rf_guided.array <- rf_guided.matrix

#' @rdname rf_guided
#' @export
rf_guided.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}
