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
