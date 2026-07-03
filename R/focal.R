#' Focal statistics over a moving window
#'
#' Computes a summary statistic of the cells in a moving window centred on
#' each cell. All statistics share the same multi-threaded engine; see
#' [rf_mean()], [rf_median()] and [rf_gaussian()] for the dedicated smoothing
#' filters and [rf_lee()] and friends for adaptive speckle filters.
#'
#' @inheritParams rf_params
#' @param stat Statistic to compute for each window. One of `"mean"`,
#'   `"median"`, `"min"`, `"max"`, `"range"` (max minus min), `"sd"` (sample
#'   standard deviation), `"sum"` or `"mode"`.
#' @details
#' `"mode"` returns the most frequent value using exact floating-point
#' equality, with ties resolved to the lowest value. It is intended for
#' categorical data encoded as numbers and is not meaningful for continuous
#' data.
#' @returns An object of the same class and dimensions as `x` (dimnames are
#'   preserved), containing the filtered values as doubles.
#' @seealso [rf_mean()], [rf_median()], [rf_gaussian()], [rf_set_threads()]
#' @examples
#' m <- matrix(as.numeric(1:20), nrow = 4)
#' rf_focal(m, window = 3L, stat = "sd")
#' rf_focal(m, window = c(3L, 5L), stat = "max", edge = "nearest")
#' @export
rf_focal <- function(x, ...) {
  UseMethod("rf_focal")
}

#' @rdname rf_focal
#' @export
rf_focal.matrix <- function(x, window = 3L,
                            stat = c(
                              "mean", "median", "min", "max",
                              "range", "sd", "sum", "mode"
                            ),
                            edge = c("shrink", "reflect", "nearest", "constant"),
                            edge_value = 0,
                            na_policy = c("omit", "propagate"),
                            ...) {
  stat <- rlang::arg_match(stat)
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  rf_dispatch(x, rf_focal_rs, list(stat), window, edge, edge_value, na_policy)
}

#' @rdname rf_focal
#' @export
rf_focal.array <- rf_focal.matrix

#' @rdname rf_focal
#' @export
rf_focal.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}
