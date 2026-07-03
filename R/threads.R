#' Control the number of threads used by rustyfilters
#'
#' All filters process windows in parallel via rayon. On load the thread count
#' defaults to all available cores, unless overridden by
#' `option(rustyfilters.threads = n)` or the `RUSTYFILTERS_NUM_THREADS`
#' environment variable (checked in that order), or capped at 2 when R is
#' running checks with `_R_CHECK_LIMIT_CORES_` set. Set to 1 for sequential
#' processing with zero threading overhead.
#'
#' @param n Integer scalar. Number of threads. Must be >= 1.
#' @returns `rf_set_threads()` invisibly returns the previous thread count;
#'   `rf_get_threads()` returns the current thread count as an integer scalar.
#' @examples
#' old <- rf_set_threads(2L)
#' rf_get_threads()
#' rf_set_threads(old)
#' @rdname rf_threads
#' @export
rf_set_threads <- function(n = 1L) {
  n <- vctrs::vec_cast(n, integer())
  vctrs::vec_assert(n, size = 1L)
  if (is.na(n) || n < 1L) {
    cli::cli_abort("{.arg n} must be a positive integer scalar.")
  }
  old <- rf_get_threads_rs()
  rf_set_threads_rs(n)
  invisible(old)
}

#' @rdname rf_threads
#' @export
rf_get_threads <- function() {
  rf_get_threads_rs()
}
