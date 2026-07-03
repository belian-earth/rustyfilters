# nocov start
default_threads <- function() {
  n <- getOption("rustyfilters.threads")
  if (is.null(n)) {
    env <- Sys.getenv("RUSTYFILTERS_NUM_THREADS", unset = "")
    if (nzchar(env)) n <- suppressWarnings(as.integer(env))
  }
  if (is.null(n)) {
    # CRAN and friends limit checks to 2 cores; respect that. Otherwise 0L,
    # which the Rust side resolves to all available cores.
    limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = "")
    n <- if (nzchar(limit) && !identical(tolower(limit), "false")) 2L else 0L
  }
  n <- suppressWarnings(as.integer(n))
  if (is.na(n) || n < 0L) n <- 1L
  n
}

.onLoad <- function(libname, pkgname) {
  rf_set_threads_rs(default_threads())
}
# nocov end
