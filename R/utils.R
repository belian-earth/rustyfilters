# Argument validation helpers shared by the public filter functions. Each
# aborts with a cli error attributed to the public caller (`call`), so the
# user sees e.g. "in `rf_lee()`" rather than an internal frame.

# A numeric matrix or 3-D array; integer storage is cast to double.
check_input <- function(x, call = rlang::caller_env()) {
  ok <- is.numeric(x) &&
    (is.matrix(x) || (is.array(x) && length(dim(x)) == 3L))
  if (!ok) {
    cli::cli_abort(
      "{.arg x} must be a numeric matrix or 3-D array,
       not {.obj_type_friendly {x}}.",
      call = call
    )
  }
  if (!is.double(x)) storage.mode(x) <- "double"
  x
}

# One or two odd positive integers; a scalar is recycled to c(rows, cols).
check_window <- function(window, call = rlang::caller_env()) {
  window <- vctrs::vec_cast(window, integer(), call = call)
  bad <- !length(window) %in% 1:2 || anyNA(window) ||
    any(window < 1L) || any(window %% 2L == 0L)
  if (bad) {
    cli::cli_abort(
      "{.arg window} must be one or two odd positive integers.",
      call = call
    )
  }
  vctrs::vec_recycle(window, 2L, call = call)
}

# A single positive finite number (looks, damping, sigma, k, ...).
check_positive_scalar <- function(x, arg = rlang::caller_arg(x),
                                  call = rlang::caller_env()) {
  x <- vctrs::vec_cast(x, double(), x_arg = arg, call = call)
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    cli::cli_abort(
      "{.arg {arg}} must be a single positive number.",
      call = call
    )
  }
  x
}

# A single positive integer (min_count, ...).
check_count_scalar <- function(x, arg = rlang::caller_arg(x),
                               call = rlang::caller_env()) {
  x <- vctrs::vec_cast(x, integer(), x_arg = arg, call = call)
  if (length(x) != 1L || is.na(x) || x < 1L) {
    cli::cli_abort(
      "{.arg {arg}} must be a single positive integer.",
      call = call
    )
  }
  x
}

# A numeric kernel matrix with odd dimensions and finite weights.
check_kernel <- function(kernel, call = rlang::caller_env()) {
  ok <- is.matrix(kernel) && is.numeric(kernel) &&
    all(dim(kernel) %% 2L == 1L) && all(is.finite(kernel))
  if (!ok) {
    cli::cli_abort(
      "{.arg kernel} must be a numeric matrix with odd dimensions and no
       missing values.",
      call = call
    )
  }
  storage.mode(kernel) <- "double"
  kernel
}

# A single number used to pad the `constant` edge policy; NA is allowed and
# behaves as missing data under `na_policy = "omit"`.
check_edge_value <- function(edge_value, call = rlang::caller_env()) {
  edge_value <- vctrs::vec_cast(edge_value, double(), call = call)
  if (length(edge_value) != 1L) {
    cli::cli_abort(
      "{.arg edge_value} must be a single number.",
      call = call
    )
  }
  # The Rust boundary takes a plain f64; NA and NaN are equivalent there.
  if (is.na(edge_value)) edge_value <- NaN
  edge_value
}
