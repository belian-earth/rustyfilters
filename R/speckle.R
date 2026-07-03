# SAR speckle filters. All are adaptive local-statistics filters: with m the
# window mean, v its population variance, ci = sqrt(v)/m the observed
# variation coefficient and cu = 1/sqrt(looks) the speckle variation
# coefficient, each filter blends the window mean and the centre pixel by how
# much the local variation exceeds pure speckle. Unlike the smoothing
# filters, a cell whose centre value is NA stays NA: the estimators
# reconstruct the signal at the observed pixel.

#' Lee speckle filter
#'
#' The Lee (1980) minimum mean-square-error filter: each pixel becomes
#' `m + W * (x - m)` with weight `W = v / (v + (m * cu)^2)`, so homogeneous
#' areas take the window mean while high-contrast features stay close to the
#' observed value.
#'
#' @inheritParams rf_params
#' @inherit rf_params return
#' @param looks Single positive number: the effective number of looks (ENL)
#'   of the intensity image. Controls the assumed speckle strength
#'   `cu = 1 / sqrt(looks)`; single-look SAR intensity is `looks = 1`.
#' @details
#' Cells whose centre value is `NA` remain `NA`. Speckle filters expect
#' intensity (linear power) data, not dB.
#' @references Lee, J.-S. (1980). Digital image enhancement and noise
#'   filtering by use of local statistics. *IEEE Transactions on Pattern
#'   Analysis and Machine Intelligence*, 2(2), 165-168.
#' @seealso [rf_enhanced_lee()], [rf_lee_sigma()], [rf_kuan()], [rf_frost()],
#'   [rf_gamma_map()], [rf_set_threads()]
#' @examples
#' set.seed(1)
#' truth <- matrix(rep(c(1, 4), each = 200), 20, 20)
#' speckled <- truth * rexp(400)
#' rf_lee(speckled, window = 7L, looks = 1)
#' @export
rf_lee <- function(x, ...) {
  UseMethod("rf_lee")
}

#' @rdname rf_lee
#' @export
rf_lee.matrix <- function(x, window = 7L, looks = 1,
                          edge = c("shrink", "reflect", "nearest", "constant"),
                          edge_value = 0,
                          na_policy = c("omit", "propagate"),
                          ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  looks <- check_positive_scalar(looks)
  rf_dispatch(x, rf_lee_rs, list(looks), window, edge, edge_value, na_policy)
}

#' @rdname rf_lee
#' @export
rf_lee.array <- rf_lee.matrix

#' @rdname rf_lee
#' @export
rf_lee.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}

#' Enhanced Lee speckle filter
#'
#' The enhanced Lee filter (Lopes et al. 1990) splits each window into three
#' regimes by its variation coefficient `ci`: homogeneous windows
#' (`ci <= cu`) take the window mean, heterogeneous windows blend mean and
#' centre with an exponential damping weight, and point targets
#' (`ci >= cmax = sqrt(1 + 2/looks)`) are preserved untouched.
#'
#' @inheritParams rf_lee
#' @inherit rf_params return
#' @param damping Single positive number controlling how quickly the blend
#'   moves towards the observed value as heterogeneity grows. Larger values
#'   preserve more detail; `1` is the conventional default.
#' @references Lopes, A., Touzi, R., & Nezry, E. (1990). Adaptive speckle
#'   filters and scene heterogeneity. *IEEE Transactions on Geoscience and
#'   Remote Sensing*, 28(6), 992-1000.
#' @seealso [rf_lee()], [rf_gamma_map()], [rf_frost()]
#' @examples
#' set.seed(1)
#' speckled <- matrix(rexp(400), 20, 20)
#' rf_enhanced_lee(speckled, window = 7L, looks = 1, damping = 1)
#' @export
rf_enhanced_lee <- function(x, ...) {
  UseMethod("rf_enhanced_lee")
}

#' @rdname rf_enhanced_lee
#' @export
rf_enhanced_lee.matrix <- function(x, window = 7L, looks = 1, damping = 1,
                                   edge = c("shrink", "reflect", "nearest", "constant"),
                                   edge_value = 0,
                                   na_policy = c("omit", "propagate"),
                                   ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  looks <- check_positive_scalar(looks)
  damping <- check_positive_scalar(damping)
  rf_dispatch(
    x, rf_enhanced_lee_rs, list(looks, damping),
    window, edge, edge_value, na_policy
  )
}

#' @rdname rf_enhanced_lee
#' @export
rf_enhanced_lee.array <- rf_enhanced_lee.matrix

#' @rdname rf_enhanced_lee
#' @export
rf_enhanced_lee.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}

#' Lee sigma speckle filter
#'
#' The classic two-sigma Lee filter (1983): each pixel is replaced by the
#' mean of the window pixels lying within `x * (1 +/- k * cu)`, the range a
#' pure-speckle observation of the centre value would span. If fewer than
#' `min_count` pixels qualify, the full window mean is used instead
#' (suppressing isolated dark/bright noise).
#'
#' @inheritParams rf_lee
#' @inherit rf_params return
#' @param k Single positive number: the sigma multiplier defining the
#'   inclusion range. `2` covers roughly 95.5% of pure speckle.
#' @param min_count Single positive integer: minimum number of in-range
#'   pixels below which the filter falls back to the full window mean.
#' @references Lee, J.-S. (1983). Digital image smoothing and the sigma
#'   filter. *Computer Vision, Graphics, and Image Processing*, 24(2),
#'   255-269.
#' @seealso [rf_lee()], [rf_enhanced_lee()]
#' @examples
#' set.seed(1)
#' speckled <- matrix(rexp(400), 20, 20)
#' rf_lee_sigma(speckled, window = 7L, looks = 1)
#' @export
rf_lee_sigma <- function(x, ...) {
  UseMethod("rf_lee_sigma")
}

#' @rdname rf_lee_sigma
#' @export
rf_lee_sigma.matrix <- function(x, window = 7L, looks = 1, k = 2,
                                min_count = 3L,
                                edge = c("shrink", "reflect", "nearest", "constant"),
                                edge_value = 0,
                                na_policy = c("omit", "propagate"),
                                ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  looks <- check_positive_scalar(looks)
  k <- check_positive_scalar(k)
  min_count <- check_count_scalar(min_count)
  rf_dispatch(
    x, rf_lee_sigma_rs, list(looks, k, min_count),
    window, edge, edge_value, na_policy
  )
}

#' @rdname rf_lee_sigma
#' @export
rf_lee_sigma.array <- rf_lee_sigma.matrix

#' @rdname rf_lee_sigma
#' @export
rf_lee_sigma.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}

#' Frost speckle filter
#'
#' The Frost (1982) filter convolves each window with a damped exponential
#' kernel `exp(-B * d)` where `d` is the distance from the centre and
#' `B = damping * ci^2` adapts to the local variation coefficient:
#' homogeneous windows average broadly, heterogeneous windows concentrate
#' weight on the centre pixel.
#'
#' @inheritParams rf_params
#' @inherit rf_params return
#' @param damping Single positive number scaling the kernel decay. Larger
#'   values preserve more edges and smooth less; `2` is a common default.
#' @references Frost, V. S., Stiles, J. A., Shanmugan, K. S., & Holtzman,
#'   J. C. (1982). A model for radar images and its application to adaptive
#'   digital filtering of multiplicative noise. *IEEE Transactions on
#'   Pattern Analysis and Machine Intelligence*, 4(2), 157-166.
#' @seealso [rf_lee()], [rf_enhanced_lee()], [rf_gamma_map()]
#' @examples
#' set.seed(1)
#' speckled <- matrix(rexp(400), 20, 20)
#' rf_frost(speckled, window = 7L, damping = 2)
#' @export
rf_frost <- function(x, ...) {
  UseMethod("rf_frost")
}

#' @rdname rf_frost
#' @export
rf_frost.matrix <- function(x, window = 7L, damping = 2,
                            edge = c("shrink", "reflect", "nearest", "constant"),
                            edge_value = 0,
                            na_policy = c("omit", "propagate"),
                            ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  damping <- check_positive_scalar(damping)
  rf_dispatch(x, rf_frost_rs, list(damping), window, edge, edge_value, na_policy)
}

#' @rdname rf_frost
#' @export
rf_frost.array <- rf_frost.matrix

#' @rdname rf_frost
#' @export
rf_frost.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}

#' Kuan speckle filter
#'
#' The Kuan (1985) filter: like [rf_lee()] but derived without linearising
#' the multiplicative noise model, giving weight
#' `W = (1 - cu^2/ci^2) / (1 + cu^2)` (clamped to `[0, 1]`).
#'
#' @inheritParams rf_lee
#' @inherit rf_params return
#' @references Kuan, D. T., Sawchuk, A. A., Strand, T. C., & Chavel, P.
#'   (1985). Adaptive noise smoothing filter for images with
#'   signal-dependent noise. *IEEE Transactions on Pattern Analysis and
#'   Machine Intelligence*, 7(2), 165-177.
#' @seealso [rf_lee()], [rf_gamma_map()]
#' @examples
#' set.seed(1)
#' speckled <- matrix(rexp(400), 20, 20)
#' rf_kuan(speckled, window = 7L, looks = 1)
#' @export
rf_kuan <- function(x, ...) {
  UseMethod("rf_kuan")
}

#' @rdname rf_kuan
#' @export
rf_kuan.matrix <- function(x, window = 7L, looks = 1,
                           edge = c("shrink", "reflect", "nearest", "constant"),
                           edge_value = 0,
                           na_policy = c("omit", "propagate"),
                           ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  looks <- check_positive_scalar(looks)
  rf_dispatch(x, rf_kuan_rs, list(looks), window, edge, edge_value, na_policy)
}

#' @rdname rf_kuan
#' @export
rf_kuan.array <- rf_kuan.matrix

#' @rdname rf_kuan
#' @export
rf_kuan.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}

#' Gamma-MAP speckle filter
#'
#' The Gamma maximum a posteriori filter (Lopes et al. 1990), assuming
#' Gamma-distributed signal and speckle. Homogeneous windows (`ci <= cu`)
#' take the window mean, point targets (`ci >= cmax`) are preserved, and in
#' between the MAP estimate solves
#' `out = (B * m + sqrt(m^2 * B^2 + 4 * alpha * looks * m * x)) / (2 * alpha)`
#' with `alpha = (1 + cu^2) / (ci^2 - cu^2)` and `B = alpha - looks - 1`.
#'
#' @inheritParams rf_lee
#' @inherit rf_params return
#' @references Lopes, A., Nezry, E., Touzi, R., & Laur, H. (1990). Maximum a
#'   posteriori speckle filtering and first order texture models in SAR
#'   images. *IGARSS 1990*, 2409-2412.
#' @seealso [rf_lee()], [rf_enhanced_lee()], [rf_frost()]
#' @examples
#' set.seed(1)
#' speckled <- matrix(rexp(400), 20, 20)
#' rf_gamma_map(speckled, window = 7L, looks = 1)
#' @export
rf_gamma_map <- function(x, ...) {
  UseMethod("rf_gamma_map")
}

#' @rdname rf_gamma_map
#' @export
rf_gamma_map.matrix <- function(x, window = 7L, looks = 1,
                                edge = c("shrink", "reflect", "nearest", "constant"),
                                edge_value = 0,
                                na_policy = c("omit", "propagate"),
                                ...) {
  edge <- rlang::arg_match(edge)
  na_policy <- rlang::arg_match(na_policy)
  looks <- check_positive_scalar(looks)
  rf_dispatch(x, rf_gamma_map_rs, list(looks), window, edge, edge_value, na_policy)
}

#' @rdname rf_gamma_map
#' @export
rf_gamma_map.array <- rf_gamma_map.matrix

#' @rdname rf_gamma_map
#' @export
rf_gamma_map.default <- function(x, ...) {
  cli::cli_abort(
    "{.arg x} must be a numeric matrix or 3-D array,
     not {.obj_type_friendly {x}}."
  )
}
