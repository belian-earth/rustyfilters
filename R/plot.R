#' Plot a matrix or array layer as an image
#'
#' A small convenience around [graphics::image()] used throughout the
#' package examples: correct raster orientation (row 1 at the top), square
#' pixels, no axes, a percentile stretch that keeps skewed data (like SAR
#' backscatter) readable, and an Inferno palette from
#' [grDevices::hcl.colors()].
#'
#' @param x A numeric matrix or 3-D array.
#' @param layer Which layer of a 3-D array to plot.
#' @param palette A palette name understood by [grDevices::hcl.colors()].
#' @param stretch Two probabilities: values outside these quantiles are
#'   clipped before plotting. Use `c(0, 1)` for no stretch.
#' @param ... Passed on to [graphics::image()] (e.g. `main`).
#' @returns `x`, invisibly.
#' @examples
#' rf_plot(volcano, main = "Maungawhau")
#' rf_plot(s1_sthelens, main = "Mount St Helens, Sentinel-1 VV")
#' @export
rf_plot <- function(x, layer = 1L, palette = "Inferno", stretch = c(0.02, 0.98),
                    ...) {
  if (is.array(x) && length(dim(x)) == 3L) {
    x <- x[, , layer]
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort(
      "{.arg x} must be a numeric matrix or 3-D array,
       not {.obj_type_friendly {x}}."
    )
  }
  stretch <- vctrs::vec_cast(stretch, double())
  if (length(stretch) != 2L || anyNA(stretch) || any(stretch < 0 | stretch > 1)) {
    cli::cli_abort("{.arg stretch} must be two probabilities.")
  }
  zlim <- stats::quantile(x, sort(stretch), na.rm = TRUE, names = FALSE)
  if (!all(is.finite(zlim)) || zlim[1] == zlim[2]) {
    zlim <- range(x, finite = TRUE)
    if (zlim[1] == zlim[2]) zlim <- zlim + c(-0.5, 0.5)
  }
  z <- pmin(pmax(x, zlim[1]), zlim[2])
  # image() draws z[i, j] with i along x and j upwards: transpose and flip
  # so matrix row 1 is the top row of the plot.
  z <- t(z)
  z <- z[, rev(seq_len(ncol(z))), drop = FALSE]
  graphics::image(
    x = seq_len(nrow(z)), y = seq_len(ncol(z)), z = z,
    col = grDevices::hcl.colors(256, palette),
    asp = 1, axes = FALSE, xlab = "", ylab = "", useRaster = TRUE,
    ...
  )
  invisible(x)
}
