# Out-of-core block engine for GDALRaster input. Rasters are processed as
# full-width row bands with a halo of `radius` rows above and below each
# band. Every filter in the package is local (an output cell depends only on
# input cells within a fixed radius), so filtering the padded band and
# cropping the halo reproduces the whole-image result at interior seams
# exactly; the user's edge policy fires only at true raster edges, where no
# halo exists. The moments fast path restarts its sliding sums at each band,
# so block results can differ from the in-memory run in the last float ulps;
# for a fixed block height they remain fully reproducible.

# Read rows [y0, y0 + nrows) of every band as a (rows, cols, bands) array.
# gdalraster masks nodata to NA on read.
gdal_read_rows <- function(x, y0, nrows) {
  dm <- x$dim()
  v <- lapply(seq_len(dm[[3L]]), function(b) {
    as.numeric(x$read(
      band = b, xoff = 0L, yoff = y0, xsize = dm[[1L]], ysize = nrows,
      out_xsize = dm[[1L]], out_ysize = nrows
    ))
  })
  a <- array(unlist(v, use.names = FALSE), dim = c(dm[[1L]], nrows, dm[[3L]]))
  aperm(a, c(2L, 1L, 3L))
}

# Write a (rows, cols, bands) array to rows [y0, y0 + nrows) of `dst`.
gdal_write_rows <- function(dst, a, y0) {
  d <- dim(a)
  for (b in seq_len(d[[3L]])) {
    dst$write(
      band = b, xoff = 0L, yoff = y0, xsize = d[[2L]], ysize = d[[1L]],
      rasterData = as.numeric(t(a[, , b]))
    )
  }
  invisible(NULL)
}

# Create the Float64 output dataset with the source geometry, open in update
# mode with NaN nodata on every band.
gdal_create_like <- function(src_file, filename, nbands) {
  gdalraster::rasterFromRaster(src_file, filename,
    fmt = "GTiff", dtName = "Float64"
  )
  dst <- methods::new(gdalraster::GDALRaster, filename, read_only = FALSE)
  for (b in seq_len(nbands)) {
    dst$setNoDataValue(b, NaN)
  }
  dst
}

# Number of rows per band: a quarter of the memory budget, leaving headroom
# for the halo, the filtered copy and any per-filter temporaries.
gdal_block_rows <- function(dm, budget) {
  bytes_per_row <- dm[[1L]] * dm[[3L]] * 8
  rows <- floor((budget / 4) / bytes_per_row)
  as.integer(min(dm[[2L]], max(rows, 64L)))
}

# Stream `fun` over row bands of `x` with a halo of `radius` rows.
# `block_prep(y0, nrows)` may supply extra per-band arguments to `fun` (the
# cross-guided filter reads the matching guide window with it).
gdal_block_apply <- function(x, fun, ..., radius, filename, block_rows,
                             block_prep = NULL) {
  dm <- x$dim()
  dst <- gdal_create_like(x$getFilename(), filename, dm[[3L]])
  y0 <- 0L
  while (y0 < dm[[2L]]) {
    nr_block <- min(block_rows, dm[[2L]] - y0)
    top <- min(radius, y0)
    bottom <- min(radius, dm[[2L]] - (y0 + nr_block))
    a <- gdal_read_rows(x, y0 - top, nr_block + top + bottom)
    extra <- if (is.null(block_prep)) {
      list()
    } else {
      block_prep(y0 - top, nr_block + top + bottom)
    }
    out <- rlang::exec(fun, a, ..., !!!extra)
    keep <- seq.int(top + 1L, top + nr_block)
    gdal_write_rows(dst, out[keep, , , drop = FALSE], y0)
    y0 <- y0 + nr_block
  }
  dst
}

# ---------------------------------------------------------------------------
# Exact streaming order statistic: the k-th smallest valid value of one band,
# found by histogram refinement. Each round bins the current value range,
# descends into the bin containing rank k, and finishes with an exact sort
# once few enough values remain, so only O(block) memory is ever held while
# the result equals sort(values)[k] exactly.

gdal_band_kth <- function(x, band, k, lo, hi, n_in_range, block_rows) {
  dm <- x$dim()
  each_block <- function(f) {
    y0 <- 0L
    while (y0 < dm[[2L]]) {
      nr_block <- min(block_rows, dm[[2L]] - y0)
      v <- as.numeric(x$read(
        band = band, xoff = 0L, yoff = y0, xsize = dm[[1L]], ysize = nr_block,
        out_xsize = dm[[1L]], out_ysize = nr_block
      ))
      f(v[!is.na(v)])
      y0 <- y0 + nr_block
    }
  }
  hi_closed <- TRUE
  in_range <- function(v) {
    v[v >= lo & (if (hi_closed) v <= hi else v < hi)]
  }
  collect_limit <- 2^22
  nbins <- 65536L
  repeat {
    if (lo == hi || n_in_range <= collect_limit) {
      vals <- numeric(0)
      each_block(function(v) vals <<- c(vals, in_range(v)))
      return(sort(vals)[k])
    }
    breaks <- seq(lo, hi, length.out = nbins + 1L)
    counts <- numeric(nbins)
    each_block(function(v) {
      bins <- findInterval(in_range(v), breaks,
        rightmost.closed = TRUE, all.inside = TRUE
      )
      counts <<- counts + tabulate(bins, nbins = nbins)
    })
    j <- which(cumsum(counts) >= k)[1L]
    if (is.na(j)) {
      # Numerical shortfall (should not happen): fall back to the top bin.
      j <- nbins
    }
    if (j > 1L) k <- k - sum(counts[seq_len(j - 1L)])
    n_in_range <- counts[[j]]
    hi_closed <- hi_closed && j == nbins
    lo <- breaks[[j]]
    hi <- breaks[[j + 1L]]
  }
}

# Per-band exact 98th percentile, matching the Rust definition:
# sort(valid)[min(max(floor(0.98 * n), 1), n)].
gdal_band_z98 <- function(x, band, block_rows) {
  dm <- x$dim()
  n <- 0
  lo <- Inf
  hi <- -Inf
  y0 <- 0L
  while (y0 < dm[[2L]]) {
    nr_block <- min(block_rows, dm[[2L]] - y0)
    v <- as.numeric(x$read(
      band = band, xoff = 0L, yoff = y0, xsize = dm[[1L]], ysize = nr_block,
      out_xsize = dm[[1L]], out_ysize = nr_block
    ))
    v <- v[!is.na(v)]
    n <- n + length(v)
    if (length(v)) {
      lo <- min(lo, min(v))
      hi <- max(hi, max(v))
    }
    y0 <- y0 + nr_block
  }
  if (n == 0) {
    return(Inf)
  }
  k <- min(max(floor(0.98 * n), 1), n)
  gdal_band_kth(x, band, k, lo, hi, n, block_rows)
}
