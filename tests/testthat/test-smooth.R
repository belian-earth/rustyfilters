# Direct 2-D normalized-convolution Gaussian reference.
ref_gaussian <- function(x, sigma, window, edge = "shrink", edge_value = 0,
                         na_policy = "omit") {
  w <- rep(as.integer(window), length.out = 2L)
  h <- (w - 1L) %/% 2L
  k1 <- function(n) {
    i <- seq(-(n - 1) / 2, (n - 1) / 2)
    v <- exp(-i^2 / (2 * sigma^2))
    v / sum(v)
  }
  kk <- outer(k1(w[1]), k1(w[2]))
  out <- matrix(NA_real_, nrow(x), ncol(x), dimnames = dimnames(x))
  for (cc in seq_len(ncol(x))) {
    for (rr in seq_len(nrow(x))) {
      num <- 0
      den <- 0
      bad <- FALSE
      for (dc in -h[2]:h[2]) {
        for (dr in -h[1]:h[1]) {
          ri <- ref_map_idx(rr + dr, nrow(x), edge)
          ci <- ref_map_idx(cc + dc, ncol(x), edge)
          outside <- is.na(ri) || is.na(ci)
          if (outside && edge != "constant") next
          v <- if (outside) edge_value else x[ri, ci]
          if (is.na(v)) {
            if (na_policy == "propagate") bad <- TRUE
            next
          }
          kw <- kk[dr + h[1] + 1, dc + h[2] + 1]
          num <- num + kw * v
          den <- den + kw
        }
      }
      out[rr, cc] <- if (bad || den <= 1e-12) NA_real_ else num / den
    }
  }
  out
}

test_that("gaussian matches the 2-D reference on a clean matrix", {
  m <- ref_mat()
  for (edge in c("shrink", "reflect", "nearest", "constant")) {
    expect_equal(
      rf_gaussian(m, sigma = 1, window = 5L, edge = edge, edge_value = 2),
      ref_gaussian(m, sigma = 1, window = 5L, edge = edge, edge_value = 2),
      info = edge
    )
  }
})

test_that("gaussian matches the 2-D reference with NAs in both policies", {
  m <- ref_mat(na = 10)
  for (edge in c("shrink", "reflect", "nearest", "constant")) {
    for (pol in c("omit", "propagate")) {
      expect_equal(
        rf_gaussian(m, sigma = 1.3, window = 3L, edge = edge, na_policy = pol),
        ref_gaussian(m, sigma = 1.3, window = 3L, edge = edge, na_policy = pol),
        info = paste(edge, pol)
      )
    }
  }
})

test_that("gaussian supports rectangular windows", {
  m <- ref_mat()
  expect_equal(
    rf_gaussian(m, sigma = 2, window = c(3L, 5L)),
    ref_gaussian(m, sigma = 2, window = c(3L, 5L))
  )
})

test_that("gaussian window defaults to 2 * ceiling(3 * sigma) + 1", {
  m <- ref_mat(nr = 15, nc = 15)
  expect_equal(rf_gaussian(m, sigma = 1), rf_gaussian(m, sigma = 1, window = 7L))
  expect_equal(rf_gaussian(m, sigma = 0.5), rf_gaussian(m, sigma = 0.5, window = 5L))
})

test_that("gaussian leaves a constant field unchanged", {
  m <- matrix(3.5, 6, 6)
  expect_equal(rf_gaussian(m, sigma = 1), m)
})

test_that("gaussian validates sigma", {
  m <- ref_mat(nr = 3, nc = 3)
  expect_error(rf_gaussian(m, sigma = 0), "positive")
  expect_error(rf_gaussian(m, sigma = -1), "positive")
  expect_error(rf_gaussian(m, sigma = c(1, 2)), "positive")
})

test_that("rf_mean and rf_median delegate to the focal engine", {
  m <- ref_mat(na = 8)
  expect_equal(rf_mean(m, window = 5L), ref_focal(m, 5L, "mean"))
  expect_equal(rf_median(m, window = 3L), ref_focal(m, 3L, "median"))
  expect_equal(
    rf_median(m, window = 3L, na_policy = "propagate"),
    ref_focal(m, 3L, "median", na_policy = "propagate")
  )
})

test_that("median averages the two middle values in even windows", {
  m <- matrix(as.numeric(1:4), 2, 2)
  out <- rf_median(m, window = 3L)
  # Every shrink window is the full matrix: median(1:4) = 2.5.
  expect_equal(out, matrix(2.5, 2, 2))
})

test_that("rf_mean matches terra::focal", {
  skip_if_not_installed("terra")
  m <- ref_mat(nr = 10, nc = 12)
  ours <- rf_mean(m, window = 3L, na_policy = "propagate")
  r <- terra::rast(m)
  theirs <- terra::as.matrix(
    terra::focal(r, w = 3, fun = "mean", na.rm = FALSE),
    wide = TRUE
  )
  # terra leaves edge cells NA when the window overhangs without expand;
  # compare the interior.
  expect_equal(ours[2:9, 2:11], theirs[2:9, 2:11])
})
