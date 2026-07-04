# R reference for the bilateral filter (all edges, both NA policies).
ref_bilateral <- function(x, sigma_d, sigma_r, window, edge = "shrink",
                          edge_value = 0, na_policy = "omit") {
  w <- rep(as.integer(window), length.out = 2L)
  h <- (w - 1L) %/% 2L
  out <- matrix(NA_real_, nrow(x), ncol(x))
  for (cc in seq_len(ncol(x))) {
    for (rr in seq_len(nrow(x))) {
      x0 <- x[rr, cc]
      if (is.na(x0)) next
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
          wgt <- exp(-(dr^2 + dc^2) / (2 * sigma_d^2)) *
            exp(-(v - x0)^2 / (2 * sigma_r^2))
          num <- num + wgt * v
          den <- den + wgt
        }
      }
      out[rr, cc] <- if (bad || den <= 0) NA_real_ else num / den
    }
  }
  out
}

# R reference for the guided filter via the naive focal-mean oracle.
ref_guided <- function(x, guide = x, window, eps, edge = "shrink",
                       na_policy = "omit") {
  both_na <- is.na(x) | is.na(guide)
  gm <- guide
  pm <- x
  gm[both_na] <- NA
  pm[both_na] <- NA
  bm <- function(v) ref_focal(v, window, "mean", edge = edge, na_policy = na_policy)
  m_i <- bm(gm)
  m_p <- bm(pm)
  m_ii <- bm(gm * gm)
  m_ip <- bm(gm * pm)
  a <- (m_ip - m_i * m_p) / (m_ii - m_i^2 + eps)
  b <- m_p - a * m_i
  bm(a) * guide + bm(b)
}

test_that("bilateral matches the reference", {
  m <- ref_mat(nr = 8, nc = 9)
  for (edge in c("shrink", "reflect", "nearest", "constant")) {
    expect_equal(
      rf_bilateral(m, sigma_d = 1, sigma_r = 0.8, window = 5L, edge = edge),
      ref_bilateral(m, 1, 0.8, 5L, edge = edge),
      info = edge
    )
  }
})

test_that("bilateral handles NAs in both policies", {
  m <- ref_mat(nr = 8, nc = 9, na = 7)
  for (pol in c("omit", "propagate")) {
    expect_equal(
      rf_bilateral(m, sigma_d = 1, sigma_r = 0.8, window = 3L, na_policy = pol),
      ref_bilateral(m, 1, 0.8, 3L, na_policy = pol),
      info = pol
    )
  }
})

test_that("bilateral preserves a step edge that rf_mean blurs", {
  truth <- matrix(rep(c(0, 10), each = 40), 8)
  withr::with_seed(3, m <- truth + rnorm(80, sd = 0.5))
  bi <- rf_bilateral(m, sigma_d = 1.5, sigma_r = 2, window = 5L)
  bl <- rf_mean(m, window = 5L)
  # Contrast across the step survives the bilateral, not the boxcar.
  expect_gt(mean(bi[, 6] - bi[, 5]), 8)
  expect_lt(mean(bl[, 6] - bl[, 5]), 5)
  expect_lt(sd(bi - truth), sd(m - truth))
})

test_that("bilateral defaults are derived and validated", {
  m <- ref_mat()
  expect_equal(
    rf_bilateral(m),
    rf_bilateral(m, sigma_r = sd(m), window = 2L * ceiling(2 * 1.5) + 1L)
  )
  expect_error(rf_bilateral(m, sigma_d = 0), "positive")
  expect_error(rf_bilateral(matrix(1, 3, 3)), "supply it explicitly")
})

test_that("guided filter matches the reference, self- and cross-guided", {
  m <- ref_mat(nr = 9, nc = 10)
  g <- ref_mat(nr = 9, nc = 10, seed = 99)
  expect_equal(
    rf_guided(m, window = 3L, eps = 0.04),
    ref_guided(m, window = 3L, eps = 0.04)
  )
  expect_equal(
    rf_guided(m, guide = g, window = 5L, eps = 0.1),
    ref_guided(m, guide = g, window = 5L, eps = 0.1)
  )
})

test_that("guided filter handles NAs and preserves NA centres", {
  m <- ref_mat(nr = 9, nc = 10, na = 8)
  expect_equal(
    rf_guided(m, window = 3L, eps = 0.04),
    ref_guided(m, window = 3L, eps = 0.04)
  )
  out <- rf_guided(m, window = 3L, eps = 0.04)
  expect_true(all(is.na(out[is.na(m)])))
})

test_that("guided filter validates inputs", {
  m <- ref_mat(nr = 5, nc = 5)
  expect_error(rf_guided(m, guide = matrix(1, 2, 2)), "same dimensions")
  expect_error(rf_guided(m, eps = -1), "positive")
  expect_error(rf_guided(matrix(1, 4, 4)), "supply it explicitly")
})

test_that("large eps degenerates the guided filter towards double smoothing", {
  m <- ref_mat(nr = 6, nc = 6)
  out <- rf_guided(m, window = 3L, eps = 1e12)
  twice <- rf_mean(rf_mean(m, 3L), 3L)
  expect_equal(out, twice, tolerance = 1e-6)
})
