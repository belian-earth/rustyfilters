# Direct R reference for kernel filtering (cross-correlation).
ref_convolve <- function(x, kernel, normalize = FALSE, edge = "shrink",
                         edge_value = 0, na_policy = "omit") {
  h <- (dim(kernel) - 1L) %/% 2L
  out <- matrix(NA_real_, nrow(x), ncol(x), dimnames = dimnames(x))
  for (cc in seq_len(ncol(x))) {
    for (rr in seq_len(nrow(x))) {
      num <- 0
      den <- 0
      n <- 0L
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
          kw <- kernel[dr + h[1] + 1, dc + h[2] + 1]
          num <- num + kw * v
          den <- den + kw
          n <- n + 1L
        }
      }
      out[rr, cc] <- if (bad || n == 0L) {
        NA_real_
      } else if (normalize) {
        if (abs(den) > 1e-12) num / den else NA_real_
      } else {
        num
      }
    }
  }
  out
}

test_that("rf_convolve matches the reference for an asymmetric kernel", {
  m <- ref_mat()
  k <- matrix(c(1, 2, 0, -1, 3, 0.5, 0, -2, 1), 3)
  for (edge in c("shrink", "reflect", "nearest", "constant")) {
    expect_equal(
      rf_convolve(m, k, edge = edge, edge_value = 1),
      ref_convolve(m, k, edge = edge, edge_value = 1),
      info = edge
    )
  }
})

test_that("rf_convolve handles NAs in both policies", {
  m <- ref_mat(na = 8)
  k <- matrix(c(0, 1, 0, 1, -4, 1, 0, 1, 0), 3)
  for (pol in c("omit", "propagate")) {
    expect_equal(
      rf_convolve(m, k, na_policy = pol),
      ref_convolve(m, k, na_policy = pol),
      info = pol
    )
  }
})

test_that("rectangular and normalized kernels work", {
  m <- ref_mat(na = 4)
  k <- matrix(1, 3, 5)
  expect_equal(
    rf_convolve(m, k / 15, normalize = TRUE),
    ref_convolve(m, k / 15, normalize = TRUE)
  )
  # A normalized box kernel is the moving mean.
  expect_equal(
    rf_convolve(m, k, normalize = TRUE),
    rf_mean(m, window = c(3L, 5L))
  )
})

test_that("identity kernel returns the input", {
  m <- ref_mat()
  k <- matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0), 3)
  expect_equal(rf_convolve(m, k), m)
})

test_that("kernel validation errors are informative", {
  m <- ref_mat(nr = 3, nc = 3)
  expect_error(rf_convolve(m, matrix(1, 2, 2)), "odd dimensions")
  expect_error(rf_convolve(m, matrix(c(1, NA, 1), 1, 3)), "odd dimensions")
  expect_error(rf_convolve(m, "no"), "odd dimensions")
  expect_error(rf_convolve(m, matrix(1, 3, 3), normalize = "yes"), "TRUE")
})

test_that("sobel detects a vertical step and is zero on flats", {
  m <- matrix(rep(c(1, 1, 5, 5, 5), each = 5), 5)
  gx <- rf_sobel(m, direction = "x")
  expect_equal(gx[, 4:5], matrix(0, 5, 2))
  expect_true(all(gx[, 2:3] == 16))
  expect_equal(rf_sobel(m, direction = "y"), matrix(0, 5, 5))
  expect_equal(rf_sobel(m), abs(gx))
})

test_that("sobel x and y are transposes of one another", {
  m <- ref_mat(nr = 8, nc = 8)
  expect_equal(rf_sobel(m, direction = "x"), t(rf_sobel(t(m), direction = "y")))
})

test_that("laplacian is zero on a linear ramp and flags a spike", {
  m <- matrix(as.numeric(1:25), 5)
  expect_equal(rf_laplacian(m)[2:4, 2:4], matrix(0, 3, 3))
  m[3, 3] <- 100
  expect_true(abs(rf_laplacian(m)[3, 3]) > 100)
  expect_error(rf_laplacian(m, neighbours = 5), "4 or 8")
})
