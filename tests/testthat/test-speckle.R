# Pure-R references for the speckle filters (shrink edge, omit policy),
# built on the helper-reference window gather. Population variance
# throughout, matching the filters' derivations.

ref_local <- function(vals) {
  vals <- vals[!is.na(vals)]
  m <- mean(vals)
  list(m = m, v = mean((vals - m)^2), n = length(vals))
}

ref_speckle <- function(x, window, fun) {
  w <- rep(as.integer(window), length.out = 2L)
  out <- matrix(NA_real_, nrow(x), ncol(x))
  for (cc in seq_len(ncol(x))) {
    for (rr in seq_len(nrow(x))) {
      x0 <- x[rr, cc]
      if (is.na(x0)) next
      vals <- ref_window(x, rr, cc, w, "shrink", 0)
      out[rr, cc] <- fun(x0, ref_local(vals))
    }
  }
  out
}

ref_lee <- function(x, window, looks) {
  cu2 <- 1 / looks
  ref_speckle(x, window, function(x0, s) {
    noise <- s$m^2 * cu2
    if (s$v + noise <= 0) return(s$m)
    s$m + (s$v / (s$v + noise)) * (x0 - s$m)
  })
}

ref_kuan <- function(x, window, looks) {
  cu2 <- 1 / looks
  ref_speckle(x, window, function(x0, s) {
    if (s$m == 0 || s$v == 0) return(s$m)
    ci2 <- s$v / s$m^2
    w <- min(max((1 - cu2 / ci2) / (1 + cu2), 0), 1)
    s$m + w * (x0 - s$m)
  })
}

ref_enhanced_lee <- function(x, window, looks, damping) {
  cu <- sqrt(1 / looks)
  cmax <- sqrt(1 + 2 / looks)
  ref_speckle(x, window, function(x0, s) {
    if (s$m == 0) return(s$m)
    ci <- sqrt(s$v) / s$m
    if (ci <= cu) return(s$m)
    if (ci >= cmax) return(x0)
    w <- exp(-damping * (ci - cu) / (cmax - ci))
    s$m * w + x0 * (1 - w)
  })
}

ref_gamma_map <- function(x, window, looks) {
  cu2 <- 1 / looks
  cu <- sqrt(cu2)
  cmax <- sqrt(1 + 2 / looks)
  ref_speckle(x, window, function(x0, s) {
    if (s$m == 0) return(s$m)
    ci2 <- s$v / s$m^2
    ci <- sqrt(ci2)
    if (ci <= cu) return(s$m)
    if (ci >= cmax) return(x0)
    alpha <- (1 + cu2) / (ci2 - cu2)
    b <- alpha - looks - 1
    dsc <- s$m^2 * b^2 + 4 * alpha * looks * s$m * x0
    (b * s$m + sqrt(max(dsc, 0))) / (2 * alpha)
  })
}

# Lee sigma and Frost need the window values / offsets, not just moments.
ref_lee_sigma <- function(x, window, looks, k, min_count) {
  w <- rep(as.integer(window), length.out = 2L)
  cu <- sqrt(1 / looks)
  out <- matrix(NA_real_, nrow(x), ncol(x))
  for (cc in seq_len(ncol(x))) {
    for (rr in seq_len(nrow(x))) {
      x0 <- x[rr, cc]
      if (is.na(x0)) next
      vals <- ref_window(x, rr, cc, w, "shrink", 0)
      vals <- vals[!is.na(vals)]
      bounds <- sort(x0 * c(1 - k * cu, 1 + k * cu))
      inr <- vals[vals >= bounds[1] & vals <= bounds[2]]
      out[rr, cc] <- if (length(inr) < min_count) mean(vals) else mean(inr)
    }
  }
  out
}

ref_frost <- function(x, window, damping) {
  w <- rep(as.integer(window), length.out = 2L)
  h <- (w - 1L) %/% 2L
  out <- matrix(NA_real_, nrow(x), ncol(x))
  for (cc in seq_len(ncol(x))) {
    for (rr in seq_len(nrow(x))) {
      x0 <- x[rr, cc]
      if (is.na(x0)) next
      s <- ref_local(ref_window(x, rr, cc, w, "shrink", 0))
      b <- if (s$m == 0) 0 else damping * s$v / s$m^2
      num <- 0
      den <- 0
      for (dc in -h[2]:h[2]) {
        for (dr in -h[1]:h[1]) {
          ri <- rr + dr
          ci <- cc + dc
          if (ri < 1 || ri > nrow(x) || ci < 1 || ci > ncol(x)) next
          v <- x[ri, ci]
          if (is.na(v)) next
          wgt <- exp(-b * sqrt(dr^2 + dc^2))
          num <- num + wgt * v
          den <- den + wgt
        }
      }
      out[rr, cc] <- num / den
    }
  }
  out
}

# A reproducible speckled scene: two-level truth times exponential speckle.
speckled_scene <- function(nr = 12, nc = 14, na = 0) {
  withr::with_seed(7, {
    n <- nr * nc
    truth <- matrix(rep(c(1, 5), times = c(ceiling(n / 2), floor(n / 2))), nr, nc)
    m <- truth * stats::rexp(n)
    if (na > 0) m[sample.int(length(m), na)] <- NA_real_
    m
  })
}

test_that("all speckle filters match their references on a clean scene", {
  m <- speckled_scene()
  expect_equal(rf_lee(m, 5L, looks = 1), ref_lee(m, 5L, 1))
  expect_equal(rf_kuan(m, 5L, looks = 1), ref_kuan(m, 5L, 1))
  expect_equal(
    rf_enhanced_lee(m, 5L, looks = 1, damping = 1),
    ref_enhanced_lee(m, 5L, 1, 1)
  )
  expect_equal(rf_gamma_map(m, 5L, looks = 3), ref_gamma_map(m, 5L, 3))
  expect_equal(
    rf_lee_sigma(m, 5L, looks = 1, k = 2, min_count = 3L),
    ref_lee_sigma(m, 5L, 1, 2, 3)
  )
  expect_equal(rf_frost(m, 5L, damping = 2), ref_frost(m, 5L, 2))
})

test_that("speckle filters match their references with NAs", {
  m <- speckled_scene(na = 15)
  expect_equal(rf_lee(m, 5L, looks = 2), ref_lee(m, 5L, 2))
  expect_equal(rf_kuan(m, 5L, looks = 2), ref_kuan(m, 5L, 2))
  expect_equal(
    rf_enhanced_lee(m, 5L, looks = 2, damping = 1.5),
    ref_enhanced_lee(m, 5L, 2, 1.5)
  )
  expect_equal(rf_gamma_map(m, 5L, looks = 2), ref_gamma_map(m, 5L, 2))
  expect_equal(
    rf_lee_sigma(m, 5L, looks = 2, k = 2, min_count = 3L),
    ref_lee_sigma(m, 5L, 2, 2, 3)
  )
  expect_equal(rf_frost(m, 5L, damping = 2), ref_frost(m, 5L, 2))
})

test_that("a constant window passes through every filter unchanged", {
  m <- matrix(4, 9, 9)
  for (fn in list(rf_lee, rf_kuan, rf_enhanced_lee, rf_gamma_map, rf_frost, rf_lee_sigma)) {
    expect_equal(fn(m, window = 3L), m)
  }
})

test_that("point targets are preserved by enhanced Lee and Gamma-MAP", {
  m <- matrix(0.1, 9, 9)
  m[5, 5] <- 100
  expect_equal(rf_enhanced_lee(m, 7L, looks = 1)[5, 5], 100)
  expect_equal(rf_gamma_map(m, 7L, looks = 1)[5, 5], 100)
})

test_that("NA centres stay NA; neighbours of NAs are still filtered", {
  m <- speckled_scene()
  m[6, 7] <- NA
  for (fn in list(rf_lee, rf_kuan, rf_enhanced_lee, rf_gamma_map, rf_frost, rf_lee_sigma)) {
    out <- fn(m, window = 5L)
    expect_true(is.na(out[6, 7]))
    expect_false(anyNA(out[-6, -7]))
  }
})

test_that("propagate policy spreads NA over the whole window", {
  m <- speckled_scene()
  m[6, 7] <- NA
  out <- rf_lee(m, window = 5L, na_policy = "propagate")
  expect_true(all(is.na(out[4:8, 5:9])))
  expect_false(is.na(out[1, 1]))
})

test_that("speckle parameters are validated", {
  m <- speckled_scene(nr = 5, nc = 5)
  expect_error(rf_lee(m, looks = 0), "positive")
  expect_error(rf_lee(m, looks = c(1, 2)), "positive")
  expect_error(rf_enhanced_lee(m, damping = -1), "positive")
  expect_error(rf_frost(m, damping = 0), "positive")
  expect_error(rf_lee_sigma(m, k = 0), "positive")
  expect_error(rf_lee_sigma(m, min_count = 0L), "positive integer")
  expect_error(rf_gamma_map("no"), "matrix or 3-D array")
})
