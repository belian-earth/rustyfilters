test_that("a 3-D array filters each layer like the matrix method", {
  a <- array(0, c(7, 9, 3))
  a[, , 1] <- ref_mat(7, 9, na = 5, seed = 1)
  a[, , 2] <- ref_mat(7, 9, na = 0, seed = 2)
  a[, , 3] <- ref_mat(7, 9, na = 9, seed = 3)
  for (fn in list(
    function(x) rf_focal(x, 3L, "sd"),
    function(x) rf_mean(x, 5L),
    function(x) rf_median(x, 3L, edge = "reflect"),
    function(x) rf_gaussian(x, sigma = 1),
    function(x) rf_lee(x, 5L),
    function(x) rf_frost(x, 5L)
  )) {
    out <- fn(a)
    expect_equal(dim(out), dim(a))
    for (l in 1:3) {
      expect_equal(out[, , l], fn(a[, , l]))
    }
  }
})

test_that("array dimnames are preserved", {
  a <- array(as.numeric(1:24), c(2, 3, 4))
  dimnames(a) <- list(c("r1", "r2"), c("c1", "c2", "c3"), paste0("b", 1:4))
  expect_identical(dimnames(rf_mean(a)), dimnames(a))
})

test_that("higher-dimensional arrays are rejected", {
  expect_error(rf_mean(array(1, c(2, 2, 2, 2))), "matrix or 3-D array")
})

test_that("results are identical across thread counts", {
  withr::defer(rf_set_threads(rf_get_threads()))
  a <- array(stats::rnorm(7 * 9 * 2), c(7, 9, 2))
  a[c(3, 40, 90)] <- NA
  rf_set_threads(1L)
  seq_out <- list(
    rf_mean(a, 5L), rf_focal(a, 3L, "median"),
    rf_gaussian(a, sigma = 1), rf_lee(a, 5L), rf_frost(a, 3L)
  )
  rf_set_threads(4L)
  par_out <- list(
    rf_mean(a, 5L), rf_focal(a, 3L, "median"),
    rf_gaussian(a, sigma = 1), rf_lee(a, 5L), rf_frost(a, 3L)
  )
  expect_identical(seq_out, par_out)
})
