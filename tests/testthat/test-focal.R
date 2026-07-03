stats_all <- c("mean", "median", "min", "max", "range", "sd", "sum", "mode")

test_that("every stat matches the reference on a clean matrix", {
  m <- ref_mat()
  for (stat in stats_all) {
    expect_equal(
      rf_focal(m, window = 3L, stat = stat),
      ref_focal(m, window = 3L, stat = stat),
      info = stat
    )
  }
})

test_that("every stat matches the reference with NAs, omit policy", {
  m <- ref_mat(na = 12)
  for (stat in stats_all) {
    expect_equal(
      rf_focal(m, window = 3L, stat = stat),
      ref_focal(m, window = 3L, stat = stat),
      info = stat
    )
  }
})

test_that("every stat matches the reference with NAs, propagate policy", {
  m <- ref_mat(na = 12)
  for (stat in stats_all) {
    expect_equal(
      rf_focal(m, window = 3L, stat = stat, na_policy = "propagate"),
      ref_focal(m, window = 3L, stat = stat, na_policy = "propagate"),
      info = stat
    )
  }
})

test_that("rectangular windows match the reference", {
  m <- ref_mat(na = 5)
  for (stat in c("mean", "median", "sd")) {
    expect_equal(
      rf_focal(m, window = c(3L, 5L), stat = stat),
      ref_focal(m, window = c(3L, 5L), stat = stat),
      info = stat
    )
    expect_equal(
      rf_focal(m, window = c(5L, 3L), stat = stat),
      ref_focal(m, window = c(5L, 3L), stat = stat),
      info = stat
    )
  }
})

test_that("windows larger than the matrix work", {
  m <- ref_mat(nr = 3, nc = 3)
  expect_equal(
    rf_focal(m, window = 9L, stat = "mean"),
    ref_focal(m, window = 9L, stat = "mean")
  )
  expect_equal(
    rf_focal(m, window = 9L, stat = "mean", edge = "reflect"),
    ref_focal(m, window = 9L, stat = "mean", edge = "reflect")
  )
})

test_that("mode breaks ties towards the lowest value", {
  m <- matrix(c(2, 1, 1, 2, 5, 5, 3, 4, 6), 3, 3)
  out <- rf_focal(m, window = 3L, stat = "mode")
  expect_equal(out[2, 2], 1)
})

test_that("window = 1 is an identity for location stats", {
  m <- ref_mat()
  expect_equal(rf_focal(m, window = 1L, stat = "mean"), m)
  expect_equal(rf_focal(m, window = 1L, stat = "median"), m)
})

test_that("dimnames are preserved", {
  m <- ref_mat(nr = 3, nc = 3)
  dimnames(m) <- list(letters[1:3], LETTERS[1:3])
  expect_identical(dimnames(rf_focal(m)), dimnames(m))
})

test_that("integer matrices are accepted and return doubles", {
  m <- matrix(1:12, 3, 4)
  out <- rf_focal(m, stat = "mean")
  expect_type(out, "double")
  expect_equal(out, ref_focal(matrix(as.numeric(1:12), 3, 4)))
})

test_that("invalid inputs error informatively", {
  m <- ref_mat(nr = 3, nc = 3)
  expect_error(rf_focal(1:10), "matrix or 3-D array")
  expect_error(rf_focal(letters), "matrix or 3-D array")
  expect_error(rf_focal(array(1, c(2, 2, 2, 2))), "matrix or 3-D array")
  expect_error(rf_focal(m, window = 4L), "odd positive")
  expect_error(rf_focal(m, window = -3L), "odd positive")
  expect_error(rf_focal(m, window = c(3L, 3L, 3L)), "odd positive")
  expect_error(rf_focal(m, stat = "banana"), "must be one of")
  expect_error(rf_focal(m, edge = "banana"), "must be one of")
  expect_error(rf_focal(m, edge_value = c(1, 2)), "single number")
})
