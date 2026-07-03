test_that("all edge policies match the reference", {
  m <- ref_mat(nr = 5, nc = 6, na = 4)
  for (edge in c("shrink", "reflect", "nearest", "constant")) {
    expect_equal(
      rf_focal(m, window = 3L, stat = "mean", edge = edge),
      ref_focal(m, window = 3L, stat = "mean", edge = edge),
      info = edge
    )
  }
})

test_that("shrink truncates the window at the corner", {
  m <- matrix(as.numeric(1:9), 3, 3)
  out <- rf_focal(m, window = 3L, stat = "sum", edge = "shrink")
  # Corner window holds cells (1,1), (2,1), (1,2), (2,2) = 1 + 2 + 4 + 5.
  expect_equal(out[1, 1], 12)
})

test_that("constant pads with edge_value", {
  m <- matrix(1, 2, 2)
  out <- rf_focal(m, window = 3L, stat = "sum", edge = "constant", edge_value = 10)
  # Each window: 4 ones + 5 tens.
  expect_equal(out, matrix(54, 2, 2))
})

test_that("constant with an NA pad follows the NA policy", {
  m <- matrix(1, 2, 2)
  omit <- rf_focal(m,
    window = 3L, stat = "mean",
    edge = "constant", edge_value = NA
  )
  expect_equal(omit, matrix(1, 2, 2))
  prop <- rf_focal(m,
    window = 3L, stat = "mean",
    edge = "constant", edge_value = NA, na_policy = "propagate"
  )
  expect_true(all(is.na(prop)))
})

test_that("nearest repeats the closest edge cell", {
  m <- matrix(as.numeric(1:4), 2, 2)
  out <- rf_focal(m, window = 3L, stat = "min", edge = "nearest")
  expect_equal(out, matrix(1, 2, 2))
})

test_that("an all-NA window yields NA under omit", {
  m <- matrix(NA_real_, 3, 3)
  m[1, 1] <- 5
  out <- rf_focal(m, window = 3L, stat = "mean")
  expect_equal(out[3, 3], NA_real_)
  expect_equal(out[1, 1], 5)
})
