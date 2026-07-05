test_that("s1_sthelens is valid single-look intensity data", {
  expect_true(is.matrix(s1_sthelens))
  expect_identical(dim(s1_sthelens), c(300L, 300L))
  expect_false(anyNA(s1_sthelens))
  expect_true(all(s1_sthelens > 0))
})

test_that("rf_plot draws matrices and array layers invisibly", {
  grDevices::pdf(NULL)
  withr::defer(grDevices::dev.off())
  m <- ref_mat()
  expect_invisible(rf_plot(m))
  expect_identical(withVisible(rf_plot(m))$value, m)
  a <- array(stats::rnorm(60), c(5, 4, 3))
  expect_no_error(rf_plot(a, layer = 2L))
  expect_no_error(rf_plot(matrix(1, 3, 3))) # constant field
  expect_no_error(rf_plot(volcano, palette = "Viridis", stretch = c(0, 1)))
})

test_that("rf_plot validates its inputs", {
  expect_error(rf_plot(1:10), "matrix or 3-D array")
  expect_error(rf_plot(ref_mat(), stretch = c(-1, 2)), "probabilities")
  expect_error(rf_plot(ref_mat(), stretch = 0.5), "probabilities")
})
