test_that("set/get threads works and returns the previous count", {
  withr::defer(rf_set_threads(rf_get_threads()))
  rf_set_threads(1L)
  old <- rf_set_threads(2L)
  expect_equal(old, 1L)
  expect_equal(rf_get_threads(), 2L)
})

test_that("thread count is validated", {
  expect_error(rf_set_threads(0L), "positive integer")
  expect_error(rf_set_threads(NA_integer_), "positive integer")
  expect_error(rf_set_threads(c(1L, 2L)))
})

test_that("default_threads honours option, env var and check limits", {
  withr::local_options(rustyfilters.threads = 3L)
  expect_equal(default_threads(), 3L)

  withr::local_options(rustyfilters.threads = NULL)
  withr::local_envvar(RUSTYFILTERS_NUM_THREADS = "5")
  expect_equal(default_threads(), 5L)

  withr::local_envvar(
    RUSTYFILTERS_NUM_THREADS = NA,
    "_R_CHECK_LIMIT_CORES_" = "TRUE"
  )
  expect_equal(default_threads(), 2L)

  withr::local_envvar("_R_CHECK_LIMIT_CORES_" = NA)
  expect_equal(default_threads(), 0L)
})
