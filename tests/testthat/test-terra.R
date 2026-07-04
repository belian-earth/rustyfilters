skip_if_not_installed("terra")

demo_rast <- function(nlyr = 1, na = 0) {
  withr::with_seed(11, {
    a <- array(stats::rnorm(10 * 12 * nlyr), c(10, 12, nlyr))
    if (na > 0) a[sample.int(length(a), na)] <- NA
    r <- terra::rast(
      a,
      crs = "EPSG:32630",
      extent = terra::ext(500000, 501200, 4600000, 4601000)
    )
    names(r) <- paste0("band", seq_len(nlyr))
    r
  })
}

test_that("SpatRaster round-trip preserves geometry, crs and names", {
  r <- demo_rast(nlyr = 2, na = 6)
  out <- rf_lee(r, window = 5L)
  expect_s4_class(out, "SpatRaster")
  expect_equal(as.vector(terra::ext(out)), as.vector(terra::ext(r)))
  expect_true(terra::same.crs(out, r))
  expect_identical(names(out), names(r))
  expect_equal(dim(out), dim(r))
})

test_that("SpatRaster values match the array method", {
  r <- demo_rast(nlyr = 2, na = 6)
  a <- terra::as.array(r)
  for (fn in list(
    function(x, ...) rf_mean(x, 3L, ...),
    function(x, ...) rf_focal(x, 3L, "sd", ...),
    function(x, ...) rf_gaussian(x, sigma = 1, ...),
    function(x, ...) rf_lee(x, 5L, ...)
  )) {
    expect_equal(terra::as.array(fn(r)), fn(a))
  }
})

test_that("round-2 filters have working SpatRaster methods", {
  r <- demo_rast(nlyr = 2, na = 4)
  a <- terra::as.array(r)
  k <- matrix(c(0, 1, 0, 1, -4, 1, 0, 1, 0), 3)
  expect_equal(terra::as.array(rf_convolve(r, k)), rf_convolve(a, k))
  expect_equal(terra::as.array(rf_sobel(r)), rf_sobel(a))
  expect_equal(terra::as.array(rf_laplacian(r)), rf_laplacian(a))
  expect_equal(
    terra::as.array(rf_bilateral(r, sigma_d = 1, sigma_r = 1, window = 3L)),
    rf_bilateral(a, sigma_d = 1, sigma_r = 1, window = 3L)
  )
  expect_equal(
    terra::as.array(rf_lee_sigma_improved(r, 5L)),
    rf_lee_sigma_improved(a, 5L)
  )
})

test_that("rf_guided accepts a SpatRaster guide", {
  r <- demo_rast(nlyr = 1)
  g <- demo_rast(nlyr = 1)
  g <- g * 2
  out <- rf_guided(r, guide = g, window = 3L, eps = 0.1)
  expect_equal(
    terra::as.array(out),
    rf_guided(terra::as.array(r),
      guide = terra::as.array(g),
      window = 3L, eps = 0.1
    )
  )
})

test_that("filter arguments pass through to the SpatRaster method", {
  r <- demo_rast()
  a <- terra::as.array(r)
  expect_equal(
    terra::as.array(rf_median(r, window = 5L, edge = "nearest")),
    rf_median(a, window = 5L, edge = "nearest")
  )
  expect_error(rf_lee(r, looks = -1), "positive")
})
