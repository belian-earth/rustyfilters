skip_if_not_installed("gdalraster")

# A small /vsimem dataset with known values and NA holes, built from the
# gdalraster sample raster's geometry.
demo_gdal <- function(nbands = 1, na = 0) {
  src <- system.file("extdata/storml_elev.tif", package = "gdalraster")
  f <- paste0("/vsimem/", basename(tempfile("rf_test_")), ".tif")
  gdalraster::rasterFromRaster(src, f, fmt = "GTiff", dtName = "Float64",
                               nbands = nbands)
  ds <- methods::new(gdalraster::GDALRaster, f, read_only = FALSE)
  dm <- ds$dim()
  withr::with_seed(21, {
    for (b in seq_len(nbands)) {
      v <- stats::rnorm(dm[1] * dm[2])
      if (na > 0) v[sample.int(length(v), na)] <- NA
      ds$setNoDataValue(b, NaN)
      ds$write(
        band = b, xoff = 0L, yoff = 0L, xsize = dm[1], ysize = dm[2],
        rasterData = v
      )
    }
  })
  ds
}

gdal_as_array <- function(ds) {
  aperm(array(as.numeric(gdalraster::read_ds(ds)), dim = ds$dim()), c(2, 1, 3))
}

test_that("GDALRaster round-trip preserves geometry, srs and values", {
  ds <- demo_gdal(nbands = 2, na = 25)
  withr::defer(ds$close())
  a <- gdal_as_array(ds)
  out <- rf_lee(ds, window = 5L)
  withr::defer(out$close())
  expect_s4_class(out, "Rcpp_GDALRaster")
  expect_identical(out$dim(), ds$dim())
  expect_identical(out$bbox(), ds$bbox())
  expect_identical(out$getProjection(), ds$getProjection())
  expect_equal(gdal_as_array(out), rf_lee(a, window = 5L))
})

test_that("NA cells survive the GDALRaster round-trip", {
  ds <- demo_gdal(na = 40)
  withr::defer(ds$close())
  a <- gdal_as_array(ds)
  out <- rf_lee(ds, window = 5L)
  withr::defer(out$close())
  expect_identical(is.na(gdal_as_array(out)), is.na(a))
})

test_that("a sample of filters agree with the array method", {
  ds <- demo_gdal(na = 10)
  withr::defer(ds$close())
  a <- gdal_as_array(ds)
  k <- matrix(c(0, 1, 0, 1, -4, 1, 0, 1, 0), 3)
  cases <- list(
    list(function(x, ...) rf_median(x, 3L, edge = "nearest", ...)),
    list(function(x, ...) rf_gaussian(x, sigma = 1, ...)),
    list(function(x, ...) rf_convolve(x, k, ...)),
    list(function(x, ...) rf_sobel(x, ...)),
    list(function(x, ...) rf_bilateral(x, sigma_d = 1, sigma_r = 1, window = 3L, ...)),
    list(function(x, ...) rf_lee_sigma_improved(x, 5L, ...)),
    list(function(x, ...) rf_focal(x, 3L, "sd", ...))
  )
  for (case in cases) {
    fn <- case[[1]]
    out <- fn(ds)
    expect_equal(gdal_as_array(out), fn(a))
    out$close()
  }
})

test_that("rf_guided accepts a GDALRaster guide", {
  ds <- demo_gdal()
  g <- demo_gdal()
  withr::defer({
    ds$close()
    g$close()
  })
  out <- rf_guided(ds, guide = g, window = 3L, eps = 0.1)
  withr::defer(out$close())
  expect_equal(
    gdal_as_array(out),
    rf_guided(gdal_as_array(ds),
      guide = gdal_as_array(g),
      window = 3L, eps = 0.1
    )
  )
})

test_that("filename argument writes the result to the given path", {
  ds <- demo_gdal()
  withr::defer(ds$close())
  f <- file.path(withr::local_tempdir(), "filtered.tif")
  out <- rf_mean(ds, window = 3L, filename = f)
  withr::defer(out$close())
  expect_true(file.exists(f))
  expect_identical(out$getFilename(), f)
})

test_that("filter arguments are validated for GDALRaster input", {
  ds <- demo_gdal()
  withr::defer(ds$close())
  expect_error(rf_lee(ds, looks = -1), "positive")
})
