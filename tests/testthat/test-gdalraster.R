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

test_that("block streaming matches whole-raster results", {
  ds <- demo_gdal(nbands = 2, na = 60)
  withr::defer(ds$close())
  # Awkward block heights: not dividing nrows, and smaller than the halo.
  cases <- list(
    list(function(d, ...) rf_lee(d, window = 7L, ...), tol = 1e-12),
    list(function(d, ...) rf_median(d, window = 5L, ...), tol = 0),
    list(function(d, ...) rf_gaussian(d, sigma = 1.5, ...), tol = 0),
    list(function(d, ...) rf_guided(d, window = 5L, eps = 0.1, ...), tol = 1e-12),
    list(function(d, ...) rf_frost(d, window = 5L, ...), tol = 1e-12),
    list(function(d, ...) rf_sobel(d, ...), tol = 0)
  )
  for (case in cases) {
    fn <- case[[1]]
    whole <- fn(ds)
    for (rows in c(13L, 3L)) {
      blocked <- fn(ds, by_block = TRUE, block_rows = rows)
      if (case$tol == 0) {
        expect_identical(gdal_as_array(blocked), gdal_as_array(whole))
      } else {
        expect_equal(gdal_as_array(blocked), gdal_as_array(whole),
          tolerance = case$tol
        )
      }
      blocked$close()
    }
    whole$close()
  }
})

test_that("blocked improved Lee sigma is exactly the whole-raster result", {
  ds <- demo_gdal(nbands = 2, na = 30)
  withr::defer(ds$close())
  whole <- rf_lee_sigma_improved(ds, window = 7L)
  withr::defer(whole$close())
  blocked <- rf_lee_sigma_improved(ds, window = 7L, by_block = TRUE, block_rows = 17L)
  withr::defer(blocked$close())
  expect_identical(gdal_as_array(blocked), gdal_as_array(whole))
})

test_that("the streaming z98 helper is exact", {
  ds <- demo_gdal(na = 25)
  withr::defer(ds$close())
  v <- as.numeric(gdalraster::read_ds(ds))
  v <- sort(v[!is.na(v)])
  k <- min(max(floor(0.98 * length(v)), 1), length(v))
  expect_identical(gdal_band_z98(ds, 1L, block_rows = 11L), v[k])
})

test_that("streaming defaults to a GeoTIFF tempfile on disk", {
  ds <- demo_gdal()
  withr::defer(ds$close())
  out <- rf_mean(ds, window = 3L, by_block = TRUE, block_rows = 20L)
  withr::defer(out$close())
  f <- out$getFilename()
  expect_false(startsWith(f, "/vsimem"))
  expect_true(file.exists(f))
  expect_match(f, "\\.tif$")
})

test_that("the memory budget option triggers streaming automatically", {
  ds <- demo_gdal()
  withr::defer(ds$close())
  withr::local_options(rustyfilters.block_memory = 1024)
  out <- rf_mean(ds, window = 3L)
  withr::defer(out$close())
  # Auto-streamed: output went to a disk tempfile, not /vsimem.
  expect_false(startsWith(out$getFilename(), "/vsimem"))
  withr::local_options(rustyfilters.block_memory = NULL)
  out2 <- rf_mean(ds, window = 3L)
  withr::defer(out2$close())
  expect_true(startsWith(out2$getFilename(), "/vsimem"))
  expect_equal(gdal_as_array(out), gdal_as_array(out2), tolerance = 1e-12)
})

test_that("cross-guidance streams the guide alongside the input", {
  ds <- demo_gdal()
  g <- demo_gdal()
  withr::defer({
    ds$close()
    g$close()
  })
  whole <- rf_guided(ds, guide = g, window = 3L, eps = 0.1)
  withr::defer(whole$close())
  blocked <- rf_guided(ds,
    guide = g, window = 3L, eps = 0.1,
    by_block = TRUE, block_rows = 15L
  )
  withr::defer(blocked$close())
  expect_equal(gdal_as_array(blocked), gdal_as_array(whole), tolerance = 1e-12)
  # An in-memory array guide is row-sliced per block the same way.
  blocked2 <- rf_guided(ds,
    guide = gdal_as_array(g), window = 3L, eps = 0.1,
    by_block = TRUE, block_rows = 15L
  )
  withr::defer(blocked2$close())
  expect_equal(gdal_as_array(blocked2), gdal_as_array(whole), tolerance = 1e-12)
})

test_that("streaming arguments are validated", {
  ds <- demo_gdal()
  withr::defer(ds$close())
  expect_error(rf_mean(ds, by_block = "yes"), "TRUE")
  expect_error(rf_mean(ds, by_block = TRUE, block_rows = 0L), "positive integer")
})
