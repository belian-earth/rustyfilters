# Build the `s1_sthelens` package dataset: a 300 x 300 patch of Sentinel-1A
# radiometrically terrain corrected (RTC) gamma-nought backscatter (VV
# polarisation, linear power, 10 m pixels) centred on the crater of Mount
# St Helens, USA (lon -122.19, lat 46.20; UTM zone 10N).
#
# Source scene: S1A_IW_GRDH_1SDV_20240928T020214_20240928T020239_055859_06D3B3
# (acquired 2024-09-28, ascending), sentinel-1-rtc collection on Microsoft
# Planetary Computer. Contains modified Copernicus Sentinel data (2024).
#
# Values are rounded to 4 significant digits to aid compression; this is far
# below the radiometric accuracy of the product.

library(gdalraster)

search <- jsonlite::fromJSON(system(paste(
  "curl -s -X POST",
  "'https://planetarycomputer.microsoft.com/api/stac/v1/search'",
  "-H 'Content-Type: application/json'",
  "-d '{\"collections\":[\"sentinel-1-rtc\"],",
  "\"intersects\":{\"type\":\"Point\",\"coordinates\":[-122.19,46.20]},",
  "\"datetime\":\"2024-09-28T00:00:00Z/2024-09-28T23:59:59Z\",\"limit\":1}'"
), intern = TRUE), simplifyVector = FALSE)
href <- search$features[[1]]$assets$vv$href
token <- jsonlite::fromJSON(
  "https://planetarycomputer.microsoft.com/api/sas/v1/token/sentinel-1-rtc"
)$token

ds <- new(GDALRaster, paste0("/vsicurl/", href, "?", token))
gt <- ds$getGeoTransform()
pt <- transform_xy(matrix(c(-122.19, 46.20), 1),
  srs_from = "EPSG:4326", srs_to = ds$getProjection()
)
px <- floor((pt[1] - gt[1]) / gt[2])
py <- floor((pt[2] - gt[4]) / gt[6])
n <- 300L
v <- ds$read(
  band = 1, xoff = px - n / 2, yoff = py - n / 2,
  xsize = n, ysize = n, out_xsize = n, out_ysize = n
)
ds$close()

# Row-major GDAL read -> matrix with row 1 = north.
s1_sthelens <- signif(matrix(as.numeric(v), ncol = n, byrow = TRUE), 4)
stopifnot(!anyNA(s1_sthelens), all(s1_sthelens > 0))

usethis::use_data(s1_sthelens, overwrite = TRUE, compress = "xz")
