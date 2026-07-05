#' Sentinel-1 backscatter over Mount St Helens
#'
#' A 300 x 300 numeric matrix of Sentinel-1A radiometrically terrain
#' corrected (RTC) gamma-nought backscatter, VV polarisation, in linear
#' power units: single-look SAR intensity as the speckle filters expect it.
#' The 3 km x 3 km patch (10 m pixels, UTM zone 10N, row 1 northmost)
#' is centred on the crater of Mount St Helens, USA (approximately
#' 46.20 N, 122.19 W), acquired on 2024-09-28.
#'
#' Values are rounded to 4 significant digits, well below the product's
#' radiometric accuracy. The retrieval script is in `data-raw/` in the
#' package sources.
#'
#' @format A 300 x 300 double matrix, linear gamma-nought power (unitless,
#'   all values positive).
#' @source Scene
#'   `S1A_IW_GRDH_1SDV_20240928T020214_20240928T020239_055859_06D3B3`,
#'   `sentinel-1-rtc` collection, Microsoft Planetary Computer. Contains
#'   modified Copernicus Sentinel data (2024).
#' @examples
#' op <- par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
#' rf_plot(s1_sthelens, main = "Sentinel-1 VV")
#' rf_plot(rf_lee(s1_sthelens, window = 7L), main = "Lee filtered")
#' par(op)
"s1_sthelens"
