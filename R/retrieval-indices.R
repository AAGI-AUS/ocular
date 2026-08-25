
## Retrieval metadata: Sentinel-2 spectral indices
s2_index_list <- list(
  NDVI = list(
    assets = c("B08", "B04"),
    fun = function(r) (r$B08 - r$B04) / (r$B08 + r$B04),
    res = 10 ),

  DVI = list(
    assets = c("B08", "B04"),
    fun = function(r) r$B08 - r$B04,
    res = 10 ),

  kNDVI = list(
    assets = c("B08", "B04"),
    fun = function(r) tanh( ((r$B08 - r$B04) / (r$B08 + r$B04))^2 ),
    res = 10 ),

  GNDVI = list(
    assets = c("B08", "B03"),
    fun = function(r) (r$B08 - r$B03) / (r$B08 + r$B03),
    res = 10 ),

  BNDVI = list(
    assets = c("B08", "B02"),
    fun = function(r) (r$B08 - r$B02) / (r$B08 + r$B02),
    res = 10 ),

  GBNDVI = list(
    assets = c("B08", "B03", "B02"),
    fun = function(r) (r$B08 - (r$B03 + r$B02)) / (r$B08 + (r$B03 + r$B02)),
    res = 10 ),

  GRNDVI = list(
    assets = c("B08", "B03", "B04"),
    fun = function(r) (r$B08 - (r$B03 + r$B04)) / (r$B08 + (r$B03 + r$B04)),
    res = 10 ),

  NIRv = list(
    assets = c("B08", "B04"),
    fun = function(r) ((r$B08 - r$B04) / (r$B08 + r$B04)) * r$B08,
    res = 10 ),

  NDWI = list(
    assets = c("B03", "B08"),
    fun = function(r) (r$B03 - r$B08) / (r$B03 + r$B08),
    res = 10 ),

  NDMI = list(
    assets = c("B08", "B11"),
    fun = function(r) (r$B08 - r$B11) / (r$B08 + r$B11),
    res = 20 ),

  EVI = list(
    assets = c("B08", "B04", "B02"),
    fun = function(r) 2.5 * (r$B08 - r$B04) / (r$B08 + 6 * r$B04 - 7.5 * r$B02 + 1),
    res = 10 ),

  SAVI = list(
    assets = c("B08", "B04"),
    fun = function(r) (r$B08 - r$B04) / (r$B08 + r$B04 + 0.5) * 1.5,
    res = 10 ),

  MSAVI = list(
    assets = c("B08", "B04"),
    fun = function(r) (2 * r$B08 + 1 - sqrt((2 * r$B08 + 1)^2 - 8 * (r$B08 - r$B04))) / 2,
    res = 10 ),

  NDRE = list(
    assets = c("B08", "B05"),
    fun = function(r) (r$B08 - r$B05) / (r$B08 + r$B05),
    res = 20 ),

  MTCI = list(
    assets = c("B06", "B05", "B04"),
    fun = function(r) (r$B06 - r$B05) / (r$B05 - r$B04),
    res = 20 ),

  CIre_gitelson = list(
    assets = c("B07", "B05"),
    fun = function(r) (r$B07 / r$B05) - 1,
    res = 20 ),

  CIRE = list(
    assets = c("B08", "B05"),
    fun = function(r) (r$B08 / r$B05) - 1,
    res = 20 ),

  PSRI = list(
    assets = c("B06", "B04", "B02"),
    fun = function(r) (r$B04 - r$B02) / r$B06,
    res = 20 ),

  BSI = list(
    assets = c("B11", "B08", "B04", "B02"),
    fun = function(r) ((r$B11 + r$B04) - (r$B08 + r$B02)) / ((r$B11 + r$B04) + (r$B08 + r$B02)),
    res = 20 ),

  SR = list(
    assets = c("B08", "B04"),
    fun = function(r) r$B08 / r$B04,
    res = 10 ),

  EVI2 = list(
    assets = c("B08", "B04"),
    fun = function(r) 2.5 * (r$B08 - r$B04) / (r$B08 + 2.4 * r$B04 + 1),
    res = 10 ),

  OSAVI = list(
    assets = c("B08", "B04"),
    fun = function(r) (r$B08 - r$B04) / (r$B08 + r$B04 + 0.16),
    res = 10 ),

  ARVI = list(
    assets = c("B08", "B04", "B02"),
    fun = function(r) (r$B08 - (2 * r$B04 - r$B02)) / (r$B08 + (2 * r$B04 - r$B02)),
    res = 10 ),

  GARI = list(
    assets = c("B08", "B03", "B02", "B04"),
    fun = function(r) (r$B08 - (r$B03 - (r$B02 - r$B04))) / (r$B08 + (r$B03 - (r$B02 - r$B04))),
    res = 10 ),

  CIG = list(
    assets = c("B08", "B03"),
    fun = function(r) (r$B08 / r$B03) - 1,
    res = 10 ),

  MNDWI = list(
    assets = c("B03", "B11"),
    fun = function(r) (r$B03 - r$B11) / (r$B03 + r$B11),
    res = 20 ),

  NBR = list(
    assets = c("B08", "B12"),
    fun = function(r) (r$B08 - r$B12) / (r$B08 + r$B12),
    res = 20 ),

  NBR2 = list(
    assets = c("B11", "B12"),
    fun = function(r) (r$B11 - r$B12) / (r$B11 + r$B12),
    res = 20 ),

  BAI = list(
    assets = c("B04", "B08"),
    fun = function(r) 1 / ((0.1 - r$B04)^2 + (0.06 - r$B08)^2),
    res = 10 ),

  MCARI = list(
    assets = c("B05", "B04", "B03"),
    fun = function(r) ((r$B05 - r$B04) - 0.2 * (r$B05 - r$B03)) * (r$B05 / r$B04),
    res = 20 ),

  IRECI = list(
    assets = c("B07", "B04", "B05", "B06"),
    fun = function(r) (r$B07 - r$B04) / (r$B05 / r$B06),
    res = 20 )
)

## Landsat spectral indices list
landsat_index_list <- list(
  NDVI = list(
    assets = c("nir08", "red"),
    fun = function(r) (r$nir08 - r$red) / (r$nir08 + r$red),
    res = 30 ),

  DVI = list(
    assets = c("nir08", "red"),
    fun = function(r) r$nir08 - r$red,
    res = 30 ),

  kNDVI = list(
    assets = c("nir08", "red"),
    fun = function(r) tanh( ((r$nir08 - r$red) / (r$nir08 + r$red))^2 ),
    res = 30 ),

  GNDVI = list(
    assets = c("nir08", "green"),
    fun = function(r) (r$nir08 - r$green) / (r$nir08 + r$green),
    res = 30 ),

  BNDVI = list(
    assets = c("nir08", "blue"),
    fun = function(r) (r$nir08 - r$blue) / (r$nir08 + r$blue),
    res = 30 ),

  GBNDVI = list(
    assets = c("nir08", "green", "blue"),
    fun = function(r) (r$nir08 - (r$green + r$blue)) / (r$nir08 + (r$green + r$blue)),
    res = 30 ),

  GRNDVI = list(
    assets = c("nir08", "green", "red"),
    fun = function(r) (r$nir08 - (r$green + r$red)) / (r$nir08 + (r$green + r$red)),
    res = 30 ),

  NIRv = list(
    assets = c("nir08", "red"),
    fun = function(r) ((r$nir08 - r$red) / (r$nir08 + r$red)) * r$nir08,
    res = 30 ),

  NDWI = list(
    assets = c("green", "nir08"),
    fun = function(r) (r$green - r$nir08) / (r$green + r$nir08),
    res = 30 ),

  NDMI = list(
    assets = c("nir08", "swir16"),
    fun = function(r) (r$nir08 - r$swir16) / (r$nir08 + r$swir16),
    res = 30 ),

  EVI = list(
    assets = c("nir08", "red", "blue"),
    fun = function(r) 2.5 * (r$nir08 - r$red) / (r$nir08 + 6 * r$red - 7.5 * r$blue + 1),
    res = 30 ),

  SAVI = list(
    assets = c("nir08", "red"),
    fun = function(r) (r$nir08 - r$red) / (r$nir08 + r$red + 0.5) * 1.5,
    res = 30 ),

  MSAVI = list(
    assets = c("nir08", "red"),
    fun = function(r) (2 * r$nir08 + 1 - sqrt((2 * r$nir08 + 1)^2 - 8 * (r$nir08 - r$red))) / 2,
    res = 30 ),

  BSI = list(
    assets = c("swir16", "nir08", "red", "blue"),
    fun = function(r) ((r$swir16 + r$red) - (r$nir08 + r$blue)) / ((r$swir16 + r$red) + (r$nir08 + r$blue)),
    res = 30 ),

  SR = list(
    assets = c("nir08", "red"),
    fun = function(r) r$nir08 / r$red,
    res = 30 ),

  EVI2 = list(
    assets = c("nir08", "red"),
    fun = function(r) 2.5 * (r$nir08 - r$red) / (r$nir08 + 2.4 * r$red + 1),
    res = 30 ),

  OSAVI = list(
    assets = c("nir08", "red"),
    fun = function(r) (r$nir08 - r$red) / (r$nir08 + r$red + 0.16),
    res = 30 ),

  ARVI = list(
    assets = c("nir08", "red", "blue"),
    fun = function(r) (r$nir08 - (2 * r$red - r$blue)) / (r$nir08 + (2 * r$red - r$blue)),
    res = 30 ),

  GARI = list(
    assets = c("nir08", "green", "blue", "red"),
    fun = function(r) (r$nir08 - (r$green - (r$blue - r$red))) / (r$nir08 + (r$green - (r$blue - r$red))),
    res = 30 ),

  CIG = list(
    assets = c("nir08", "green"),
    fun = function(r) (r$nir08 / r$green) - 1,
    res = 30 ),

  MNDWI = list(
    assets = c("green", "swir16"),
    fun = function(r) (r$green - r$swir16) / (r$green + r$swir16),
    res = 30 ),

  NBR = list(
    assets = c("nir08", "swir22"),
    fun = function(r) (r$nir08 - r$swir22) / (r$nir08 + r$swir22),
    res = 30 ),

  NBR2 = list(
    assets = c("swir16", "swir22"),
    fun = function(r) (r$swir16 - r$swir22) / (r$swir16 + r$swir22),
    res = 30 ),

  BAI = list(
    assets = c("red", "nir08"),
    fun = function(r) 1 / ((0.1 - r$red)^2 + (0.06 - r$nir08)^2),
    res = 30 )

)

## ---------------------------------------------------------------------------
## MCD43A4 (modis-43A4-061) -- daily 500 m NBAR (BRDF-normalised reflectance).
## Every index is computed from NBAR reflectance bands. Band map: Band1 red,
## Band2 NIR, Band3 blue, Band7 SWIR2 (used as MIR for NBR). Band order per
## index mirrors landsat_index_list so the fusion functions can pair fine and
## coarse bands positionally. res = 500.
## Confirm the asset names ("Nadir_Reflectance_BandN") against the live
## modis-43A4-061 collection before production use.
## ---------------------------------------------------------------------------
mcd43a4_index_list <- list(
  NDVI = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band1"),
    fun    = function(r) (r[[1L]] - r[[2L]]) / (r[[1L]] + r[[2L]]),
    res    = 500 ),

  EVI = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band1", "Nadir_Reflectance_Band3"),
    fun    = function(r) 2.5 * (r[[1L]] - r[[2L]]) / (r[[1L]] + 6 * r[[2L]] - 7.5 * r[[3L]] + 1),
    res    = 500 ),

  DVI = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band1"),
    fun    = function(r) r[[1L]] - r[[2L]],
    res    = 500 ),

  kNDVI = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band1"),
    fun    = function(r) tanh( ((r[[1L]] - r[[2L]]) / (r[[1L]] + r[[2L]]))^2 ),
    res    = 500 ),

  BNDVI = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band3"),
    fun    = function(r) (r[[1L]] - r[[2L]]) / (r[[1L]] + r[[2L]]),
    res    = 500 ),

  NIRv = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band1"),
    fun    = function(r) ((r[[1L]] - r[[2L]]) / (r[[1L]] + r[[2L]])) * r[[1L]],
    res    = 500 ),

  SAVI = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band1"),
    fun    = function(r) (r[[1L]] - r[[2L]]) / (r[[1L]] + r[[2L]] + 0.5) * 1.5,
    res    = 500 ),

  MSAVI = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band1"),
    fun    = function(r) (2 * r[[1L]] + 1 - sqrt((2 * r[[1L]] + 1)^2 - 8 * (r[[1L]] - r[[2L]]))) / 2,
    res    = 500 ),

  SR = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band1"),
    fun    = function(r) r[[1L]] / r[[2L]],
    res    = 500 ),

  EVI2 = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band1"),
    fun    = function(r) 2.5 * (r[[1L]] - r[[2L]]) / (r[[1L]] + 2.4 * r[[2L]] + 1),
    res    = 500 ),

  OSAVI = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band1"),
    fun    = function(r) (r[[1L]] - r[[2L]]) / (r[[1L]] + r[[2L]] + 0.16),
    res    = 500 ),

  ARVI = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band1", "Nadir_Reflectance_Band3"),
    fun    = function(r) (r[[1L]] - (2 * r[[2L]] - r[[3L]])) / (r[[1L]] + (2 * r[[2L]] - r[[3L]])),
    res    = 500 ),

  NBR = list(
    assets = c("Nadir_Reflectance_Band2", "Nadir_Reflectance_Band7"),
    fun    = function(r) (r[[1L]] - r[[2L]]) / (r[[1L]] + r[[2L]]),
    res    = 500 ),

  BAI = list(
    assets = c("Nadir_Reflectance_Band1", "Nadir_Reflectance_Band2"),
    fun    = function(r) 1 / ((0.1 - r[[1L]])^2 + (0.06 - r[[2L]])^2),
    res    = 500 )
)
