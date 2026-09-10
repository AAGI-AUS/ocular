test_that("named bands retain their identities when indices are computed", {
  r <- terra::rast(nrows=2, ncols=2, xmin=0, xmax=2, ymin=0, ymax=2)
  red <- terra::setValues(r, c(.1, .2, .3, .4))
  nir <- terra::setValues(r, c(.5, .6, .7, .8))
  names(red) <- "downloaded_red_asset"
  names(nir) <- "downloaded_nir_asset"
  bands <- list(B04=red, B08=nir)
  got <- ocular:::.computeIndexFromBands(bands, ocular:::s2_index_list$NDVI)
  expect_s4_class(got, "SpatRaster")
  expect_equal(terra::values(got, mat=FALSE), c(2/3, .5, .4, 1/3),
               tolerance=1e-12)
  mask <- terra::setValues(r, c(1, 0, 0, 0))
  masked <- lapply(bands, function(b) terra::mask(b, mask, maskvalue=1))
  got <- ocular:::.computeIndexFromBands(masked, ocular:::s2_index_list$NDVI)
  expect_true(is.na(terra::values(got, mat=FALSE)[1]))
  expect_equal(terra::values(got, mat=FALSE)[-1], c(.5, .4, 1/3),
               tolerance=1e-12)
  ls <- ocular:::.computeIndexFromBands(list(red=red, nir08=nir),
                                      ocular:::landsat_index_list$NDVI)
  expect_equal(terra::values(ls,mat=FALSE), c(2/3, .5, .4, 1/3),
               tolerance=1e-12)
})

test_that("date names do not become raster merge arguments", {
  left <- terra::rast(nrows=2, ncols=2, xmin=0, xmax=2, ymin=0, ymax=2,
                      crs="EPSG:4326", vals=1:4)
  right <- terra::rast(nrows=2, ncols=2, xmin=2, xmax=4, ymin=0, ymax=2,
                       crs="EPSG:4326", vals=5:8)
  rasters <- setNames(list(left, right), c("2024-06-01", "2024-06-17"))
  got <- ocular:::.mergeRasterTiles(rasters)
  expect_s4_class(got, "SpatRaster")
  expect_equal(as.vector(terra::ext(got)), c(xmin=0,xmax=4,ymin=0,ymax=2))
  expect_equal(terra::values(got, mat=FALSE), c(1,2,5,6,3,4,7,8))
})

test_that("STAC retrieval stacks named assets before masking and mosaicking", {
  for (source in c("sentinel-2", "landsat-8", "mcd43a4")) {
    local({
      sensor <- source
      cfg <- ocular:::.stacCfg(sensor)
      lookup <- switch(sensor,
        "sentinel-2"=ocular:::s2_index_list$NDVI,
        "landsat-8"=ocular:::landsat_index_list$NDVI,
        "mcd43a4"=ocular:::mcd43a4_index_list$NDVI)
      r <- terra::rast(nrows=2,ncols=2,xmin=0,xmax=2,ymin=0,ymax=2,
                       crs="EPSG:4326")
      raw <- c(10000, 20000)
      if (sensor != "landsat-8") raw <- c(1000, 5000)
      red_name <- if (sensor == "sentinel-2") "B04" else
        if (sensor == "landsat-8") "red" else "Nadir_Reflectance_Band1"
      values <- setNames(lapply(lookup$assets, function(an)
        terra::setValues(r, if(an == red_name) raw[1] else raw[2])),
        lookup$assets)
      if (!is.null(cfg$mask_asset))
        values[[cfg$mask_asset]] <- terra::setValues(r,
          c(if(sensor == "sentinel-2") 9 else 8, 0, 0, 0))
      assets <- setNames(lapply(names(values), function(an)
        list(href=paste0("fixture-",an))), names(values))
      items <- list(features=lapply(c("2024-06-01", "2024-06-17"),
        function(date) list(properties=list(datetime=paste0(date,"T00:00:00Z"),
          platform="landsat-8"), assets=assets)))
      original_rast <- terra::rast
      testthat::local_mocked_bindings(rast=function(x, ...) {
        if (missing(x)) return(original_rast(...))
        if (is.character(x) && length(x)==1L &&
            startsWith(x,"/vsicurl/fixture-")) {
          out <- values[[sub("/vsicurl/fixture-","",x,fixed=TRUE)]]
          names(out) <- "downloaded_asset"
          return(out)
        }
        original_rast(x,...)
      }, .package="terra")
      testthat::local_mocked_bindings(.fetchStacItems=function(...) items,
                                      .package="ocular")
      bbox <- sf::st_bbox(c(xmin=0,ymin=0,xmax=2,ymax=2), crs=4326)
      mask_classes <- if(sensor == "sentinel-2") 9L else
        if(sensor == "landsat-8") 3L else NULL
      got <- ocular:::.fetchStac(bbox,"2024-06-01","2024-06-18",
        "NDVI",50,scl_classes=mask_classes,source=sensor)
      expect_length(got,2L)
      scaled <- raw*cfg$scale_factor + cfg$scale_offset
      expected <- rep(diff(scaled)/sum(scaled),4)
      if (!is.null(mask_classes)) expected[1] <- NA_real_
      for (scene in got) {
        expect_s4_class(scene$index,"SpatRaster")
        expect_equal(terra::values(scene$index,mat=FALSE),expected,
                     tolerance=1e-12)
      }
    })
  }
})
