# We start by defining a few functions. Data preparation starts in the next section. Don't forget to change the working directory.
########## Function definitions ###############

listSDStoImport <- function(searchString, rawDataFilesLocation, dayOffset, dayRange, collectionDates) {
  temperatureFiles <- list.files(path = rawDataFilesLocation, pattern = searchString, full.names = TRUE)
  subFiles <- sapply(paste("A2012", dayOffset + dayRange, sep = ""), grep, x = temperatureFiles, value = TRUE)
  temperatures <- lapply(subFiles, MODIS::getSds)
  splitTemperatures <- split(temperatures, f = factor(substr(subFiles, start = 0, stop = gregexpr(pattern = ".h2", text = subFiles[[1]])[[1]] - 1)))
  names(splitTemperatures) <- collectionDates
  splitTemperatures
}

uniformiseLandCover <- function(landCoverPointsList) {
  landCovers <- do.call("c", lapply(landCoverPointsList, function(x) colnames(x@data)))
  uniqueLandCovers <- unique(landCovers)
  landCoverIndices <- as.numeric(substr(uniqueLandCovers, start = 10, stop = 100))
  uniqueLandCovers <- uniqueLandCovers[order(landCoverIndices)]
  lapply(landCoverPointsList, function(landCoverPoints) {
    if (length(missingCols <- setdiff(uniqueLandCovers, colnames(landCoverPoints@data))) > 0) {
      landCoverPoints@data[missingCols] <- 0
      landCoverPoints@data <- landCoverPoints@data[ , uniqueLandCovers]
    }
    landCoverPoints
  })
}

produceLandCover <- function(landCoverFiles, regionPolygon) { 
  landCoverRasters <- lapply(landCoverFiles, function(filename) {
    landCoverSds <- MODIS::getSds(filename)
    landCover <- raster::raster(readGDAL(landCoverSds$SDS4gdal[2], as.is = TRUE)) # Based on land type classification 2: https://lpdaac.usgs.gov/products/mcd12q1v006/
    landCover
  })
  landCover <- do.call(raster::merge, landCoverRasters)
  smallerRaster <- raster::crop(x = landCover, y = regionPolygon) 
  spObject <- raster::rasterToPoints(smallerRaster, spatial = TRUE)
  indiaValuesIndex <- sp::over(x = spObject, y = regionPolygon)
  pointsInIndia <- subset(spObject, subset = !is.na(indiaValuesIndex))
  raster::values(smallerRaster) <- rep(NA, raster::ncell(smallerRaster))
  output <- raster::rasterize(x = pointsInIndia, y = smallerRaster, field = "layer")
  output
}

prepareDataForMRAinla <- function(temperatures, elevations, landCover, satelliteNamesVec, collectionDatesPOSIX, completeDateVector = collectionDatesPOSIX) {
  if ("RasterLayer" %in% class(temperatures[[1]])) {
    temperaturePoints <- lapply(temperatures, FUN = raster::rasterToPoints, spatial = TRUE)
  } else {
    temperaturePoints <- temperatures
  }

  satelliteNamesList <- lapply(seq_along(satelliteNamesVec), function(dayIndex) {
    rep(satelliteNamesVec[[dayIndex]], nrow(temperaturePoints[[dayIndex]]@coords))
  })
  satellite <- do.call("c", satelliteNamesList)
  satellite <- as.numeric(factor(x = satellite, levels = c("Terra", "Aqua"))) - 1
  numTimePoints <- length(completeDateVector)

  timeValues <- do.call("c", lapply(seq_along(collectionDatesPOSIX), function(x) rep(collectionDatesPOSIX[[x]], length(temperaturePoints[[x]]))))

  timeLevels <- as.numeric(factor(as.character(timeValues), levels = as.character(completeDateVector))) - 1

  timeModelMatrix <- t(sapply(timeLevels, function(x) {
    unitVector <- rep(0, numTimePoints - 1)
    unitVector[x] <- 1 
    unitVector
  }))
  colnames(timeModelMatrix) <- paste("time", 2:numTimePoints, sep = "")

  funToGetLandCoverPoints <- function(tempPoints) {
    tempPointsReproj <- sp::spTransform(tempPoints, crs(landCover))
    landCoverAtPoints <- raster::extract(landCover, tempPointsReproj)
    landCoverValues <- sort(unique(landCoverAtPoints))
    columnNames <- paste("landCover", landCoverValues, sep = "")
    landCoverMatrix <- t(sapply(landCoverAtPoints, function(x) {
      unitVec <- numeric(length(columnNames))
      unitVec[match(x, landCoverValues)] <- 1
      unitVec
    }))
    colnames(landCoverMatrix) <- columnNames
    sp::SpatialPointsDataFrame(coords = tempPoints@coords, data = as.data.frame(landCoverMatrix), proj4string = raster::crs(tempPoints)) 
  }

  landCoverPoints <- lapply(temperaturePoints, FUN = funToGetLandCoverPoints)
  landCoverPoints <- uniformiseLandCover(landCoverPoints)

  elevationPoints <- lapply(temperaturePoints, function(tempPoints) {
    tempPoints <- sp::spTransform(tempPoints, crs(elevations[[1]]))
    elevationValues <- rep(0, length(tempPoints))
    lapply(elevations, function(elevationRaster) {
      extractedValues <- raster::extract(elevationRaster, tempPoints)
      elevationValues[!is.na(extractedValues)] <<- extractedValues[!is.na(extractedValues)]
      NULL
    })
    sp::SpatialPointsDataFrame(coords = tempPoints@coords, data = data.frame(elevation = elevationValues), proj4string = crs(tempPoints))
  })

  latitudePoints <- lapply(temperaturePoints, function(x) {
    sp::SpatialPointsDataFrame(x@coords, data = data.frame(latitude = x@coords[, 2]), proj4string = raster::crs(x))
  })
  combinedData <- do.call("cbind", lapply(list(landCoverPoints, latitudePoints, elevationPoints), function(x) do.call("rbind", lapply(x, function(y) y@data))))
  combinedData <- cbind(combinedData, timeModelMatrix, Aqua = satellite)
  if ("RasterLayer" %in% class(temperatures[[1]])) {
    combinedData <- cbind(do.call("rbind", lapply(temperaturePoints, function(y) y@data)), combinedData)
    colnames(combinedData)[[1]] <- "y"
  }

  coordinates <- do.call("rbind", lapply(temperaturePoints, function(x) x@coords))
  rownames(coordinates) <- as.character(1:nrow(coordinates))
  missingLandCoverOrElevation <- (rowSums(combinedData[ , grep(colnames(combinedData), pattern = "landCover", value = TRUE)]) == 0) | is.na(combinedData[, "elevation"])

  spacetime::STIDF(sp = SpatialPoints(coordinates[!missingLandCoverOrElevation, ], proj4string = crs(temperaturePoints[[1]])), time = timeValues[!missingLandCoverOrElevation], data = as.data.frame(combinedData[!missingLandCoverOrElevation,]))
}

funToCreateRaster <- function(temperatureSdsList, polygonBound) {
  extractionFun <- function(x) {
    tempGrid <- readGDAL(x$SDS4gdal[1], as.is = TRUE)
    hourGrid <- readGDAL(x$SDS4gdal[3], as.is = TRUE)
    tempGrid$band1 <- tempGrid$band1 * 0.02 - 273.15 # See https://gis.stackexchange.com/questions/72524/how-do-i-convert-the-lst-values-on-the-modis-lst-image-to-degree-celsius
    # There's a 0.02 scaling factor applied to values in file to get the temperatures.
    # The -273.15 brings temperatures back in Celsius
    hourGrid@data[,1] <- hourGrid@data[,1] * 0.1

    list(temperatureRaster = raster(tempGrid), hourRaster = raster(hourGrid))
  }
  tempAndTimeRasters <- lapply(temperatureSdsList, extractionFun)

  createRaster <- function(rasterName) {
    rasterList <- lapply(tempAndTimeRasters, function(x) x[[rasterName]])
    mergedRasters <- do.call(raster::merge, rasterList)
    smallerRaster <- crop(x = mergedRasters, y = polygonBound)
    spObject <- rasterToPoints(smallerRaster, spatial = TRUE)
    polygonValuesIndex <- over(x = spObject, y = polygonBound)
    pointsInPolygon <- subset(spObject, subset = !is.na(polygonValuesIndex))
    values(smallerRaster) <- rep(NA, ncell(smallerRaster))
    if (sum(!is.na(polygonValuesIndex)) == 0) {
      return(smallerRaster)
    }
    rasterize(x = pointsInPolygon, y = smallerRaster, field = "layer")
  }
  rasterNames <- c("temperatureRaster", "hourRaster")
  tempAndTime <- lapply(rasterNames, FUN = createRaster)
  names(tempAndTime) <- rasterNames
  tempAndTime
}

funToGetDailyRastersAndSatelliteName <- function(var1, splitTemperaturesBySatellite, MaharashtraPolygonOtherCRS) {
  aquaRasters <- funToCreateRaster(splitTemperaturesBySatellite$Aqua[[var1]], polygonBound = MaharashtraPolygonOtherCRS)
  terraRasters <- funToCreateRaster(splitTemperaturesBySatellite$Terra[[var1]], polygonBound = MaharashtraPolygonOtherCRS)
  if (sum(!is.na(values(aquaRasters$temperatureRaster))) >= sum(!is.na(values(terraRasters$temperatureRaster)))) {
    cat("Returning Aqua!\n")
    c(aquaRasters, satellite = "Aqua")
  } else {
    cat("Returning Terra!\n")
    c(terraRasters, satellite = "Terra")
  }
}

###### END: function definitions #################################

library(sp)
library(raster)
library(MODIS)
library(rgdal)
library(rgeos)
library(spacetime)
library(MRAinla)
library(RhpcBLASctl)
library(geoR)
library(maptools)
library(doParallel)
library(mapmisc)

RandomFields::RFoptions(cores = 1)
blas_set_num_threads(1)
omp_set_num_threads(1)

# The working directory should have subfolder "data".
# Change the following line for your system.

setwd("/home/luc/INLAMRAfiles/INLAMRApaper1/realData")

######## Importing data from MODIS #########################

# This is the folder where the data files downloaded from EarthData are.
# We suggest storing the files in the "data" subfolder.

rawDataFilesLocation <- "data/"

######## Importing MODIS data ####################################

# Naming convention: nnnnnnn.Ayyyyddd.h00v00.vvv.yyyydddhhmmss.
# nnnnnnn: Product name
# Ayyyyddd: Sampling date, year (yyyy), then day (ddd), between 1 and 365.
# h00v00: Identifies the grid tile (see https://lpdaac.usgs.gov/dataset_discovery/modis)
# vvv: Data version
# yyyydddhhmmss: Date data were processed, year, day, hour, minute, second.
dayOffset <- 121
dayRange <- 18:31
collectionDates <- paste("May", dayRange, "_2012", sep = "")
collectionDatesPOSIX <- as.POSIXct(paste("2012-05-", dayRange, sep = ""))

splitTemperaturesBySatellite <- lapply(c(Terra = "MOD11A1.A2012", Aqua = "MYD11A1.A2012"), FUN = listSDStoImport, rawDataFilesLocation = rawDataFilesLocation, dayOffset = dayOffset, dayRange = dayRange, collectionDates = collectionDates)


westMaharashtraPolygonEdges <- rbind(c(19.55, 72.76), c(19.55, 74), c(18.35, 74), c(18.35, 72.76), c(19.55, 72.76))
westMaharashtraPolygonEdges <- westMaharashtraPolygonEdges[ , 2:1]
westMaharashtraPolygon <- SpatialPolygons(Srl = list(Polygons(list(Polygon(coords = westMaharashtraPolygonEdges)), ID = "Mumbai")))
crs(westMaharashtraPolygon) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
westMaharashtraPolygonOtherCRS <- spTransform(westMaharashtraPolygon, CRSobj = CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"))

indiaTemperaturesAndTimes <- lapply(seq_along(splitTemperaturesBySatellite$Aqua), FUN = funToGetDailyRastersAndSatelliteName, splitTemperaturesBySatellite = splitTemperaturesBySatellite, MaharashtraPolygonOtherCRS = westMaharashtraPolygonOtherCRS)

indiaTemperatures <- lapply(indiaTemperaturesAndTimes, function(x) x$temperatureRaster)
satellitePerDay <- sapply(indiaTemperaturesAndTimes, function(x) x$satellite)

landCoverFiles <- list.files(rawDataFilesLocation, pattern = "MCD*", full.names = TRUE)
landCover <- produceLandCover(landCoverFiles, regionPolygon = westMaharashtraPolygonOtherCRS)

# Getting elevation data

elevationFiles <- list.files(path = rawDataFilesLocation, pattern = "*dem.tif", full.names = TRUE)
elevation <- lapply(elevationFiles, raster)

######## END: Importing MODIS data ####################################

######## Generating validation datasets ###############################

may21rasterAddedMissing <- indiaTemperatures[[4]]
may28raster <- indiaTemperatures[[length(indiaTemperatures) - 3]]
indicesForKnownRemovedValues <- which(is.na(values(may28raster)) & !is.na(values(may21rasterAddedMissing)))
values(may21rasterAddedMissing) <- replace(values(may21rasterAddedMissing), indicesForKnownRemovedValues, NA)

mainDataCompleteMap <- prepareDataForMRAinla(landCover = landCover, elevations = elevation, temperatures = c(indiaTemperatures[1:3], may21rasterAddedMissing, indiaTemperatures[5:7]), collectionDatesPOSIX = collectionDatesPOSIX[1:7], satelliteNamesVec = satellitePerDay[1:7])

missingRaster <- indiaTemperatures[[4]]
values(missingRaster) <- NA
values(missingRaster) <- replace(values(missingRaster), indicesForKnownRemovedValues, values(indiaTemperatures[[4]])[indicesForKnownRemovedValues])

testDataMay21 <- prepareDataForMRAinla(landCover = landCover, elevations = elevation, temperatures = list(missingRaster), collectionDatesPOSIX = collectionDatesPOSIX[4], satelliteNamesVec = satellitePerDay[[4]], completeDateVector = collectionDatesPOSIX[1:7])

## landCover = 0 (Water) will be the reference category.

mainDataCompleteMap@data <- subset(mainDataCompleteMap@data, select = -landCover0)
testDataMay21@data <- subset(testDataMay21@data, select = -landCover0)
testDataMay21@data <- subset(testDataMay21@data, select = -y)

# We need to jitter the training data...

set.seed(10)
mainDataCompleteMap@sp@coords <- geoR::jitter2d(mainDataCompleteMap@sp@coords, max = 0.00001)

######## END: Generating validation datasets ###############################

save(mainDataCompleteMap, file = "data/mainDataCompleteMap_May18_24.Rdata", compress = TRUE)
save(testDataMay21, file = "data/testDataMay21_May18_24.Rdata", compress = TRUE)

### END: Setting up datasets
