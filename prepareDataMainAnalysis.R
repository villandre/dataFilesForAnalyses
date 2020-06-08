# We start by defining a few functions. The script starts afterwards: don't forget to change the working directory.
############ Function definitions ###############

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
    unitVector[x] <- 1 # When x takes value 0, the vector remains all 0s, which is what we want.
    unitVector
  }))
  colnames(timeModelMatrix) <- paste("time", 2:numTimePoints, sep = "")

  landCoverPoints <- lapply(temperaturePoints, function(tempPoints) {
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
  })
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
    tempGrid <- rgdal::readGDAL(x$SDS4gdal[1], as.is = TRUE)
    hourGrid <- rgdal::readGDAL(x$SDS4gdal[3], as.is = TRUE)
    tempGrid$band1 <- tempGrid$band1 * 0.02 - 273.15 # See https://gis.stackexchange.com/questions/72524/how-do-i-convert-the-lst-values-on-the-modis-lst-image-to-degree-celsius
    # There's a 0.02 scaling factor applied to values in file to get the temperatures.
    # The -273.15 brings temperatures back in Celsius
    hourGrid@data[,1] <- hourGrid@data[,1] * 0.1

    list(temperatureRaster = raster::raster(tempGrid), hourRaster = raster::raster(hourGrid))
  }
  tempAndTimeRasters <- lapply(temperatureSdsList, extractionFun)

  createRaster <- function(rasterName) {
    rasterList <- lapply(tempAndTimeRasters, function(x) x[[rasterName]])
    mergedRasters <- do.call(raster::merge, rasterList)
    smallerRaster <- raster::crop(x = mergedRasters, y = polygonBound)
    spObject <- raster::rasterToPoints(smallerRaster, spatial = TRUE)
    polygonValuesIndex <- sp::over(x = spObject, y = polygonBound)
    pointsInPolygon <- subset(spObject, subset = !is.na(polygonValuesIndex))
    values(smallerRaster) <- rep(NA, ncell(smallerRaster))
    if (sum(!is.na(polygonValuesIndex)) == 0) {
      return(smallerRaster)
    }
    raster::rasterize(x = pointsInPolygon, y = smallerRaster, field = "layer")
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

produceTestData <- function(indiaTemperatures, landCover, elevation, collectionDatesPOSIX, boundaryPolygon, satelliteNamesVec, dayIndex) {
  day28raster <- indiaTemperatures[[dayIndex]]
  raster::values(day28raster) <- replace(values(day28raster), is.na(values(day28raster)), -50)
  rasterCellsMidpoints <- raster::rasterToPoints(day28raster, spatial = TRUE)
  indiaValuesIndex <- sp::over(x = rasterCellsMidpoints, y = boundaryPolygon)
  pointsInIndia <- subset(rasterCellsMidpoints, subset = !is.na(indiaValuesIndex))
  missingIndianValueIndices <- pointsInIndia@data$layer == -50
  missingPoints <- subset(pointsInIndia, subset = missingIndianValueIndices)

  landCoverInMissingZones <- raster::extract(x = landCover, y = missingPoints)
  missingPointsNoWaterNoNA <- missingPoints[which(!is.na(landCoverInMissingZones) & !(landCoverInMissingZones == 0)), ]

  emptyRaster <- day28raster
  raster::values(emptyRaster) <- rep(NA, ncell(emptyRaster))
  missingRaster <- raster::rasterize(x = missingPointsNoWaterNoNA, y = emptyRaster, field = "layer")

  testDataMay28 <- prepareDataForMRAinla(landCover = landCover, elevations = elevation, temperatures = list(missingRaster), collectionDatesPOSIX = collectionDatesPOSIX[length(collectionDatesPOSIX) - 3], satelliteNamesVec = satelliteNamesVec, completeDateVector = collectionDatesPOSIX)
  testDataMay28
}

########## END: Function definitions ##################
#######################################################

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

# MODIS naming convention: nnnnnnn.Ayyyyddd.h00v00.vvv.yyyydddhhmmss.
# nnnnnnn: Product name
# Ayyyyddd: Sampling date, year (yyyy), then day (ddd), between 1 and 365.
# h00v00: Identifies the grid tile (see https://lpdaac.usgs.gov/dataset_discovery/modis)
# vvv: Data version
# yyyydddhhmmss: Date data were processed, year, day, hour, minute, second.
dayOffset <- 121 # Days are listed numerically, from 1 to 366 (leap year). Day 121 corresponds to April 30.
dayRange <- 25:31 # We extract data for May 25 until May 31.
collectionDates <- paste("May", dayRange, "_2012", sep = "")
collectionDatesPOSIX <- as.POSIXct(paste("2012-05-", dayRange, sep = ""))

splitTemperaturesBySatellite <- lapply(c(Terra = "MOD11A1.A2012", Aqua = "MYD11A1.A2012"), FUN = listSDStoImport, rawDataFilesLocation = rawDataFilesLocation, dayOffset = dayOffset, dayRange = dayRange, collectionDates = collectionDates)

indiaPolygons <- raster::getData(country = "IND", level = 2)

# These are the latitude-longitude boundaries that cicumscribe the region we picked for the analysis, that corresponds roughly to western Maharashtra.
MaharashtraPolygonEdges <- rbind(c(21, 72.76), c(21, 76.5), c(17.0, 76.5), c(17.0, 72.6), c(21, 72.6))
MaharashtraPolygonEdges <- MaharashtraPolygonEdges[ , 2:1]
MaharashtraPolygon <- SpatialPolygons(Srl = list(Polygons(list(Polygon(coords = MaharashtraPolygonEdges)), ID = "Mumbai")))
crs(MaharashtraPolygon) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
MaharashtraPolygonOtherCRS <- spTransform(MaharashtraPolygon, CRSobj = CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"))

# On each day, we extract either the Aqua or Terra satellite data, depending on which has fewer missing values.

indiaTemperaturesAndTimes <- lapply(seq_along(splitTemperaturesBySatellite$Aqua), FUN = funToGetDailyRastersAndSatelliteName, splitTemperaturesBySatellite = splitTemperaturesBySatellite, MaharashtraPolygonOtherCRS = MaharashtraPolygonOtherCRS)

indiaTemperatures <- lapply(indiaTemperaturesAndTimes, function(x) x$temperatureRaster)
names(indiaTemperatures) <- as.character(collectionDatesPOSIX)
satellitePerDay <- sapply(indiaTemperaturesAndTimes, function(x) x$satellite)

# We convert the first raster to the standard coordinates and plot the result...

landCoverFiles <- list.files(rawDataFilesLocation, pattern = "MCD*", full.names = TRUE)

landCover <- produceLandCover(landCoverFiles, regionPolygon = MaharashtraPolygonOtherCRS)

# Getting elevation data

elevationFiles <- list.files(path = rawDataFilesLocation, pattern = "*dem.tif", full.names = TRUE)
elevation <- lapply(elevationFiles, raster)

######## END: Importing data from MODIS #########################

######## Creating datasets for main analysis ####################
mainDataCompleteMap <- prepareDataForMRAinla(landCover = landCover, elevations = elevation, temperatures = indiaTemperatures, collectionDatesPOSIX = collectionDatesPOSIX, satelliteNamesVec = satellitePerDay)

# The test data in the main analysis consists of all missing LSTs on May 28th (excluding oceanic tiles). The goal is to complete the map.

testDataMay28 <- produceTestData(indiaTemperatures = indiaTemperatures, boundaryPolygon = MaharashtraPolygonOtherCRS, dayIndex = length(indiaTemperatures) - 3, satelliteNamesVec = satellitePerDay[[4]], landCover = landCover, elevation = elevation, collectionDatesPOSIX = collectionDatesPOSIX)

# landCover = 0 (water) will be the reference category.

mainDataCompleteMap@data <- subset(mainDataCompleteMap@data, select = -landCover0)
if ("landCover0" %in% colnames(testDataMay28@data)) {
  testDataMay28@data <- subset(testDataMay28@data, select = -landCover0)
}

testDataMay28@data <- subset(testDataMay28@data, select = -y)

# We very lightly jitter the training data...

set.seed(10)
mainDataCompleteMap@sp@coords <- geoR::jitter2d(mainDataCompleteMap@sp@coords, max = 0.00001)

######## END: Creating datasets for main analysis ####################

save(mainDataCompleteMap, file = "data/mainDataCompleteMap_May25_31_Larger.Rdata", compress = TRUE)
save(testDataMay28, file = "data/testDataMay28_May25_31_Larger.Rdata", compress = TRUE)
