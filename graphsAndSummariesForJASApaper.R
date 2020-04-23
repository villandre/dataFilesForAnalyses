########## Function definitions ###############

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

# completeDateVector is for building the time design matrix when data in temperatures were not sampled across the entirety
# of the period of interest, i.e. training data were sampled on days 1-7, but test data are only for day 4.

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
    sp::SpatialPointsDataFrame(coords = tempPoints@coords, data = as.data.frame(landCoverMatrix), proj4string = raster::crs(tempPoints)) # A line in data with only zeros corresponds to a missing value.
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

produceLandCover <- function(landCoverFiles) {
  landCoverRasters <- lapply(landCoverFiles, function(filename) {
    landCoverSds <- getSds(filename)
    landCover <- raster(readGDAL(landCoverSds$SDS4gdal[2], as.is = TRUE)) # Based on land type classification 2: https://lpdaac.usgs.gov/products/mcd12q1v006/
    landCover
  })
  landCover <- do.call(raster::merge, landCoverRasters)
  landCover
}

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
library(xtable)

RandomFields::RFoptions(cores = 1)
blas_set_num_threads(1)
omp_set_num_threads(1)

# Input the correct working directory for your system
# **The working directory should have two subfolders: "data" for datasets and "outputFiles" for model fitting outputs and graphs.**

setwd("/home/luc/INLAMRAfiles/INLAMRApaper1/realData")

# Input in variable rawDataFilesLocation the folder name where MODIS data files are found

# rawDataFilesLocation <- "data/"
rawDataFilesLocation <- "/store/luc/rawDataFiles"

######### Loading datasets necessary for analyses ########################

trainingDataValidationName <- load("data/mainDataCompleteMap_May18_24.Rdata")
trainingDataValidation <- get(trainingDataValidationName)
rm(list = trainingDataValidationName)

testDataValidationName <- load("data/testDataMay21_May18_24.Rdata")
testDataValidation <- get(testDataValidationName)
rm(list = testDataValidationName)

indiaAnalysisValidationName <- load("outputFiles/INLAMRA_validationAnalysis.Rdata")
indiaAnalysisValidation <- get(indiaAnalysisValidationName)

if (indiaAnalysisValidationName != "indiaAnalysisValidation") {
  rm(list = indiaAnalysisValidationName)
}

SPDEresultName <- load("outputFiles/inlaFitForValidationAnalysis.Rdata")
SPDEresult <- get(SPDEresultName)
rm(list = SPDEresultName)

predictionDataMainName <- load("data/testDataMay28_May25_31_Larger.Rdata")
predictionDataMain <- get(predictionDataMainName)
rm(list = predictionDataMainName)

trainingDataMainName <- load("data/mainDataCompleteMap_May25_31_Larger.Rdata")
trainingDataMain <- get(trainingDataMainName)
rm(list = trainingDataMainName)

indiaAnalysisMainName <- load("outputFiles/indiaAnalysisMay28evenLarger_May25_31_VersionApril13.Rdata")
indiaAnalysisMain <- get(indiaAnalysisMainName)
rm(list = indiaAnalysisMainName)

######### END: Loading datasets necessary for analyses ########################

######## Obtaining basic LST plots ########

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

splitTemperaturesBySatellite <- lapply(c(Terra = "MOD11A1.A2012", Aqua = "MYD11A1.A2012"), function(searchString) {
  temperatureFiles <- list.files(path = rawDataFilesLocation, pattern = searchString, full.names = TRUE)
  # temperatureFiles <- list.files(path = "data", pattern = searchString, full.names = TRUE)
  subFiles <- sapply(paste("A2012", dayOffset + dayRange, sep = ""), grep, x = temperatureFiles, value = TRUE)
  temperatures <- lapply(subFiles, getSds)
  splitTemperatures <- split(temperatures, f = factor(substr(subFiles, start = 0, stop = gregexpr(pattern = ".h2", text = subFiles[[1]])[[1]] - 1)))
  names(splitTemperatures) <- collectionDates
  splitTemperatures
})

indiaPolygons <- raster::getData(country = "IND", level = 2)
indiaPolygonsOtherCRS <- spTransform(indiaPolygons, CRSobj = crs(readGDAL(splitTemperaturesBySatellite$Aqua[[1]][[1]]$SDS4gdal[1], as.is = TRUE)))

largeWesternMahaPolygonEdges <- rbind(c(21, 72.76), c(21, 76.5), c(17.0, 76.5), c(17.0, 72.6), c(21, 72.6))
largeWesternMahaPolygonEdges <- largeWesternMahaPolygonEdges[ , 2:1]
largeWesternMahaPolygon <- SpatialPolygons(Srl = list(Polygons(list(Polygon(coords = largeWesternMahaPolygonEdges)), ID = "Mumbai")))
crs(largeWesternMahaPolygon) <- crs(indiaPolygons)
largeWesternMahaPolygonOtherCRS <- spTransform(largeWesternMahaPolygon, CRSobj = crs(readGDAL(splitTemperaturesBySatellite$Aqua[[1]][[1]]$SDS4gdal[1], as.is = TRUE)))

lastWeekIndices <- collectionDatesPOSIX >= as.POSIXct("2012-05-25")
indiaTemperaturesAndTimesLastWeek <- lapply(seq_along(splitTemperaturesBySatellite$Aqua)[lastWeekIndices], function(var1) {
  aquaRasters <- funToCreateRaster(splitTemperaturesBySatellite$Aqua[[var1]], polygonBound = largeWesternMahaPolygonOtherCRS)
  terraRasters <- funToCreateRaster(splitTemperaturesBySatellite$Terra[[var1]], polygonBound = largeWesternMahaPolygonOtherCRS)
  if (sum(!is.na(values(aquaRasters$temperatureRaster))) >= sum(!is.na(values(terraRasters$temperatureRaster)))) {
    cat("Returning Aqua!\n")
    c(aquaRasters, satellite = "Aqua")
  } else {
    cat("Returning Terra!\n")
    c(terraRasters, satellite = "Terra")
  }
})

indiaTemperaturesLastWeek <- lapply(indiaTemperaturesAndTimesLastWeek, function(x) x$temperatureRaster)
names(indiaTemperaturesLastWeek) <- as.character(collectionDatesPOSIX[lastWeekIndices])
indiaTimesLastWeek <- lapply(indiaTemperaturesAndTimesLastWeek, function(x) x$hourRaster)
satellitePerDayLastWeek <- sapply(indiaTemperaturesAndTimesLastWeek, function(x) x$satellite)

citiesNames <- c("Mumbai City", "Pune", "Nashik", "Ahmadnagar", "Aurangabad")
cityPoints <- mapmisc::geocode(x = citiesNames)
cityPointsOtherCRS <- spTransform(cityPoints, CRSobj = crs(readGDAL(splitTemperaturesBySatellite$Aqua[[1]][[1]]$SDS4gdal[1], as.is = TRUE)))
indiaExtent <- extent(c(xmin = min(largeWesternMahaPolygonEdges[ , 1]), xmax = max(largeWesternMahaPolygonEdges[ , 1]), ymin = min(largeWesternMahaPolygonEdges[ , 2]), ymax = max(largeWesternMahaPolygonEdges[ , 2])))
indiaRasterReprojected <- raster(x = indiaExtent, nrows = nrow(indiaTemperaturesLastWeek[[1]]), ncols = ncol(indiaTemperaturesLastWeek[[1]]), crs = crs(cityPoints))

lapply(c("2012-05-27", "2012-05-28", "2012-05-29"), function(dateName) {

  rasterReprojected <- projectRaster(from = indiaTemperaturesLastWeek[[dateName]],
                                     to = indiaRasterReprojected)
  ecol <- mapmisc::colourScale(values(rasterReprojected), col = "Spectral", breaks = 10, rev = TRUE, style = "equal", dec = 1, opacity = 0.5)
  filename <- paste("outputFiles/temperatures", dateName, ".jpg", sep = "")
  jpeg(filename = filename, width = 1200, height = 1200)
  raster::plot(rasterReprojected, legend = FALSE, cex.axis = 2.5, col = ecol$col, breaks = ecol$breaks)
  plot(indiaPolygons, add = TRUE)
  # raster::plot(indiaTemperaturesLastWeek[[dateName]], legend.only = TRUE, legend.width = 4, axis.args = list(cex.axis = 3))
  mapmisc::legendBreaks("bottomleft", ecol, title = "LST (Celsius)", bg = "white", cex = 2.5)
  plot(cityPoints, pch = 15, col = "red", cex = 5, add = TRUE)
  text(x = cityPoints@coords[ , 1], y = cityPoints@coords[ , 2], labels = replace(citiesNames, 1, "Mumbai"), offset = 2, cex = 4, pos = 4)
  scaleBar(crs = crs(indiaPolygons), pos = "topleft", cex = 3, pt.cex = 2.2, title.cex = 3.5)
  dev.off()
  NULL
})

######## END: Obtaining basic LST plots ########

######## Obtaining plots for presenting validation task ########

smallWesternMahaPolygonEdges <- rbind(c(19.55, 72.76), c(19.55, 74), c(18.35, 74), c(18.35, 72.76), c(19.55, 72.76))
smallWesternMahaPolygonEdges <- smallWesternMahaPolygonEdges[ , 2:1]
smallWesternMahaPolygon <- SpatialPolygons(Srl = list(Polygons(list(Polygon(coords = smallWesternMahaPolygonEdges)), ID = "Mumbai")))
crs(smallWesternMahaPolygon) <- crs(indiaPolygons)
smallWesternMahaPolygonOtherCRS <- spTransform(smallWesternMahaPolygon, CRSobj = crs(readGDAL(splitTemperaturesBySatellite$Aqua[[1]][[1]]$SDS4gdal[1], as.is = TRUE)))

secondToLastWeekIndices <- (collectionDatesPOSIX >= as.POSIXct("2012-05-18")) & (collectionDatesPOSIX < as.POSIXct("2012-05-25"))
indiaTemperaturesAndTimesSecondToLastWeek <- lapply(seq_along(splitTemperaturesBySatellite$Aqua)[secondToLastWeekIndices], function(var1) {
 aquaRasters <- funToCreateRaster(splitTemperaturesBySatellite$Aqua[[var1]], polygonBound = smallWesternMahaPolygonOtherCRS)
 terraRasters <- funToCreateRaster(splitTemperaturesBySatellite$Terra[[var1]], polygonBound = smallWesternMahaPolygonOtherCRS)
 if (sum(!is.na(values(aquaRasters$temperatureRaster))) >= sum(!is.na(values(terraRasters$temperatureRaster)))) {
   cat("Returning Aqua!\n")
   c(aquaRasters, satellite = "Aqua")
 } else {
   cat("Returning Terra!\n")
   c(terraRasters, satellite = "Terra")
 }
})
names(indiaTemperaturesAndTimesSecondToLastWeek) <- collectionDatesPOSIX[secondToLastWeekIndices]

indiaTemperaturesSecondToLastWeek <- lapply(indiaTemperaturesAndTimesSecondToLastWeek, function(x) x$temperatureRaster)
names(indiaTemperaturesSecondToLastWeek) <- names(indiaTemperaturesAndTimesSecondToLastWeek)
indiaTimesSecondToLastWeek <- lapply(indiaTemperaturesAndTimesSecondToLastWeek, function(x) x$hourRaster)
satellitePerDaySecondToLastWeek <- sapply(indiaTemperaturesAndTimesSecondToLastWeek, function(x) x$satellite)
smallRasterExtent <- extent(c(xmin = min(smallWesternMahaPolygonEdges[ , 1]), xmax = max(smallWesternMahaPolygonEdges[ , 1]), ymin = min(smallWesternMahaPolygonEdges[ , 2]), ymax = max(smallWesternMahaPolygonEdges[ , 2])))
emptyRaster <- raster(x = smallRasterExtent, nrows = 120, ncols = 120, crs = crs(indiaPolygons))

######## Producing plot to represent validation test set ########

may21raster <- may21rasterAddedMissing <- indiaTemperaturesSecondToLastWeek[["2012-05-21"]]
may28raster <- crop(x = indiaTemperaturesLastWeek[["2012-05-28"]], y = extent(may21raster))
indicesForKnownRemovedValues <- which(is.na(values(may28raster)) & !is.na(values(may21rasterAddedMissing)))
values(may21rasterAddedMissing) <- replace(values(may21rasterAddedMissing), indicesForKnownRemovedValues, NA)
may21rasterAddedMissingReproj <- projectRaster(from = may21rasterAddedMissing, to = emptyRaster)
may21rasterReproj <- projectRaster(from = may21raster, to = emptyRaster)

ecolBasicValid <- mapmisc::colourScale(values(may21rasterReproj), col = "Spectral", breaks = 10, rev = TRUE, style = "equal", dec = 1, opacity = 0.5)
jpeg("outputFiles/may21rasterOriginal.jpg", width = 1200, height = 1200)
raster::plot(may21rasterReproj, legend = FALSE, cex.axis = 2.5, col = ecolBasicValid$col, breaks = ecolBasicValid$breaks)
raster::plot(indiaPolygons, add = TRUE)
raster::plot(cityPoints, pch = 15, col = "red", cex = 5, add = TRUE)
text(x = cityPoints@coords[,1], y = cityPoints@coords[,2], labels = replace(citiesNames, 1, "Mumbai"), offset = 2, cex = 4, pos = 4)
mapmisc::legendBreaks("bottomleft", ecolBasicValid, title = "LST (Celsius)", bg = "white", cex = 2)
scaleBar(crs = crs(indiaPolygons), pos = "topleft", cex = 3, pt.cex = 2.2, title.cex = 3.5)
dev.off()

jpeg("outputFiles/may21rasterAddedMissing.jpg", width = 1200, height = 1200)
raster::plot(may21rasterAddedMissingReproj, legend = FALSE, cex.axis = 2.5, col = ecolBasicValid$col, breaks = ecolBasicValid$breaks)
raster::plot(indiaPolygons, add = TRUE)
raster::plot(cityPoints, pch = 15, col = "red", cex = 5, add = TRUE)
text(x = cityPoints@coords[,1], y = cityPoints@coords[,2], labels = replace(citiesNames, 1, "Mumbai"), offset = 2, cex = 4, pos = 4)
mapmisc::legendBreaks("bottomleft", ecolBasicValid, title = "LST (Celsius)", bg = "white", cex = 2)
scaleBar(crs = crs(indiaPolygons), pos = "topleft", cex = 3, pt.cex = 2.2, title.cex = 3.5)
dev.off()
######## END: Producing plot to represent validation test set ########

######## Checking validation error ########

temperaturesMay21 <- indiaTemperaturesSecondToLastWeek[["2012-05-21"]]
temperaturesMay21reproj <- raster(x = extent(testDataValidation@sp), nrows = nrow(temperaturesMay21), ncols = ncol(temperaturesMay21), crs = crs(testDataValidation@sp))
temperaturesMay21reproj <- projectRaster(from = temperaturesMay21, to = temperaturesMay21reproj)

recordedTemperaturesInMissingZone <- extract(x = temperaturesMay21reproj, y = testDataValidation@sp)

differences <- indiaAnalysisValidation$predMoments$Mean - recordedTemperaturesInMissingZone
MSPE <- mean(differences^2, na.rm = TRUE)
MedSPE <- median(differences^2, na.rm = TRUE)

quantile(abs(differences), probs = 0.9)

# What is the coordinate of the tile that produced the largest absolute difference?
testDataSpReproj <- spTransform(testDataValidation@sp, CRSobj = crs(indiaPolygons))
testDataSpReproj@coords[which.max(abs(differences)),]

spObjectReprojected <- sp::spTransform(testDataValidation@sp, CRSobj = crs(indiaPolygons))

spObjectReprojected@coords[which.min(indiaAnalysisValidation$predMoments$Mean), ]

fieldValues <- differences[!is.na(differences)]
sqDiffRaster <- raster::rasterize(x = spObjectReprojected, y = emptyRaster, field = fieldValues^2)
diffRaster <- raster::rasterize(x = spObjectReprojected, y = emptyRaster, field = fieldValues)

######## END: Checking validation error ########

####### SPDE results analysis ############

library(INLA)

trainingDataValidation@sp <- sp::spTransform(trainingDataValidation@sp, CRSobj = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
testDataValidation@sp <- sp::spTransform(testDataValidation@sp, CRSobj = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

timeVecTraining <- (as.numeric(time(trainingDataValidation)) - min(as.numeric(time(trainingDataValidation))))/(3600*24) + 1
timeVecTest <- (as.numeric(time(testDataValidation)) - min(as.numeric(time(trainingDataValidation))))/(3600*24) + 1
knots <- seq(1, max(timeVecTraining), length = max(timeVecTraining))
mesh1 <- inla.mesh.1d(loc = knots, degree = 2, boundary = "free")

## generate space mesh

mesh2 <- inla.mesh.2d(loc = trainingDataValidation@sp@coords[timeVecTraining == 1, ], cutoff = 0.01, offset = c(0.1, 0.2), max.n = 2000)

# range0 and sigma0 control the prior means for the range and scale parameters.
# See Lindgren INLA tutorial page 5.
d <- 1
alpha <- 2
kappa <- 1 # = 1/(range parameter in my model)
spatialSmoothness <- alpha - d/2 # cf p.3 INLA tutorial
loghyperparaSDinMyModel <- log(10)
# range0 and sigma0 seem to be the prior means...
range0 <- sqrt(8 * spatialSmoothness)/kappa # sqrt(8 * spatial smoothness) / Kappa. In my model, I use 1 as prior mean for spatial range and fix smoothness at 1.5. This means Kappa = 1.
sigma0 <- 1
lkappa0 <- log(8 * spatialSmoothness)/2 - log(range0)
ltau0 <- 0.5*log(gamma(spatialSmoothness)/(gamma(alpha)*(4*pi)^(d/2))) - log(sigma0) - spatialSmoothness * lkappa0

## build the spatial spde
spde <- inla.spde2.matern(mesh2, B.tau = matrix(c(ltau0, -1, spatialSmoothness), 1, 3),
                          B.kappa = matrix(c(lkappa0, 0, -1), 1,3),
                          theta.prior.mean = c(0,0), theta.prior.prec = c(1/loghyperparaSDinMyModel^2, 1/loghyperparaSDinMyModel^2))

## build the space time indices
STindex <- inla.spde.make.index("space", n.spde = spde$n.spde, n.group = mesh1$m)

## Link data and process

Atraining <- inla.spde.make.A(mesh2, loc = trainingDataValidation@sp@coords, group = timeVecTraining, group.mesh = mesh1)
Atest <- inla.spde.make.A(mesh2, loc = testDataValidation@sp@coords, group = timeVecTest, group.mesh = mesh1)

stackTraining <- inla.stack(data = list(y = trainingDataValidation@data$y), A = list(Atraining, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 , 1, 1, 1, 1),
 effects = list(
  c(STindex, list(intercept = 1)),
   list(landCover3 = trainingDataValidation@data$landCover2),
   list(landCover4 = trainingDataValidation@data$landCover4),
   list(landCover5 = trainingDataValidation@data$landCover5),
   list(landCover8 = trainingDataValidation@data$landCover8),
   list(landCover9 = trainingDataValidation@data$landCover9),
   list(landCover10 = trainingDataValidation@data$landCover10),
   list(landCover11 = trainingDataValidation@data$landCover11),
   list(landCover12 = trainingDataValidation@data$landCover12),
   list(landCover13 = trainingDataValidation@data$landCover13),
   list(landCover14 = trainingDataValidation@data$landCover14),
   list(landCover15 = trainingDataValidation@data$landCover15),
   list(elevation = trainingDataValidation@data$elevation),
   list(Aqua = trainingDataValidation@data$Aqua),
   list(time2 = trainingDataValidation@data$time2),
   list(time3 = trainingDataValidation@data$time3),
   list(time4 = trainingDataValidation@data$time4),
   list(time5 = trainingDataValidation@data$time5),
   list(time6 = trainingDataValidation@data$time6),
   list(time7 = trainingDataValidation@data$time7)
  ), tag="est")

stackTest <- inla.stack(
    data = list(y = NA),
    A = list(Atest, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 , 1, 1, 1, 1),
    effects = list(
     c(STindex, list(intercept = 1)),
      list(landCover3 = testDataValidation@data$landCover2),
      list(landCover4 = testDataValidation@data$landCover4),
      list(landCover5 = testDataValidation@data$landCover5),
      list(landCover8 = testDataValidation@data$landCover8),
      list(landCover9 = testDataValidation@data$landCover9),
      list(landCover10 = testDataValidation@data$landCover10),
      list(landCover11 = testDataValidation@data$landCover11),
      list(landCover12 = testDataValidation@data$landCover12),
      list(landCover13 = testDataValidation@data$landCover13),
      list(landCover14 = testDataValidation@data$landCover14),
      list(landCover15 = testDataValidation@data$landCover15),
      list(elevation = testDataValidation@data$elevation),
      list(Aqua = testDataValidation@data$Aqua),
      list(time2 = testDataValidation@data$time2),
      list(time3 = testDataValidation@data$time3),
      list(time4 = testDataValidation@data$time4),
      list(time5 = testDataValidation@data$time5),
      list(time6 = testDataValidation@data$time6),
      list(time7 = testDataValidation@data$time7)
     ),
    tag = 'predictions')

combinedStack <- inla.stack(stackTraining, stackTest)

stackIndex <- inla.stack.index(combinedStack, "predictions")$data
preds <- SPDEresult$summary.linear.predictor

differencesSPDE <- preds$mean[stackIndex] - recordedTemperaturesInMissingZone
MSPE_SPDE <- mean(differencesSPDE^2, na.rm = TRUE) # 4.94362
MedSPE_SPDE <- median(differencesSPDE^2, na.rm = TRUE) # 2.015934

smallRasterExtent <- raster::extent(c(xmin = min(smallWesternMahaPolygonEdges[ , 1]), xmax = max(smallWesternMahaPolygonEdges[ , 1]), ymin = min(smallWesternMahaPolygonEdges[ , 2]), ymax = max(smallWesternMahaPolygonEdges[ , 2])))

fieldValuesSPDE <- differencesSPDE[!is.na(differencesSPDE)]
spObjectReprojectedSPDE <- sp::spTransform(testDataValidation@sp[!is.na(differencesSPDE)], CRSobj = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
emptyRaster <- raster(x = smallRasterExtent, nrows = 120, ncols = 120, crs = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
diffRasterSPDE <- raster::rasterize(x = spObjectReprojectedSPDE, y = emptyRaster, field = fieldValuesSPDE)

citiesNames <- c("Mumbai", "Pune", "Nashik")
cityPoints <- mapmisc::geocode(x = c("Mumbai City", "Pune", "Nashik"))

ecolSPDE <- mapmisc::colourScale(values(diffRasterSPDE), col = "Spectral", breaks = 10, rev = TRUE, style = "equal", dec = 1, opacity = 0.5)
ecolSPDE <- ecolSPDE[c("col", "breaks")]

ecolINLAMRA <- mapmisc::colourScale(values(diffRaster), col = "Spectral", breaks = 10, rev = TRUE, style = "equal", dec = 1, opacity = 0.5)
ecolINLAMRA <- ecolINLAMRA[c("col", "breaks")]

newBoundaries <- range(c(ecolINLAMRA$breaks, ecolSPDE$breaks))
combinedBreaks <- seq(from = newBoundaries[[1]], to = newBoundaries[[2]], length.out = length(ecolINLAMRA$breaks))

jpeg(file = "outputFiles/predErrorOnMay21withTimeCovarSPDE.jpg", width = 1200, height = 1200)
  raster::plot(diffRasterSPDE, legend = FALSE, cex.axis = 2.5, col = ecolINLAMRA$col, breaks = combinedBreaks)
  plot(indiaPolygons, add = TRUE)
  plot(cityPoints, pch = 15, col = "red", cex = 5, add = TRUE)
  text(x = cityPoints@coords[,1], y = cityPoints@coords[,2], labels = citiesNames, offset = 2, cex = 5, pos = 4)
  mapmisc::legendBreaks("bottomleft", list(breaks = round(combinedBreaks, digits = 1), col = ecolINLAMRA$col), title = "Diff. (Celsius)", bg = "white", cex = 2)
  mapmisc::scaleBar(crs = crs(indiaPolygons), pos = "topleft", cex = 3, pt.cex = 2.2, title.cex = 3.5)
dev.off()

jpeg(file = "outputFiles/predErrorOnMay21withTimeCovar.jpg", width = 1200, height = 1200)
raster::plot(diffRaster, legend = FALSE, cex.axis = 2.5, col = ecolINLAMRA$col, breaks = combinedBreaks)
plot(indiaPolygons, add = TRUE)
mapmisc::legendBreaks("bottomleft", list(breaks = round(combinedBreaks, digits = 1), col = ecolINLAMRA$col), title = "Diff. (Celsius)", bg = "white", cex = 2)
# raster::plot(diffRaster, legend.only = TRUE, legend.width = 4, axis.args = list(cex.axis = 3))
plot(cityPoints, pch = 15, col = "red", cex = 5, add = TRUE)
text(x = cityPoints@coords[,1], y = cityPoints@coords[,2], labels = replace(citiesNames, 1, "Mumbai"), offset = 2, cex = 5, pos = 4)
mapmisc::scaleBar(crs = crs(indiaPolygons), pos = "topleft", cex = 3, pt.cex = 2.2, title.cex = 3.5)
dev.off()

####### END: SPDE results analysis ############

######## Results for main analysis ########

spAndTimeCorr <- lapply(list(space = list(paraName = c("space.rho", "space.smoothness"), distances = c(close = 1, far = 10)), time = list(paraName = c("time.rho", "time.smoothness"), distances = c(close = 1, far = 7))), function(component) {
  sapply(component$distances, function(distanceMeasure) {
    maternCov(d = distanceMeasure, rho = exp(indiaAnalysisMain$hyperMarginalMoments[component$paraName[[1]], "Mean"]), smoothness = exp(indiaAnalysisMain$hyperMarginalMoments[component$paraName[[2]], "Mean"]), scale = 1)
  })
})

# Latitude was excluded because we did not have enough variation...
predictionDataMain@data <- subset(predictionDataMain@data, select = -latitude)

######### Preparing parameter and hyperparameter moments table ########

hyperMoments <- subset(indiaAnalysisMain$hyperMarginalMoments[c("space.rho", "time.rho", "scale"), ], select = -Skewness)

covariateMoments <- indiaAnalysisMain$FEmarginalMoments

combinedMoments <- rbind(hyperMoments, covariateMoments)

xtable::xtable(combinedMoments, caption = "Mean and standard deviation of hyperparameters and fixed effects posteriors", digits = c(0, 3, 3, 3, 3))
######## END: Preparing parameter and hyperparameter moments table ########

table((indiaAnalysisMain$predMoments$Mean - min(trainingDataMain@data[, "y"])) < -1)
extremeValueIndices <- which(indiaAnalysisMain$predMoments$Mean < (min(trainingDataMain@data[, "y"]) - 1))
mostExtremeIndex <- which.min(indiaAnalysisMain$predMoments$Mean)
predictionDataMain@data[mostExtremeIndex, ]
predictionOrder <- order(indiaAnalysisMain$predMoments$Mean)
predictionDataMain@sp@coords[predictionOrder[1:2],]

extremeDataset <- predictionDataMain[extremeValueIndices]
extremeDataset@data <- data.frame(Temperature = indiaAnalysisMain$predMoments$Mean[extremeValueIndices])

jpeg("outputFiles/extremeValuePositions.jpeg", width = 1000, height = 1000)
stplot(extremeDataset)
dev.off()

######## Plotting predictions ########

indiaPolygons <- raster::getData(country = "IND", level = 2)
trainingDataMainReproject <- trainingDataMain
trainingDataMainReproject@sp <- sp::spTransform(trainingDataMain@sp, crs(indiaPolygons))
predictionDataMainReproject <- predictionDataMain
predictionDataMainReproject@sp <- sp::spTransform(predictionDataMain@sp, crs(indiaPolygons))

plot(indiaAnalysisMain,
   filename = "outputFiles/predictionsIndiaDataMay28mainAnalysisJoint_Larger.jpg",
   type = "joint",
   polygonsToOverlay = indiaPolygons,
   control =  list(
                graphicsEngine = jpeg,
                controlForScaleBar = list(
                  pos = "topleft",
                  cex = 2,
                  pt.cex = 1.5
                ),
                controlForRasterLegend = list(
                  pos = "bottomleft",
                  title = "LST (Celsius)",
                  bg = "white",
                  cex = 2
                ),
                controlForRasterColourScale = list(
                  col = "Spectral",
                  breaks = 10,
                  rev = TRUE,
                  style = "equal",
                  dec = 1,
                  opacity = 0.5
                ),
                controlForRasterPlot = list(
                  cex.axis = 2,
                  cex.main = 3
                ),
                resolutionInMeters = 1000,
                trim = 2,
                timesToPlot = unique(time(predictionDataMain))
              ),
   width = 1100,
   height = 1200
)

plot(indiaAnalysisMain,
   filename = "outputFiles/predictionsIndiaDataMay28mainAnalysisSDs_Larger.jpg",
   type = "SD",
   polygonsToOverlay = indiaPolygons,
   control =  list(
                graphicsEngine = jpeg,
                controlForScaleBar = list(
                  pos = "topleft",
                  cex = 2,
                  pt.cex = 1.5
                ),
                controlForRasterLegend = list(
                  pos = "bottomleft",
                  title = "LST (Celsius)",
                  bg = "white",
                  cex = 2
                ),
                controlForRasterColourScale = list(
                  col = "Spectral",
                  breaks = 10,
                  rev = TRUE,
                  style = "equal",
                  dec = 1,
                  opacity = 0.5
                ),
                controlForRasterPlot = list(
                  cex.axis = 2,
                  cex.main = 3
                ),
                resolutionInMeters = 1000,
                timesToPlot = unique(time(predictionDataMain))
              ),
   width = 1100,
   height = 1200
)

plot(indiaAnalysisMain,
   filename = "outputFiles/predictionsIndiaDataMay28mainAnalysisTraining_Larger.jpg",
   type = "training",
   polygonsToOverlay = indiaPolygons,
   control =  list(
                graphicsEngine = jpeg,
                controlForScaleBar = list(
                  pos = "topleft",
                  cex = 2,
                  pt.cex = 1.5
                ),
                controlForRasterLegend = list(
                  pos = "bottomleft",
                  title = "LST (Celsius)",
                  bg = "white",
                  cex = 2
                ),
                controlForRasterColourScale = list(
                  col = "Spectral",
                  breaks = 10,
                  rev = TRUE,
                  style = "equal",
                  dec = 1,
                  opacity = 0.5
                ),
                controlForRasterPlot = list(
                  cex.axis = 2,
                  cex.main = 3
                ),
                resolutionInMeters = 1000,
                timesToPlot = unique(time(predictionDataMain))
              ),
   width = 1100,
   height = 1200
)

######## END: Plotting predictions ########
