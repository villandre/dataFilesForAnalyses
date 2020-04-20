library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(spacetime)
library(MRAinla)
library(RhpcBLASctl)

blas_set_num_threads(1)
omp_set_num_threads(8)

setwd("/home/luc/INLAMRAfiles/INLAMRApaper1/realData")

## Setting up parameters

mainDataName <- load("data/mainDataCompleteMap_May25_31_Larger.Rdata")
mainData <- get(mainDataName)
rm(ls = mainDataName)

predDataName <- load("data/testDataMay28_May25_31_Larger.Rdata")
predData <- get(predDataName)
rm(ls = predDataName)

mainData@sp <- spTransform(x = mainData@sp, CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
predData@sp <- spTransform(x = predData@sp, CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

## Not enough latitude variation to make it worth adjusting for it...

mainData@data <- subset(mainData@data, select = -latitude)
predData@data <- subset(predData@data, select = -latitude)

mainDataMeanElevation <- mean(mainData@data$elevation)
mainData@data$elevation <- mainData@data$elevation - mainDataMeanElevation

predData@data$elevation <- predData@data$elevation - mainDataMeanElevation

##################################
# Log-scale:
hyperStart <- list(
  space = c(rho = 0),
  time = c(rho = 0),
  scale = 0)

errorSD <- 0.5 # Based on https://landval.gsfc.nasa.gov/Results.php?TitleID=mod11_valsup10
fixedEffSD <- 10

fixedHyperValues <- list(
  space = c(smoothness = log(1.5)),
  time = c(smoothness = log(0.5)),
  errorSD = log(errorSD),
  fixedEffSD = log(fixedEffSD)
)

logHyperpriorSD <- 2

hyperNormalList <- list(
  space = list(
    smoothness = c(mu = fixedHyperValues$space[["smoothness"]], sigma = logHyperpriorSD),
    rho = c(mu = hyperStart$space[["rho"]], sigma = logHyperpriorSD)),
  time = list(
    smoothness = c(mu = fixedHyperValues$time[["smoothness"]], sigma = logHyperpriorSD),
    rho = c(mu = hyperStart$time[["rho"]], sigma = logHyperpriorSD)),
  scale = c(mu = hyperStart$scale, sigma = logHyperpriorSD * 2), # Prior should be more vague, see Lindgren INLA tutorial p. 12
  errorSD = c(mu = fixedHyperValues$errorSD , sigma = logHyperpriorSD),
  fixedEffSD = c(mu = fixedHyperValues$fixedEffSD, sigma = logHyperpriorSD)
)

indiaAnalysis <- INLAMRA(
  responseVec = mainData@data[ , "y"],
  covariateFrame = subset(mainData@data, select = -y),
  spatialCoordMat = mainData@sp@coords,
  timePOSIXorNumericVec = time(mainData),
  predCovariateFrame = predData@data,
  predSpatialCoordMat = predData@sp@coords,
  predTimePOSIXorNumericVec = time(predData),
  spatialRangeList = list(start = hyperStart$space[["rho"]], hyperpars = hyperNormalList$space$rho),
  spatialSmoothnessList = list(start = fixedHyperValues$space[["smoothness"]]),
  timeRangeList = list(start = hyperStart$time[["rho"]], hyperpars = hyperNormalList$time$rho),
  timeSmoothnessList = list(start = fixedHyperValues$time[["smoothness"]]),
  scaleList = list(start = hyperStart$scale, hyperpars = hyperNormalList$scale),
  errorSDlist = list(start = fixedHyperValues$errorSD),
  fixedEffSDlist = list(start = fixedHyperValues$fixedEffSD),
  control = list(
    Mlon = 6,
    Mlat = 5,
    Mtime = 1,
    numKnotsRes0 = 8,
    numIterOptim = 20,
    fileToSaveOptOutput = "outputFiles/optimOutputMay28evenLarger.Rdata",
    folderToSaveISpoints = "outputFiles/ISpointsMay28evenLarger",
    tipKnotsThinningRate = 1/3
  )
)

save(indiaAnalysis, file = "outputFiles/INLAMRA_mainAnalysisResults.Rdata", compress = TRUE)
