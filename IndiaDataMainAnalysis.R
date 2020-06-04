# Script to perform the main analysis. We suggest running it in BATCH mode.
 
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(spacetime)
library(MRAinla)
library(RhpcBLASctl)

# This ensures that no unwanted parallelisation occurs.

blas_set_num_threads(1)
omp_set_num_threads(1)

# The working directory should be changed before running the code on another machine.

setwd("/home/luc/INLAMRAfiles/INLAMRApaper1/realData")

## Loading datasets for analyses

mainDataName <- load("data/mainDataCompleteMap_May25_31_Larger.Rdata")
mainData <- get(mainDataName)
rm(ls = mainDataName)

predDataName <- load("data/testDataMay28_May25_31_Larger.Rdata")
predData <- get(predDataName)
rm(ls = predDataName)

## Switching to the standard longitude/latitude projection

mainData@sp <- spTransform(x = mainData@sp, CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
predData@sp <- spTransform(x = predData@sp, CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

## Not enough latitude variation to make it worth adjusting for it: we remove the covariate.

mainData@data <- subset(mainData@data, select = -latitude)
predData@data <- subset(predData@data, select = -latitude)

## We center the elevation covariate.

mainDataMeanElevation <- mean(mainData@data$elevation)
mainData@data$elevation <- mainData@data$elevation - mainDataMeanElevation

predData@data$elevation <- predData@data$elevation - mainDataMeanElevation

## We define starting values for hyperparameters (on the logarithmic scale).

hyperStart <- list(
  space = c(rho = 0),
  time = c(rho = 0),
  scale = 0)

## The following hyperparameter values are fixed.

errorSD <- 0.5 # Based on https://landval.gsfc.nasa.gov/Results.php?TitleID=mod11_valsup10
fixedEffSD <- 10

fixedHyperValues <- list(
  space = c(smoothness = log(1.5)),
  time = c(smoothness = log(0.5)),
  errorSD = log(errorSD),
  fixedEffSD = log(fixedEffSD)
)

logHyperpriorSD <- 2

## These are the values for the mean and standard deviation parameters for the normal hyperpriors.

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

## We now run the main analysis.

indiaAnalysis <- INLAMRA(
  responseVec = mainData@data[ , "y"],
  covariateFrame = subset(mainData@data, select = -y), # The covariate frame should not include the response variable: we remove it.
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
    numOpenMPthreads = 12L,
    fileToSaveOptOutput = "outputFiles/optimOutputMay28evenLarger.Rdata", # Comment out this line and the next if you don't want intermediate results saved to the hard drive (these intermediate results allow the user to stop and resume the function.) 
    folderToSaveISpoints = "outputFiles/ISpointsMay28evenLarger",
    tipKnotsThinningRate = 1/3,
    spaceJitterMax = 0, # Already jittered in spatial coordinates.
    timeJitterMaxInDecimalDays = 0 # No need to jitter time.
  )
)

## Results are saved in subdirectory outputFiles. Please ensure that it exists or change the file argument.

save(indiaAnalysis, file = "outputFiles/INLAMRA_mainAnalysis.Rdata", compress = TRUE)
