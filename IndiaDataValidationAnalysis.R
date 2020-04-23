# Script to run the validation task. We suggest running it in BATCH mode.

library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(spacetime)
library(MRAinla)
library(RhpcBLASctl)

blas_set_num_threads(1)
omp_set_num_threads(1)

 # Edit the following line for your system
 # Make sure to create an "outpufFiles" subfolder.
 # Put the data files in a "data" subfolder.

setwd("/home/luc/INLAMRAfiles/INLAMRApaper1/realData")

## Importing data

# The data files should be in subfolder data/.
validationTestDataName <- load("data/testDataMay21_May18_24.Rdata")
validationTestData <- get(validationTestDataName)
rm(ls = validationTestDataName)

validationTrainingDataName <- load("data/mainDataCompleteMap_May18_24.Rdata")
validationTrainingData <- get(validationTrainingDataName)
rm(ls = validationTrainingDataName)

validationTrainingData@sp <- spTransform(x = validationTrainingData@sp, CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
validationTestData@sp <- spTransform(x = validationTestData@sp, CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

## Not enough latitude variation to make it worth adjusting for it...

validationTrainingData@data <- subset(validationTrainingData@data, select = -latitude)
validationTestData@data <- subset(validationTestData@data, select = -latitude)

mainDataMeanElevation <- mean(validationTrainingData@data$elevation)
validationTrainingData@data$elevation <- validationTrainingData@data$elevation - mainDataMeanElevation

validationTestData@data$elevation <- validationTestData@data$elevation - mainDataMeanElevation

fixedEffSD <- 10
errorSD <- 0.5 # Based on https://landval.gsfc.nasa.gov/Results.php?TitleID=mod11_valsup10

# Log-scale:
hyperStart <- list(
  space = c(rho = 0),
  time = c(rho = 0),
  scale = 0)

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
    rho = c(mu = hyperStart$space["rho"], sigma = logHyperpriorSD)),
  time = list(
    smoothness = c(mu = fixedHyperValues$time[["smoothness"]], sigma = logHyperpriorSD),
    rho = c(mu = hyperStart$time["rho"], sigma = logHyperpriorSD)),
  scale = c(mu = hyperStart$scale, sigma = logHyperpriorSD * 2),
  errorSD = c(mu = fixedHyperValues$errorSD , sigma = logHyperpriorSD),
  fixedEffSD = c(mu = fixedHyperValues$fixedEffSD, sigma = logHyperpriorSD))

set.seed(10)

indiaAnalysisValidation <- INLAMRA(
  responseVec = validationTrainingData@data[ , "y"],
  covariateFrame = subset(validationTrainingData@data, select = -y),
  spatialCoordMat = validationTrainingData@sp@coords,
  timePOSIXorNumericVec = time(validationTrainingData),
  predCovariateFrame = validationTestData@data,
  predSpatialCoordMat = validationTestData@sp@coords,
  predTimePOSIXorNumericVec = time(validationTestData),
  spatialRangeList = list(start = hyperStart$space[["rho"]], hyperpars = hyperNormalList$space$rho),
  spatialSmoothnessList = list(start = fixedHyperValues$space[["smoothness"]]),
  timeRangeList = list(start = hyperStart$time[["rho"]], hyperpars = hyperNormalList$time$rho),
  timeSmoothnessList = list(start = fixedHyperValues$time[["smoothness"]]),
  scaleList = list(start = hyperStart$scale, hyperpars = hyperNormalList$scale),
  errorSDlist = list(start = fixedHyperValues$errorSD),
  fixedEffSDlist = list(start = fixedHyperValues$fixedEffSD),
   control = list(
     Mlon = 4,
     Mlat = 5,
     Mtime = 0,
     numValuesForIS = 100,
     numKnotsRes0 = 8,
     numIterOptim = 20,
     numOpenMPthreads = 12L,
     #fileToSaveOptOutput = "outputFiles/optimOutputValidation.Rdata", # Not essential, but convenient
     #folderToSaveISpoints = "outputFiles/ISpointsValidation", # Not essential, but convenient
     tipKnotsThinningRate = 0.5,
     spaceJitterMax = 0, # Already jittered in spatial coordinates.
     timeJitterMaxInDecimalDays = 0 # No need to jitter time.
   )
)

# The working directory should have subfolder /outputFiles.

save(indiaAnalysisValidation, file = "outputFiles/INLAMRA_validationAnalysis.Rdata", compress = TRUE)
