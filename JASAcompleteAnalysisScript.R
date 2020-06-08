# IMPORTANT: Read instructions at https://github.com/villandre/dataFilesForAnalyses/blob/master/README.md to learn how to set up the necessary directory structure and make the necessary edits to the script files. The most important step is to modify the working directories specified in *all the scripts*. The scripts were written to work independently one from the another.
# 
# Running the script could take several days. We recommend running it in BATCH mode, i.e. input
# 
# R CMD BATCH --no-restore --no-save JASAcompleteAnalysisScript.R &
# 
# in your console. On a remote server, running it in a terminal multiplexer like 'tmux' could also help ensure the code runs to completion.  
# 
# We ran all simulations in Ubuntu Bionic, but the code should work properly under any Linux system with OpenMP support. 

# **CHANGE THE FOLLOWING WORKING DIRECTORY.**
 
setwd("/home/luc/INLAMRAfiles/INLAMRApaper1/realData")

# **NOTE: THE CODE FILES ARE ASSUMED TO BE IN A SUBFOLDER CALLED "code".** Modify this if it is not the case.

################################
# This section is OPTIONAL. We create here all the datasets that will be used in the upcoming analyses, which are already available as Rdata objects on Github. This step requires downloading the LST, land cover, and elevation data from EarthData using the MODISdownload.sh script beforehand. Connecting to the EarthData server requires creating a profile with EarthData, which can be done for free. By default, the files have to be stored in the "data" subfolder, but this can be changed by modifying variable "rawDataFilesLocation" in the scripts.

rm(list = ls()) # To prevent potential interactions...
source("code/prepareDataMainAnalysis.R")
rm(list = ls()) # To prevent potential interactions...
source("code/prepareDataValidationAnalysis.R")
################################

# We start with the main analysis.
# WARNING: The following script will only run properly on a machine with a large memory bank, i.e. > 100 GB. **A memory shortage might result in a segmentation fault/bad write, which could cause R to crash**.

rm(list = ls()) # To prevent potential interactions...
source("code/IndiaDataMainAnalysis.R")

# We then run the validation analysis.

rm(list = ls()) # To prevent potential interactions...
source("code/IndiaDataValidationAnalysis.R")

# Finally, we produce all summaries/statistics presented in the paper 
# REMINDER: This step will also require running the MODISdownload.sh script beforehand.

rm(list = ls()) # To prevent potential interactions...
source("code/graphsAndSummariesForJASApaper.R", echo = TRUE)