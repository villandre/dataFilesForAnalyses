# IMPORTANT: Read instructions at https://github.com/villandre/dataFilesForAnalyses/blob/master/README.md to learn how to set up the necessary directory structure and make the necessary edits to the script files.

# **CHANGE THE FOLLOWING WORKING DIRECTORY.**
 
setwd("/home/luc/INLAMRAfiles/INLAMRApaper1/realData")

# Running the following script could take several days. We recommend running it in BATCH mode, i.e. input
# 
# R CMD BATCH --no-restore --no-save JASAcompleteAnalysisScript.R &
# 
# in your console. On a remote server, running it in a terminal multiplexer like 'tmux' could also help ensure the code runs to completion.  
# 
# We ran all simulations in Ubuntu Bionic, but the code should work properly under any Linux system with OpenMP support. 

# We start with the main analysis (WARNING: The following code will only run properly on a machine with a large memory bank, i.e. > 100 GB. A memory shortage might result in a segmentation fault/bad write.)

# **NOTE THAT THE CODE FILES ARE ASSUMED TO BE IN A SUBFOLDER CALLED "code".** Modify this if it is not the case.

source("code/IndiaDataMainAnalysis.R")

# We then run the validation analysis.

rm(list = ls()) # To prevent potential interactions...

source("code/IndiaDataValidationAnalysis.R")

# Finally, we produce all summaries/statistics presented in the paper (REMINDER: This step will require running the MODISdownload.sh script beforehand.)

rm(list = ls()) # To prevent potential interactions...

source("code/graphsAndSummariesForJASApaper.R", echo = TRUE)