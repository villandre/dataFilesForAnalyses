# dataFilesForAnalyses
Data and script files to perform analyses in studies submitted for publication. For now, the repository
contains material for the 2020 paper introducing INLA-MRA.

Prior to running the analyses, a working directory must be set up, and its name must be entered in the call to setwd in IndiaDataMainAnalysis.R and IndiaDataValidationAnalysis.R. The code should work straight away as long as all required libraries are installed and **subfolders "data" and "outputFiles" are created within the working directory.**

The script files assume that **data files are in a subfolder called "data"**. This can of course be changed by modifying the folder name in the two calls to "load".  

The INLAMRA function will save intermediate results (to allow for interruptions/restarts) in a subfolder called "outputFiles". This behaviour can be disabled by commenting out the lines starting with "fileToSaveOptOutput =" and "folderToSaveISpoints =".

**The final output will also be saved in the "outputFiles" subfolder.**

The script in graphsAndSummariesForJASApaper.R should produce all the graphs and the table in the paper.
