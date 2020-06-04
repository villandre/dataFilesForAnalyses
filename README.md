# dataFilesForAnalyses
Data and script files to perform analyses in studies submitted for publication. For now, the repository
contains material for the 2020 paper introducing INLA-MRA.

Prior to running the analyses, a working directory must be set up, and its name must be entered in the call to setwd in IndiaDataMainAnalysis.R and IndiaDataValidationAnalysis.R. The code should work straight away as long as all required libraries are installed and **subfolders "data" and "outputFiles" are created within the working directory.**

The script files assume that **data files are in a subfolder called "data"**. This can of course be changed by modifying the folder name in the two calls to "load".  

Installing the MRAinla R library is best done by running (in the R console)

remotes::install_github(repo = "https://github.com/villandre/MRAinla")

The calls to the INLAMRA function in IndiaDataMainAnalysis.R and IndiaDataValidationAnalysis.R will save intermediate results (to allow for interruptions/restarts) in a subfolder called "outputFiles". This behaviour can be disabled by commenting out the lines starting with "fileToSaveOptOutput =" and "folderToSaveISpoints =".

**The final output will also be saved in the "outputFiles" subfolder**, so make sure that the subfolder exists or that the line has been edited to refer to an existing directory. 

The script in graphsAndSummariesForJASApaper.R should produce all the graphs, tables, and summary statistics presented in the paper.

The code has been tested in Ubuntu Bionic (18.04). With suitable edits, it should work properly in other Linux distributions and macOS. Unfortunately, we have been unable to test it in Windows.
