# dataFilesForAnalyses
Data and script files to perform analyses in studies submitted for publication. For now, the repository
contains material for the 2020 paper introducing INLA-MRA.

Prior to running any of the analyses, a working directory must be set up, and its name must be entered in the calls to `setwd` in `IndiaDataMainAnalysis.R`, `IndiaDataValidationAnalysis.R`, `graphsAndSummariesForJASApaper.R` and `JASAcompleteAnalysisScript.R` (as well as the `prepareData*.R` files, should they be required). The code should work straight away as long as all required libraries are installed and **subfolders `data`, `code`, and `outputFiles` are created within the working directory.**

The scripts assume that **all data files are in a subfolder called `data`**. This can be changed by modifying the folder name in the two calls to `load`, and the `rawDataFilesLocation` variable in the `prepare*.R` and `graphsAndSummariesForJASApaper.R` files.  

Installing the `MRAinla` library is best done by calling (in the R console)

remotes::install_github(repo = "https://github.com/villandre/MRAinla")

The calls to the `INLAMRA` function in `IndiaDataMainAnalysis.R` and `IndiaDataValidationAnalysis.R` will save intermediate results (to allow for interruptions/restarts) in a subfolder called `outputFiles`. This behaviour can be disabled by commenting out the lines starting with `fileToSaveOptOutput =` and `folderToSaveISpoints =`.

**The final output will also be saved in the `outputFiles` subfolder**, so make sure that the subfolder exists or that the lines with the calls to `save` have been edited to refer to an existing directory. 

The script in `graphsAndSummariesForJASApaper.R` should produce all the graphs, tables, and summary statistics presented in the paper.

The code has been tested in Ubuntu Bionic (18.04), and should therefore work properly in other Linux distributions and macOS.
