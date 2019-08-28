
source("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/mFam Aggregation/CommunityToMassBank.R")

parentFolder <- "/mnt/data/IPB/Projects/2017_005_MS-databases/mFam contributions/"
aggregationFolder <- paste(parentFolder, "mFam Aggregation", sep = "")

folder <- "/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/MyProject"
ms_type <- "MS2"
accessionPrefix <- "[TODO]"
xlsxFile <- "[TODO]"
takeRecordedNames <- TRUE
applyIntensityThreshold <- FALSE

preprocessContributorToMassBankWorkflow(folder, accessionPrefix, xlsxFile, ms_type, takeRecordedNames)
runContributorToMassBankWorkflow(folder, applyIntensityThreshold, reprocess = FALSE)

## utility
#cleanContribtorDirectoryForReprocessing(folder)
#correctRecords(folder)

#aggregateSpectra_all(parentFolder, aggregationFolder)
#aggregateSpectra_ready(processedFolders, aggregationFolder, tag = "Validated")
