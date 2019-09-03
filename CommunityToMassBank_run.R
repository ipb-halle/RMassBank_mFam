#!/usr/bin/env Rscript
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

source("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/mFam Aggregation/CommunityToMassBank.R")

# Taking the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Load required libraries
library(stringr)
library(stringi)

# Set files and folders
parentFolder <- "/vol/metfamily/massbank/"
aggregationFolder <- paste(parentFolder, "mFam Aggregation", sep = "")

folder <- "/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam contributions/MyProject"
ms_type <- "MS2"
accessionPrefix <- "[TODO]"
xlsxFile <- "[TODO]"
takeRecordedNames <- TRUE
applyIntensityThreshold <- FALSE

# Run
preprocessContributorToMassBankWorkflow(folder, accessionPrefix, xlsxFile, ms_type, takeRecordedNames)
runContributorToMassBankWorkflow(folder, applyIntensityThreshold, reprocess = FALSE)

# Utility
#cleanContribtorDirectoryForReprocessing(folder)
#correctRecords(folder)

#aggregateSpectra_all(parentFolder, aggregationFolder)
#aggregateSpectra_ready(processedFolders, aggregationFolder, tag = "Validated")
