#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load required libraries
library(stringr)
library(stringi)
library(rinchi)
library(RMassBank)
#library("RMassBank", lib.loc="/home/kpeters/R/x86_64-pc-linux-gnu-library/3.5/")
#source("/tmp/RMassBank/R/leCsvAccess.R")
#source("/tmp/RMassBank/R/formulaCalculator.R")
source("/usr/src/RMassBank/R/formulaCalculator.R")
source("/usr/src/RMassBank/R/leCsvAccess.R")
library(readxl)
library(webchem)

# Options
options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

# Taking the command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Take in trailing command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    print("Error! No or not enough arguments given.")
    print("Usage: $0 accessionPrefix xlsxFile")
    quit(save="no", status=1, runLast=FALSE)
}

# Set files and folders
#parentFolder <- "/home/kpeters/msp2massbank"
parentFolder <- "/tmp/msp2massbank"
aggregationFolder <- paste(parentFolder, "mFam Aggregation", sep = "")

#folder <- "/home/kpeters/msp2massbank"
folder <- parentFolder
ms_type <- "MS2"
accessionPrefix <- as.character(args[1]) #"IH" # "[TODO]" #as.character(args[1])
xlsxFile <- as.character(args[2])
takeRecordedNames <- TRUE
applyIntensityThreshold <- FALSE

folderMsp        <- paste(folder, "converted to msp", sep = "/")
folderMassBank   <- paste(folder, "converted to MassBank", sep = "/")
folderMetaData   <- paste(folder, "meta data", sep = "/")
folderRawData    <- paste(folder, "raw data", sep = "/")
folderRMassBankData <- paste(folder, "RMassBank", sep = "/")







## accession praefixes:
## [AU, BM, BS, CA, CE, CO, EA, EQ, ET, FF, FI, FU, GL, JE, JP, KN, KO, KZ, MC, ML, MS, MT, NU, OU, PB, PR, TT, TY, UA, UF, UO, UT, WA]
knownAccessionPrefixes <- c(
  "AC", "AU", "BM", "BML", "BS", "BSU", "CA", "CE", "CO", "EA", "EQ", "ET", "FF", "FFF", "FI", "FIO", "FU", "GL", "GLS", "HB", "JE", "JEL", "JP", 
  "KN", "KNA", "KO", "KZ", "LIT", "MC", "MCH", "ML", "MS", "MSJ", "MT", "NA", "NU", "OU", "OUF", "PB", "PR", "RP", "SM", "TT", "TUE", "TY", "UA", 
  "UF", "UN", "UO", "UP", "UPA", "UT", "WA"
)

knownAdducts <- getAdductInformation("")$adductString# c("[M+H]+", "[M+Na]+", "[M+K]+", "[M]+", "[M+NH4]+", "[M+ACN+H]+", "[M+ACN+Na]+", "[M+2Na-H]+", "[2M+H]+", "[2M+K]+", "[2M+Na]+", "[2M+NH4]+", "[2M+ACN+H]+", "[M+ACN+2H]2+", "[M+2H]2+", "[M-H]-", "[M]-", "[M+FA]-") ## --> adductToIonMode
adductOverview <- function(){
  library("RMassBank")
  adductDf <- getAdductInformation("C100")
  massAdditions <- sapply(X = adductDf$addition, FUN = RMassBank:::getMonoisotopicMass)
  adductDf$massAddition <- massAdditions
  adductDf <- adductDf[adductDf$massAddition < 1000,]
  adductDf <- adductDf[order(adductDf$massAddition), ]
  
  print(adductDf)
  
  deltaMass <- 219.1743 - 234.162
  print(adductDf[abs(adductDf$massAddition - deltaMass) < 0.5, ])
}

preprocessData_compoundInformation <- function(metaDataCompoundsDf){
  ##################################################################################
  ## collect compound information
  cat("\nFetching structures and PubChem CIDs...")
  
  cat("pre...")
  smiles     <- metaDataCompoundsDf$SMILES
  inchis     <- metaDataCompoundsDf$InChI
  pubchemIDs <- metaDataCompoundsDf$`PubChem CID`
  
  pubchemIDthere <- grepl(x = pubchemIDs, pattern = "^\\d+$")
  pubchemIDs[!pubchemIDthere] <- NA
  
  smiles     <- smiles[!pubchemIDthere]
  inchis     <- inchis[!pubchemIDthere]
  
  smiles     <- unique(smiles    [smiles     != ""])
  inchis     <- unique(inchis    [inchis     != ""])
  #pubchemIDs <- unique(pubchemIDs[pubchemIDs != ""])
  inchiKeys  <- vector(mode = "character", length = 0)
  
  inchiKeyRegExPattern <- "^[A-Z]{14,14}-[A-Z]{10,10}-[A-Z]$"
  isInchiKey <- grepl(x = inchis, pattern = inchiKeyRegExPattern)
  if(any(isInchiKey)){
    inchiKeys <- inchis[isInchiKey]
    inchis <- inchis[!isInchiKey]
  }
  
  smiles     <- smiles    [!(is.na(smiles    ) | smiles     == "NA")]
  inchis     <- inchis    [!(is.na(inchis    ) | inchis     == "NA")]
  #pubchemIDs <- pubchemIDs[!(is.na(pubchemIDs) | pubchemIDs == "NA")]
  inchiKeys  <- inchiKeys [!(is.na(inchiKeys ) | inchiKeys  == "NA")]
  
  ## remove inchiKeys
  #metaDataCompoundsDf$InChI[metaDataCompoundsDf$InChI %in% inchiKeys] <- ""
  
  ## smiles / inchiKey to inchi
  cat("to inchi...")
  inchiKeyToInchi <- vector(mode = "character", length = 0)
  if(length(inchiKeys) > 0)
    #inchiKeyToInchi <- sapply(X = inchiKeys,  FUN = function(inchiKey){ webchem::cs_inchikey_inchi(inchikey = inchiKey, verbose = FALSE) })
    inchiKeyToInchi <- sapply(X = inchiKeys,  FUN = function(inchiKey){ tryCatch({webchem::cs_inchikey_inchi(inchikey = inchiKey, verbose = FALSE)}, warning = function(w) {ifelse(test = w$message=="inchikey not found... Returning NA.", yes = NA, no = warning(w))}) })
  smilesToInchi <- vector(mode = "character", length = 0)
  if(length(smiles) > 0)
    smilesToInchi   <- sapply(X = smiles,     FUN = function(smile   ){ tryCatch({webchem::cs_smiles_inchi(  smiles = smile,      verbose = FALSE)}, warning = function(w) {ifelse(test = w$message=="inchi not found... Returning NA.",    yes = NA, no = warning(w))}) })
  
  isNAtmp          <- is.na(inchiKeyToInchi) | inchiKeyToInchi == ""
  inchiKeyToInchi  <- inchiKeyToInchi[!isNAtmp]
  inchiKeys2       <- inchiKeys[!isNAtmp]
  
  isNAtmp          <- is.na(smilesToInchi) | smilesToInchi == ""
  smilesToInchi    <- smilesToInchi[!isNAtmp]
  smiles2          <- smiles[!isNAtmp]
  
  ## record results
  potentialInchiKeysFromTable <- metaDataCompoundsDf$InChI
  #metaDataCompoundsDf$InChI        [match(x = inchiKeys2, table = metaDataCompoundsDf$InChI )] <- inchiKeyToInchi
  #metaDataCompoundsDf$InChI        [match(x = smiles2,    table = metaDataCompoundsDf$SMILES, )] <- smilesToInchi
  for(idx in seq_along(inchiKeys2))
    metaDataCompoundsDf$InChI        [metaDataCompoundsDf$InChI  %in% inchiKeys2[idx]] <- inchiKeyToInchi[[idx]]
  for(idx in seq_along(smiles2))
    metaDataCompoundsDf$InChI        [metaDataCompoundsDf$SMILES %in% smiles2[idx]   ] <- smilesToInchi  [[idx]]
  
  ## inchi to PubChem
  cat("to CID")
  allInchis <- unique(c(inchis, smilesToInchi, inchiKeyToInchi))
  allInchis <- allInchis[allInchis %in% metaDataCompoundsDf$InChI]
  #inchiToPubchemId_webchem <- vector(mode = "character", length = 0)
  #inchiToPubchemId_webchem2 <- vector(mode = "character", length = 0)
  inchiToPubchemId_pubchem <- vector(mode = "character", length = 0)
  inchiToPubchemId_pubchem2 <- vector(mode = "character", length = 0)
  inchiToPubchemId_pubchem3 <- vector(mode = "character", length = 0)
  
  library("jsonlite")
  cat(".")
  if(length(inchiKeys) > 0){
    #inchiToPubchemId_webchem2 <- sapply(X = inchiKeys, FUN = function(inchiKey){webchem::cs_inchikey_csid(inchikey = inchiKey, verbose = FALSE)})
    inchiToPubchemId_pubchem2 <- sapply(X = inchiKeys, FUN = function(inchiKey){
      
      ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON?inchi=InChI=1S/C3H8/c1-3-2/h3H2,1-2H3
      ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/AADVZSXPNRLYLV/JSON
      #url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/", inchiKey, "/JSON", sep = "")
      #url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON?inchi=", inchi, sep = "")
      url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/cids/JSON?inchikey=", inchiKey, sep = "")
      #tmp <- fromJSON("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/AADVZSXPNRLYLV/JSON")[[1]]
      
      result = tryCatch({
        fromJSON(url)[[1]]
      }, error = function(e) {
        NA
      })
      if(all(!is.na(result), result==0)) result <- NA
      return(result)
    })
  }
  
  cat(".")
  if(length(smiles) > 0){
    inchiToPubchemId_pubchem3 <- sapply(X = smiles, FUN = function(smile){
      
      ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON?inchi=InChI=1S/C3H8/c1-3-2/h3H2,1-2H3
      ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/AADVZSXPNRLYLV/JSON
      #url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/", smile, "/JSON", sep = "")
      #url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON?inchi=", inchi, sep = "")
      url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/JSON?smiles=", smile, sep = "")
      #tmp <- fromJSON("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/AADVZSXPNRLYLV/JSON")[[1]]
      
      result = tryCatch({
        fromJSON(url)[[1]]
      }, error = function(e) {
        NA
      })
      if(all(!is.na(result), result==0)) result <- NA
      return(result)
    })
  }
  
  cat(".")
  if(length(allInchis) > 0){
    #inchiToPubchemId_webchem <- sapply(X = allInchis, FUN = function(inchi){webchem::cs_inchi_csid(inchi = inchi, verbose = FALSE)})
    #inchiToPubchemId_webchem[inchiToPubchemId_webchem==""] <- NA
    
    inchiToPubchemId_pubchem <- sapply(X = allInchis, FUN = function(inchi){
      
      ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON?inchi=InChI=1S/C3H8/c1-3-2/h3H2,1-2H3
      ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/AADVZSXPNRLYLV/JSON
      #url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/", inchi, "/JSON", sep = "")
      url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON?inchi=", inchi, sep = "")
      #tmp <- fromJSON("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/AADVZSXPNRLYLV/JSON")[[1]]
      
      result = tryCatch({
        fromJSON(url)[[1]]
      }, error = function(e) {
        NA
      })
      if(all(!is.na(result), result==0)) result <- NA
      return(result)
    })
  }
  
  cat("post...")
  ## from inchikey
  pubChemResultMatrix <- matrix(data = NA, nrow = nrow(metaDataCompoundsDf), ncol = 4)
  #pubChemResultMatrix[match(x = inchiKeys, table = metaDataCompoundsDf$InChI ), 1] <- unlist(unname(inchiToPubchemId_webchem2))
  pubChemResultMatrix[match(x = inchiKeys, table = potentialInchiKeysFromTable), 1] <- unlist(unname(inchiToPubchemId_pubchem2))
  
  #isNAtmp          <- is.na(inchiToPubchemId) | inchiToPubchemId == ""
  #inchiToPubchemId <- inchiToPubchemId[!isNAtmp]
  #allInchis        <- allInchis[!isNAtmp]
  
  #metaDataCompoundsDf$`PubChem CID`[match(x = allInchis, table = metaDataCompoundsDf$InChI )] <- inchiToPubchemId
  
  ## remove remaining inchiKeys
  isInchiKey <- grepl(x = metaDataCompoundsDf$InChI, pattern = inchiKeyRegExPattern)
  metaDataCompoundsDf$InChI[isInchiKey] <- NA
  
  
  ## from smiles
  pubChemResultMatrix[match(x = smiles,    table = metaDataCompoundsDf$SMILES), 2] <- unlist(unname(inchiToPubchemId_pubchem3))
  ## from inchi
  #pubChemResultMatrix[match(x = allInchis, table = metaDataCompoundsDf$InChI ), 4] <- unlist(unname(inchiToPubchemId_webchem))
  pubChemResultMatrix[match(x = allInchis, table = metaDataCompoundsDf$InChI ), 3] <- unlist(unname(inchiToPubchemId_pubchem))
  
  pubChemResultMatrix[, 4] <- pubchemIDs
  
  ## wrap
  #resultPubchemIDs <- apply(X = pubChemResultMatrix, MARGIN = 1, FUN = function(x){paste(unique(x[!is.na(x)]), collapse = "; ")})
  resultPubchemIDsList <- apply(X = pubChemResultMatrix, MARGIN = 1, FUN = function(x){
    cids <- unique(x[!is.na(x)])
    if(length(cids) == 0) return(NA)
    else return(cids)
  })
  resultPubchemIDs <- vector(mode = "list", length = length(resultPubchemIDsList))
  resultPubchemIDs2 <- vector(mode = "list", length = length(resultPubchemIDsList))
  #for(rowIdx in 3149:length(resultPubchemIDsList)){
  for(rowIdx in seq_along(resultPubchemIDsList)){
    #print(rowIdx)
    if(any(is.na(metaDataCompoundsDf$InChI[[rowIdx]]), metaDataCompoundsDf$InChI[[rowIdx]] == "", metaDataCompoundsDf$InChI[[rowIdx]] == "NA", is.na(resultPubchemIDsList[[rowIdx]]))){
      ## no inchi to compare
      resultPubchemIDs [[rowIdx]] <- resultPubchemIDsList[[rowIdx]]
    } else {
      ## inchi to compare
      inchisHere <- unlist(sapply(X = resultPubchemIDsList[[rowIdx]], FUN = function(cid){
        
        ## Too many requests or server too busy. . Returning NA.
        inchiHere <- NA
        while(is.na(inchiHere)){
          inchiHere <- tryCatch(expr = { webchem::pc_prop(cid = cid, verbose = F)[["InChI"]] },
                                error = function(e){ if(grepl(x = e$message, pattern = "Too many requests or server too busy.")){print(e$message); Sys.sleep(time = 10); NA} else {stop(e$message)} }
          )
        }
        return(inchiHere)
        #webchem::pc_prop(cid = cid, verbose = F)[["InChI"]]
      }))
      rightPubChemId <- inchisHere == metaDataCompoundsDf$InChI[[rowIdx]]
      
      resultPubchemIDs [[rowIdx]] <- resultPubchemIDsList[[rowIdx]][ rightPubChemId]
      resultPubchemIDs2[[rowIdx]] <- resultPubchemIDsList[[rowIdx]][!rightPubChemId]
    }
    
    if(length(resultPubchemIDs [[rowIdx]]) == 0) resultPubchemIDs [[rowIdx]] <- NA
    if(length(resultPubchemIDs2[[rowIdx]]) == 0) resultPubchemIDs2[[rowIdx]] <- NA
  }
  
  resultPubchemIDs  <- unlist(lapply(X = resultPubchemIDs,  FUN = function(x){paste(x, collapse = "; ")}))
  resultPubchemIDs2 <- unlist(lapply(X = resultPubchemIDs2, FUN = function(x){paste(x, collapse = "; ")}))
  #if(is.null(resultPubchemIDs)) resultPubchemIDs <- 
  
  metaDataCompoundsDf$"PubChem CID" <- resultPubchemIDs
  metaDataCompoundsDf$"PubChem CID (secondary)" <- resultPubchemIDs2
  cat("ready")
  
  return(metaDataCompoundsDf)
}
molecularFormulaToSMILES <- function(formula){
  list   <- RMassBank:::formulastring.to.list(formula)
  atoms  <- unlist(sapply(X = seq_along(list), FUN = function(idx2){rep(x = names(list)[[idx2]], times = list[[idx2]])}))
  smiles <- paste("[", atoms, "]", sep = "", collapse = "")
  return(smiles)
}
preprocessData <- function(folder, fileMetaData, fileMetaDataCompoundInformationProcessed, ms_type){
  folderRawData    <- paste(folder, "raw data", sep = "/")
  
  #install.packages('readxl')
  library('readxl')
  
  presentSheets <- readxl::excel_sheets(path = fileMetaData)
  conventionalSheets <- c("Compounds", "Chromatography", "Mass_Spectrometry")
  if(!all(conventionalSheets %in% presentSheets))  stop(paste("Sheets", paste(conventionalSheets[!(conventionalSheets %in% presentSheets)], collapse = "; "), "are not present"))
  
  metaDataCompoundsDf        <- as.data.frame(read_excel(path = fileMetaData, sheet = "Compounds"))
  metaDataChromatographyDf   <- as.data.frame(read_excel(path = fileMetaData, sheet = "Chromatography"))
  metaDataMassSpectrometryDf <- as.data.frame(read_excel(path = fileMetaData, sheet = "Mass_Spectrometry"))
  
  ## sanity checks TODO
  #if(!(ncol(metaDataCompoundsDf)        %in% c(23, 24))) stop("Sheet 'Compounds' is not conventional")
  #if(ncol(metaDataChromatographyDf)   !=  2) stop("Sheet 'Chromatography' is not conventional")
  #if(ncol(metaDataMassSpectrometryDf) !=  2) stop("Sheet 'Mass_Spectrometry' is not conventional")
  
  ## reformat and check whether the parameters are default and can therrefore not be used
  colnames(metaDataCompoundsDf) <- metaDataCompoundsDf[1, ]
  colnames(metaDataCompoundsDf)[[which(is.na(colnames(metaDataCompoundsDf)) | colnames(metaDataCompoundsDf)=="NA")]] <- "Sample"
  metaDataCompoundsDf <- metaDataCompoundsDf[-1, ]
  
  ## delete borderung white spaces
  for(colIdx in seq_len(ncol(metaDataCompoundsDf)))
    if(is.character(metaDataCompoundsDf[, colIdx]))
      metaDataCompoundsDf[, colIdx] <- trimws(metaDataCompoundsDf[, colIdx])
  
  ## TODO warn if only a subset of the information is updated?
  chromatographyPropertyNames   <- c("COLUMN_NAME", "FLOW_GRADIENT", "FLOW_RATE", "RETENTION_TIME", "SOLVENT")
  chromatographyValueNames      <- c("Symmetry C18 Column, Waters", "0min:5%, 24min:95%, 28min:95%, 28.1:5% (acetonitrile)", "0.3 ml/min", "CH3CN(0.1%HCOOH)/ H2O(0.1%HCOOH)")
  massSpectrometryPropertyNames <- c("MS_TYPE", "ION_MODE", "ACTIVATION_PARAMETER", "ACTIVATION_TIME", "AUTOMATIC_GAIN_CONTROL", "CAPILLARY_TEMPERATURE", "CAPILLARY_VOLTAGE", "COLLISION_ENERGY", "FRAGMENTATION_METHOD", "IONIZATION", "RESOLUTION_SETTING", "SPRAY_VOLTAGE", "TUBE_LENS_VOLTAGE")
  massSpectrometryValueNames    <- c("MS2", "POSITIVE", "q=0.25", "30 ms", "3.0E5", "275 C", "39 V", "35eV", "CID", "ESI", "7500", "4.5 kV", "140 V")
  chromatographyInformationCanByUsed   <- !( all(metaDataChromatographyDf$Property   %in% chromatographyPropertyNames  ) & all(metaDataChromatographyDf$Value[-4]  %in% chromatographyValueNames  ) )
  massSpectrometryInformationCanByUsed <- !( all(metaDataMassSpectrometryDf$Property %in% massSpectrometryPropertyNames) & all(metaDataMassSpectrometryDf$Value    %in% massSpectrometryValueNames) )
  
  ## create ID column
  metaDataCompoundsDf$ID <- 1:nrow(metaDataCompoundsDf)
  
  ##################################################################################
  ## validate meta data sheet content
  
  ## "File"                   "Name"                   "Synonyms"               "InChI"                  "SMILES"                
  ## "PubChem CID"            "Formula"                "Exact mass"             "CHROMATOGRAPHY"         "RT (min)"              
  ## "INSTRUMENT"             "INSTRUMENT_TYPE"        "IONIZATION"             "Collision energy"       "Adduct"                
  ## "MS/MS Acquisition mode" "Ionization mode"        "Confidence"             "Compound class"         "Sample"                
  ## "Database links"         "Authors"                "Publication"           
  
  ## validate column names
  columnNames  <- c("ID", "File", "Name", "Synonyms", "InChI", "SMILES", "PubChem CID", "Formula", "Exact mass", "CHROMATOGRAPHY", "RT (min)",
                    "INSTRUMENT", "INSTRUMENT_TYPE", "IONIZATION", "Collision energy", "Adduct", "MS/MS Acquisition mode", "Ionization mode",
                    "Confidence", "Compound class", "Sample", "Database links", "Authors", "Publication")
  #columnNames2 <- c(columnNames, "INSTRUMENT_SETTING_ID", "CHROM_SETTING_ID")
  
  if(!all(columnNames %in% colnames(metaDataCompoundsDf))) stop(paste("Column names of sheet 'Compounds' are not conventional - missingcolumns: ", paste(columnNames[!(columnNames %in% colnames(metaDataCompoundsDf))], collapse="; ")))
  column_INSTRUMENT_SETTING_ID_isThere     <- "INSTRUMENT_SETTING_ID"     %in% colnames(metaDataCompoundsDf)
  column_CHROM_SETTING_ID_isThere <- "CHROM_SETTING_ID" %in% colnames(metaDataCompoundsDf)
  
  ## validate mandatory columns
  metaDataCompoundsDfEmpty <- data.frame(is.na(metaDataCompoundsDf) | (metaDataCompoundsDf == "") | (metaDataCompoundsDf == "NA") | (metaDataCompoundsDf == "N/A"), check.names = FALSE)
  invalidIndeces <- vector(mode = "integer", length = 0)
  for(rowIdx in seq_len(nrow(metaDataCompoundsDf))){
    if(ms_type=="MS"){
      if(
        metaDataCompoundsDfEmpty[rowIdx, "File"] | 
        metaDataCompoundsDfEmpty[rowIdx, "Name"] |
        (
          metaDataCompoundsDfEmpty[rowIdx, "InChI"] & 
          metaDataCompoundsDfEmpty[rowIdx, "SMILES"] & 
          metaDataCompoundsDfEmpty[rowIdx, "PubChem CID"]
        ) | 
        metaDataCompoundsDfEmpty[rowIdx, "RT (min)"] | 
        metaDataCompoundsDfEmpty[rowIdx, "INSTRUMENT"] | 
        metaDataCompoundsDfEmpty[rowIdx, "INSTRUMENT_TYPE"] | 
        metaDataCompoundsDfEmpty[rowIdx, "IONIZATION"] | 
        metaDataCompoundsDfEmpty[rowIdx, "Collision energy"] | 
        #metaDataCompoundsDfEmpty[rowIdx, "Adduct"] | 
        metaDataCompoundsDfEmpty[rowIdx, "MS/MS Acquisition mode"] | 
        metaDataCompoundsDfEmpty[rowIdx, "Ionization mode"] | 
        metaDataCompoundsDfEmpty[rowIdx, "Confidence"] | 
        metaDataCompoundsDfEmpty[rowIdx, "Authors"]
      ){
        invalidIndeces[[length(invalidIndeces) + 1]] <- rowIdx
        cat(paste("\n### Warning ### Compound ID ", metaDataCompoundsDf[rowIdx, "ID"], " is missing mandatory information:", paste(
          c("File",
            "Name",
            "Structure",
            "RT (min)",
            "INSTRUMENT",
            "INSTRUMENT_TYPE",
            "IONIZATION",
            "Collision energy",
            #"Adduct",
            "MS/MS Acquisition mode",
            "Ionization mode",
            "Confidence",
            "Authors")[
              c(
                metaDataCompoundsDfEmpty[rowIdx, "File"],
                metaDataCompoundsDfEmpty[rowIdx, "Name"],
                (
                  metaDataCompoundsDfEmpty[rowIdx, "InChI"] & 
                    metaDataCompoundsDfEmpty[rowIdx, "SMILES"] & 
                    metaDataCompoundsDfEmpty[rowIdx, "PubChem CID"]
                ),
                metaDataCompoundsDfEmpty[rowIdx, "RT (min)"],
                metaDataCompoundsDfEmpty[rowIdx, "INSTRUMENT"],
                metaDataCompoundsDfEmpty[rowIdx, "INSTRUMENT_TYPE"],
                metaDataCompoundsDfEmpty[rowIdx, "IONIZATION"],
                metaDataCompoundsDfEmpty[rowIdx, "Collision energy"],
                #metaDataCompoundsDfEmpty[rowIdx, "Adduct"],
                metaDataCompoundsDfEmpty[rowIdx, "MS/MS Acquisition mode"],
                metaDataCompoundsDfEmpty[rowIdx, "Ionization mode"],
                metaDataCompoundsDfEmpty[rowIdx, "Confidence"],
                metaDataCompoundsDfEmpty[rowIdx, "Authors"]
              )
              ], collapse = "; "
        ), sep = ""))
      }
    } else {
      if(
        metaDataCompoundsDfEmpty[rowIdx, "File"] | 
        metaDataCompoundsDfEmpty[rowIdx, "Name"] |
        (
          metaDataCompoundsDfEmpty[rowIdx, "InChI"] & 
          metaDataCompoundsDfEmpty[rowIdx, "SMILES"] & 
          metaDataCompoundsDfEmpty[rowIdx, "PubChem CID"]
        ) | 
        metaDataCompoundsDfEmpty[rowIdx, "RT (min)"] | 
        metaDataCompoundsDfEmpty[rowIdx, "INSTRUMENT"] | 
        metaDataCompoundsDfEmpty[rowIdx, "INSTRUMENT_TYPE"] | 
        metaDataCompoundsDfEmpty[rowIdx, "IONIZATION"] | 
        metaDataCompoundsDfEmpty[rowIdx, "Collision energy"] | 
        metaDataCompoundsDfEmpty[rowIdx, "Adduct"] | 
        metaDataCompoundsDfEmpty[rowIdx, "MS/MS Acquisition mode"] | 
        metaDataCompoundsDfEmpty[rowIdx, "Ionization mode"] | 
        metaDataCompoundsDfEmpty[rowIdx, "Confidence"] | 
        metaDataCompoundsDfEmpty[rowIdx, "Authors"]
      ){
        invalidIndeces[[length(invalidIndeces) + 1]] <- rowIdx
        cat(paste("\n### Warning ### Compound ID ", metaDataCompoundsDf[rowIdx, "ID"], " is missing mandatory information:", paste(
          c("File",
            "Name",
            "Structure",
            "RT (min)",
            "INSTRUMENT",
            "INSTRUMENT_TYPE",
            "IONIZATION",
            "Collision energy",
            "Adduct",
            "MS/MS Acquisition mode",
            "Ionization mode",
            "Confidence",
            "Authors")[
              c(
                metaDataCompoundsDfEmpty[rowIdx, "File"],
                metaDataCompoundsDfEmpty[rowIdx, "Name"],
                (
                  metaDataCompoundsDfEmpty[rowIdx, "InChI"] & 
                    metaDataCompoundsDfEmpty[rowIdx, "SMILES"] & 
                    metaDataCompoundsDfEmpty[rowIdx, "PubChem CID"]
                ),
                metaDataCompoundsDfEmpty[rowIdx, "RT (min)"],
                metaDataCompoundsDfEmpty[rowIdx, "INSTRUMENT"],
                metaDataCompoundsDfEmpty[rowIdx, "INSTRUMENT_TYPE"],
                metaDataCompoundsDfEmpty[rowIdx, "IONIZATION"],
                metaDataCompoundsDfEmpty[rowIdx, "Collision energy"],
                metaDataCompoundsDfEmpty[rowIdx, "Adduct"],
                metaDataCompoundsDfEmpty[rowIdx, "MS/MS Acquisition mode"],
                metaDataCompoundsDfEmpty[rowIdx, "Ionization mode"],
                metaDataCompoundsDfEmpty[rowIdx, "Confidence"],
                metaDataCompoundsDfEmpty[rowIdx, "Authors"]
              )
              ], collapse = "; "
        ), sep = ""))
      }
    }
  }
  
  if(length(invalidIndeces) > 0){
    ## TODO add to global log
    cat(paste("\n### Warning ### Rows ", paste(invalidIndeces, collapse = "; "), " of sheet 'Compounds' are missing mandatory information:", sep = ""))
    #cat(paste("\n", paste(apply(X = metaDataCompoundsDf[invalidIndeces, c("File", "Name")], MARGIN = 1, FUN = function(row){paste(row, collapse = "\t")}), collapse = "\n"), sep = ""))
    stop("")
    metaDataCompoundsDf <- metaDataCompoundsDf[-invalidIndeces, ]
  }
  
  if(chromatographyInformationCanByUsed){
    #chromatographyInformationCanByUsed <- FALSE
    ## do not validate chromatography input
    #stop("Not implemented yet")
    if(ncol(metaDataChromatographyDf) > 2 && !column_CHROM_SETTING_ID_isThere) stop("More than one chrom setting specified without column CHROM_SETTING_ID")
    if((column_CHROM_SETTING_ID_isThere & !("CHROM_SETTING_ID" %in% metaDataChromatographyDf$Property)) | (!column_CHROM_SETTING_ID_isThere & ("CHROM_SETTING_ID" %in% metaDataChromatographyDf$Property))) stop("Columns CHROM_SETTING_ID not consistent")
    #if(!("CHROM_SETTING_ID" %in% metaDataChromatographyDf$Property)) stop("No CHROM_SETTING_ID in Property column of sheet Mass_Spectrometry")
    if(sum(duplicated(metaDataChromatographyDf[, 1])) > 0) stop("Duplicated chromatography meta data properties")
    
    rownames(metaDataChromatographyDf) <- metaDataChromatographyDf[, 1]
    metaDataChromatographyDf <- metaDataChromatographyDf[, -1, drop=FALSE]
    
    ## exclude RT as it makes no sence to specify it globally
    if(any(rownames(metaDataChromatographyDf) %in% c("RETENTION_TIME")))
      metaDataChromatographyDf <- metaDataChromatographyDf[-which(rownames(metaDataChromatographyDf) %in% c("RETENTION_TIME")), , drop=FALSE]
    
    ## do not check chromatography parameters
    ## - ac_lc[['CAPILLARY_VOLTAGE']]  <- 
    ## + ac_lc[['COLUMN_NAME']]        <- getOption("RMassBank")$annotations$lc_column
    ## - ac_lc[['COLUMN_TEMPERATURE']] <- 
    ## + ac_lc[['FLOW_GRADIENT']]      <- getOption("RMassBank")$annotations$lc_gradient
    ## + ac_lc[['FLOW_RATE']]          <- getOption("RMassBank")$annotations$lc_flow
    ## + ac_lc[['RETENTION_TIME']]     <- 
    ## + ac_lc[['SOLVENT A']]          <- getOption("RMassBank")$annotations$lc_solvent_a
    ## + ac_lc[['SOLVENT B']]          <- getOption("RMassBank")$annotations$lc_solvent_b
    
    #if(any(metaDataChromatographyDf[, 1] == "SOLVENT")) metaDataChromatographyDf[metaDataChromatographyDf[, 1] == "SOLVENT", 1] <- "SOLVENT A"
    #
    #chromProperties1 <- c("lc_column", "lc_gradient", "lc_flow", "lc_solvent_a", "lc_solvent_b")
    #chromProperties2 <- c("COLUMN_NAME", "FLOW_GRADIENT", "FLOW_RATE", "SOLVENT A", "SOLVENT B")
    #if(!all(metaDataChromatographyDf[, 1] %in% chromProperties2)) stop(paste("Chromatography property not conventional:", paste(metaDataChromatographyDf[, 1][!(metaDataChromatographyDf[, 1] %in% chromProperties2)], collapse = "; ")))
    
    row_CHROM_SETTING_ID <- ifelse(test = column_CHROM_SETTING_ID_isThere, yes = which(rownames(metaDataChromatographyDf) == "CHROM_SETTING_ID"), no = NA)
  } else row_CHROM_SETTING_ID <- NA
  if(massSpectrometryInformationCanByUsed){
    #massSpectrometryInformationCanByUsed <- FALSE
    ## do not validate MS input
    #stop("Not implemented yet")
    if(ncol(metaDataMassSpectrometryDf) > 2 && !column_INSTRUMENT_SETTING_ID_isThere) stop("More than one MS setting specified without column INSTRUMENT_SETTING_ID")
    if((column_INSTRUMENT_SETTING_ID_isThere & !("INSTRUMENT_SETTING_ID" %in% metaDataMassSpectrometryDf$Property)) | (!column_INSTRUMENT_SETTING_ID_isThere & ("INSTRUMENT_SETTING_ID" %in% metaDataMassSpectrometryDf$Property))) stop("Columns INSTRUMENT_SETTING_ID not consistent")
    #if(!("INSTRUMENT_SETTING_ID" %in% metaDataMassSpectrometryDf$Property)) stop("No INSTRUMENT_SETTING_ID in Property column of sheet Mass_Spectrometry")
    if(sum(duplicated(metaDataMassSpectrometryDf[, 1])) > 0) stop("Duplicated MS meta data properties")
    
    rownames(metaDataMassSpectrometryDf) <- metaDataMassSpectrometryDf[, 1]
    metaDataMassSpectrometryDf <- metaDataMassSpectrometryDf[, -1, drop=FALSE]
    
    ## exclude ION_MODE,  as it is already specified in the sheet 'Compounds'
    if(any(rownames(metaDataChromatographyDf) %in% c("COLLISION_ENERGY", "ION_MODE",  "IONIZATION")))
      metaDataChromatographyDf <- metaDataChromatographyDf[-which(rownames(metaDataChromatographyDf) %in% c("COLLISION_ENERGY", "ION_MODE",  "IONIZATION")), , drop=FALSE]
    
    ## do not check mass spectrometry parameters
    ## -         ACTIVATION_PARAMETER
    ## -         ACTIVATION_TIME
    ## -         AUTOMATIC_GAIN_CONTROL
    ## -         CAPILLARY_TEMPERATURE
    ## -         CAPILLARY_VOLTAGE
    ## + ac_ms[['COLLISION_ENERGY']]        <- substring(grep('AC$MASS_SPECTROMETRY: COLLISION_ENERGY',record, value = TRUE, fixed = TRUE),40)
    ## - ac_ms[['COLLISION_GAS']]           <- substring(grep('AC$MASS_SPECTROMETRY: COLLISION_GAS',record, value = TRUE, fixed = TRUE),37)
    ## - ac_ms[['DATE']]                    <- substring(grep('AC$MASS_SPECTROMETRY: DATE',record, value = TRUE, fixed = TRUE),28)
    ## - ac_ms[['DESOLVATION_GAS_FLOW']]    <- substring(grep('AC$MASS_SPECTROMETRY: DESOLVATION_GAS_FLOW',record, value = TRUE, fixed = TRUE),44)
    ## - ac_ms[['DESOLVATION_TEMPERATURE']] <- substring(grep('AC$MASS_SPECTROMETRY: DESOLVATION_TEMPERATURE',record, value = TRUE, fixed = TRUE),47)
    ## -         FRAGMENTATION_METHOD
    ## + ac_ms[['FRAGMENTATION_MODE']]      <- msmsdata@info$mode
    ## + ac_ms[['ION_MODE']]                <- mode
    ## + ac_ms[['IONIZATION']]              <- getOption("RMassBank")$annotations$ionization
    ## - ac_ms[['IONIZATION_ENERGY']]       <- substring(grep('AC$MASS_SPECTROMETRY: IONIZATION_ENERGY',record, value = TRUE, fixed = TRUE),41)
    ## - ac_ms[['LASER']]                   <- substring(grep('AC$MASS_SPECTROMETRY: LASER',record, value = TRUE, fixed = TRUE),29)
    ## - ac_ms[['MATRIX']]                  <- substring(grep('AC$MASS_SPECTROMETRY: MATRIX',record, value = TRUE, fixed = TRUE),30)
    ## - ac_ms[['MASS_ACCURACY']]           <- substring(grep('AC$MASS_SPECTROMETRY: MASS_ACCURACY',record, value = TRUE, fixed = TRUE),37)
    ## + ac_ms[['MS_TYPE']]                 <- getOption("RMassBank")$annotations$ms_type
    ## - ac_ms[['REAGENT_GAS']]             <- substring(grep('AC$MASS_SPECTROMETRY: REAGENT_GAS',record, value = TRUE, fixed = TRUE),35)
    ## + ac_ms[['RESOLUTION']]              <- msmsdata@info$res
    ## -         RESOLUTION_SETTING
    ## - ac_ms[['SCANNING']]                <- substring(grep('AC$MASS_SPECTROMETRY: SCANNING',record, value = TRUE, fixed = TRUE),32)
    ## -         SPRAY_VOLTAGE
    ## -         TUBE_LENS_VOLTAGE
    
    row_INSTRUMENT_SETTING_ID <- ifelse(test = column_INSTRUMENT_SETTING_ID_isThere, yes = which(rownames(metaDataMassSpectrometryDf) == "INSTRUMENT_SETTING_ID"), no = NA)
  } else row_INSTRUMENT_SETTING_ID <- NA
  
  cat("\nMeta data sheets are valid")
  if(chromatographyInformationCanByUsed)   cat("\nChromatography meta data is supplied") else cat("\nNo chromatography meta data supplied")
  if(massSpectrometryInformationCanByUsed) cat("\nMS meta data is supplied")             else cat("\nNo MS meta data supplied")
  
  #################################################################################################################
  ## validate files
  
  metaDataCompoundsDf$MspFile <- gsub(x = metaDataCompoundsDf$File, pattern = "\\.((abf)|(ABF)|(raw)|(RAW)|(mgf)|(MGF))$", replacement = ".msp")
  if(!any(grepl(metaDataCompoundsDf$MspFile, pattern = "\\.msp$"))) metaDataCompoundsDf$MspFile <- paste(metaDataCompoundsDf$MspFile, ".msp", sep = "")
  
  fileMspData <- list.files(path = folderRawData, pattern = "\\.msp$", full.names = TRUE, recursive = TRUE)
  if(any(grepl(x = metaDataCompoundsDf$MspFile, pattern = "\\\\"))) metaDataCompoundsDf$MspFile <- unlist(lapply(X = strsplit(x = metaDataCompoundsDf$MspFile, split = "\\\\"), FUN = tail, 1))
  metaDataCompoundsDf$MspFileThere <- metaDataCompoundsDf$MspFile %in% basename(fileMspData)
  if(!all(metaDataCompoundsDf$MspFileThere)){
    if(all(gsub(x = toupper(metaDataCompoundsDf$MspFile), pattern = "[ _-]", replacement = "") %in% gsub(x = toupper(basename(fileMspData)), pattern = "[ _-]", replacement = ""))){
      cat(paste("\nCorrecting", sum(!metaDataCompoundsDf$MspFileThere), "file names"))
      matchIndeces <- match(x = gsub(x = toupper(metaDataCompoundsDf$MspFile), pattern = "[ _-]", replacement = ""), table = gsub(x = toupper(basename(fileMspData)), pattern = "[ _-]", replacement = ""))
      cat(paste("\n    ", metaDataCompoundsDf$MspFile, " > ", basename(fileMspData)[matchIndeces], sep = ""))
      metaDataCompoundsDf$MspFile <- basename(fileMspData)[matchIndeces]
      metaDataCompoundsDf$MspFileThere <- metaDataCompoundsDf$MspFile %in% basename(fileMspData)
    }
    if(
      all(gsub(x = toupper(metaDataCompoundsDf$MspFile), pattern = "((_POS)|(_NEG))\\.MSP$", replacement = ".MSP") %in% toupper(basename(fileMspData))) & 
      all(grepl(x = fileMspData, pattern = "/((pos)|(neg))/"))
    ){
      filesPos <- fileMspData[grepl(x = fileMspData, pattern = "/pos/")]
      filesNeg <- fileMspData[grepl(x = fileMspData, pattern = "/neg/")]
      if(!all(file.rename(from = filesPos, to = gsub(x = filesPos, pattern = "\\.msp", replacement = "_pos.msp")))) stop(paste("Error renaming files", paste(filesPos, collapse="; ")))
      if(!all(file.rename(from = filesNeg, to = gsub(x = filesNeg, pattern = "\\.msp", replacement = "_neg.msp")))) stop(paste("Error renaming files", paste(filesNeg, collapse="; ")))
      fileMspData <- list.files(path = folderRawData, pattern = "\\.msp$", full.names = TRUE, recursive = TRUE)
      metaDataCompoundsDf$MspFileThere <- metaDataCompoundsDf$MspFile %in% basename(fileMspData)
    }
  }
  
  if(!all(metaDataCompoundsDf$MspFileThere)){
    stop(paste("\n### Error ###", length(unique(metaDataCompoundsDf$MspFile[!metaDataCompoundsDf$MspFileThere])), "Msp files", paste(unique(metaDataCompoundsDf$MspFile[!metaDataCompoundsDf$MspFileThere]), collapse = "; "), "are missing"))
  }
  metaDataCompoundsDf$"AbfFile" <- metaDataCompoundsDf$"File"
  metaDataCompoundsDf$"File" <- sapply(X = metaDataCompoundsDf$MspFile, FUN = function(file){
    if(!(file %in% basename(fileMspData)))  return(file)
    if(sum(basename(fileMspData) %in% file) > 1) stop(paste("File", file, "is ambiguous:", fileMspData[basename(fileMspData) %in% file]))
    else return(fileMspData[[match(x = file, table = basename(fileMspData))]])
  })
  
  ######################################################################################################
  ## TODO validate remaining fields
  library("stringr")
  numberOfOpeningBrackets <- str_count(string = metaDataCompoundsDf$SMILES, pattern = "\\(")
  numberOfClosingBrackets <- str_count(string = metaDataCompoundsDf$SMILES, pattern = "\\)")
  if(any(numberOfOpeningBrackets != numberOfClosingBrackets, na.rm = TRUE)) stop(paste("SMILES of compounds", paste(metaDataCompoundsDf$ID[numberOfOpeningBrackets != numberOfClosingBrackets], collapse = "; "), "are invalid:", paste(metaDataCompoundsDf$SMILES[numberOfOpeningBrackets != numberOfClosingBrackets], collapse = "; ")))
  
  noStructure <- apply(X = metaDataCompoundsDf, MARGIN = 1, FUN = function(row){all(row[c("InChI", "SMILES", "PubChem CID")] == "[none]")})
  if(any(noStructure))
    for(idx in which(noStructure))
      metaDataCompoundsDf[idx, c("InChI", "SMILES", "PubChem CID")] <- c("", molecularFormulaToSMILES(metaDataCompoundsDf$Formula[[idx]]), "")
  
  invalidInchis <- NULL
  inchiSet <- metaDataCompoundsDf$InChI[!(is.na(metaDataCompoundsDf$InChI) | metaDataCompoundsDf$InChI == "")]
  inchiSet <- inchiSet[!grepl(x = inchiSet, pattern = "^[A-Z]{14,14}-[A-Z]{10,10}-[A-Z]$")]
  inchiSet <- inchiSet[inchiSet != "NA"]
  for(inchi in inchiSet)
    tryCatch(expr = { rinchi::parse.inchi(inchis = inchi)[[1]] },
             error = function(e){ invalidInchis <<- c(invalidInchis, inchi) }    
    )
  
  invalidSMILES <- NULL
  smilesSet <- metaDataCompoundsDf$SMILES[!(is.na(metaDataCompoundsDf$SMILES) | metaDataCompoundsDf$SMILES == "")]
  smilesSet <- smilesSet[smilesSet != "NA"]
  for(smiles in smilesSet)
    tryCatch(expr = { rcdk::parse.smiles(smiles = smiles)[[1]] },
             error = function(e){ invalidSMILES <<- c(invalidSMILES, smiles) }    
    )
  
  if(length(invalidInchis) > 0)
    stop(paste("invalid inchis:\n", paste(invalidInchis, collapse = "\n"), sep = ""))
  if(length(invalidSMILES) > 0)
    stop(paste("invalid SMILES:\n", paste(invalidSMILES, collapse = "\n"), sep = ""))
  
  metaDataCompoundsDf$"Adduct"[metaDataCompoundsDf$"Adduct"=="[M+FA-H]-"]  <- "[M+HCOOH-H]-"
  metaDataCompoundsDf$"Adduct"[metaDataCompoundsDf$"Adduct"=="[2M+FA-H]-"] <- "[2M+HCOOH-H]-"
  metaDataCompoundsDf$"Adduct"[metaDataCompoundsDf$"Adduct"=="[M+TFA-H]-"] <- "[M+CF3CO2H-H]-"
  metaDataCompoundsDf$"Adduct"[metaDataCompoundsDf$"Adduct"=="[M+H-H2O]+"] <- "[M-H2O+H]+"
  
  if(!all(metaDataCompoundsDf$"Ionization mode" %in% c("Positive", "Negative"))) stop(paste("Invalid values in column 'Ionization mode'", paste(unique(metaDataCompoundsDf$"Ionization mode"[!(metaDataCompoundsDf$"Ionization mode" %in% c("Positive", "Negative"))]), collapse = "; ")))
  if(ms_type == "MS"){
    if(!all(is.na(metaDataCompoundsDf$"Adduct"))) stop(paste("Invalid values in column 'Adduct'", paste(unique(metaDataCompoundsDf$"Adduct"[!is.na(metaDataCompoundsDf$"Adduct")]), collapse = "; ")))
  } else
    if(!all(metaDataCompoundsDf$"Adduct" %in% knownAdducts)) stop(paste("Invalid values in column 'Adduct'", paste(unique(metaDataCompoundsDf$"Adduct"[!(metaDataCompoundsDf$"Adduct" %in% knownAdducts)]), collapse = "; ")))
  
  if(any(grepl(x = metaDataCompoundsDf$"RT (min)", pattern = ","))) metaDataCompoundsDf$"RT (min)" <- gsub(x = metaDataCompoundsDf$"RT (min)", pattern = ",", replacement = ".")
  metaDataCompoundsDf$"RT (min)"[metaDataCompoundsDf$"RT (min)"==""] <- NA
  metaDataCompoundsDf$"RT (min)"[metaDataCompoundsDf$"RT (min)"=="NA"] <- NA
  if(!all(grepl(x = metaDataCompoundsDf$"RT (min)"[!is.na(metaDataCompoundsDf$"RT (min)")], pattern = "^\\d+(\\.(\\d+)?)?$"))) stop("Invalid column 'RT (min)'")
  if(any(grepl(x = metaDataCompoundsDf$"Exact mass", pattern = ","))) metaDataCompoundsDf$"Exact mass" <- gsub(x = metaDataCompoundsDf$"Exact mass", pattern = ",", replacement = ".")
  #if(!all(grepl(x = metaDataCompoundsDf$"Exact mass", pattern = "^\\d+(\\.(\\d+)?)?$"))) stop("Invalid column 'Exact mass'")
  
  metaDataCompoundsDf$"Synonyms" <- gsub(x = metaDataCompoundsDf$"Synonyms", pattern = "\"", replacement = "")
  metaDataCompoundsDf$"Synonyms" <- gsub(x = metaDataCompoundsDf$"Synonyms", pattern = "\\\\", replacement = "")
  
  #metaDataCompoundsDf$"Synonyms"[grepl(x = metaDataCompoundsDf$"Synonyms", pattern = "<U\\+[A-Za-z0-9]{4,4}>")]
  library("stringr")
  knownUnicodeSymbols <- c(NA, "<U+03B2>", "<U+03B1>", "<U+2192>")
  presentUnicodeSymbols <- unique(unlist(str_extract_all(string = metaDataCompoundsDf$"Synonyms", pattern = "<U\\+[A-Za-z0-9]{4,4}>")))
  if(!all(presentUnicodeSymbols %in% knownUnicodeSymbols)) stop(paste("Unrecognized unicode characters", presentUnicodeSymbols[presentUnicodeSymbols %in% knownUnicodeSymbols]))
  metaDataCompoundsDf$"Synonyms" <- gsub(x = metaDataCompoundsDf$"Synonyms", pattern = "<U\\+03B1>", replacement = "\u03b1") ## alpha
  metaDataCompoundsDf$"Synonyms" <- gsub(x = metaDataCompoundsDf$"Synonyms", pattern = "<U\\+03B2>", replacement = "\u03b2") ## beta
  metaDataCompoundsDf$"Synonyms" <- gsub(x = metaDataCompoundsDf$"Synonyms", pattern = "<U\\+2192>", replacement = "\u2192") ## rightwards arrow
  
  ##################################################################################
  ## collect Chromatography
  if(column_CHROM_SETTING_ID_isThere)
    if(!all(metaDataCompoundsDf$"CHROM_SETTING_ID" %in% metaDataChromatographyDf[row_CHROM_SETTING_ID, ])) stop("Missing CHROM_SETTING_ID", metaDataCompoundsDf$"CHROM_SETTING_ID"[!(metaDataCompoundsDf$"CHROM_SETTING_ID" %in% metaDataChromatographyDf[row_CHROM_SETTING_ID, ])])
  
  if(chromatographyInformationCanByUsed){
    ## TODO value of column CHROMATOGRAPHY vs value of property COLUMN_NAME ??
    ## TODO unify subtags: https://github.com/MassBank/MassBank-web/blob/master/Documentation/MassBankRecordFormat.md#2.4.6
    
    ## add columns to df
    chromatographyColumns <- rownames(metaDataChromatographyDf)
    if(column_CHROM_SETTING_ID_isThere) chromatographyColumns <- chromatographyColumns[chromatographyColumns != "CHROM_SETTING_ID"]
    if(any(grepl(x = colnames(metaDataCompoundsDf), pattern = "^CHROMATOGRAPHY_SETTING_"))) stop("CHROMATOGRAPHY_SETTING column name pattern already there")
    
    ## TODO AC$CHROMATOGRAPHY: COLUMN_NAME ACQUITY UPLC HSS T3 Column
    ##   vs AC$CHROMATOGRAPHY: COLUMN_NAME HSS T3 column (100 × 1 mm, particle size 1.8 µm; Waters
    ##   -> AC$CHROMATOGRAPHY: COLUMN_SPECIFICATION HSS T3 column (100 × 1 mm, particle size 1.8 µm; Waters
    
    chromatographyColumns_tagged <- paste("CHROMATOGRAPHY_SETTING", chromatographyColumns, sep = "_")
    metaDataCompoundsDf[, chromatographyColumns_tagged] <- NA
    for(colIdx in seq_len(ncol(metaDataChromatographyDf))){
      values <- metaDataChromatographyDf[, colIdx]
      values[values==""] <- NA
      if(!is.na(row_CHROM_SETTING_ID))  values <- values[-row_CHROM_SETTING_ID]
      if(column_CHROM_SETTING_ID_isThere)  rowIndeces <- which(metaDataCompoundsDf$"CHROM_SETTING_ID" == metaDataChromatographyDf[row_CHROM_SETTING_ID, colIdx])
      else                                          rowIndeces <- seq_len(nrow(metaDataCompoundsDf))
      if(length(rowIndeces) == 0) stop("No cpds for certain CHROM_SETTING_ID")
      
      for(rowIdx in rowIndeces) metaDataCompoundsDf[rowIdx, chromatographyColumns_tagged] <- values
    }
    
    metaDataCompoundsDf$"CHROM_SETTING_ID" <- NULL
  }
  
  ##################################################################################
  ## collect Mass_Spectrometry
  if(column_INSTRUMENT_SETTING_ID_isThere)
    if(!all(metaDataCompoundsDf$"INSTRUMENT_SETTING_ID"[!is.na(metaDataCompoundsDf$"INSTRUMENT_SETTING_ID")] %in% metaDataMassSpectrometryDf[row_INSTRUMENT_SETTING_ID, ])) stop("Missing INSTRUMENT_SETTING_ID", metaDataCompoundsDf$"INSTRUMENT_SETTING_ID"[!is.na(metaDataCompoundsDf$"INSTRUMENT_SETTING_ID")][!(metaDataCompoundsDf$"INSTRUMENT_SETTING_ID"[!is.na(metaDataCompoundsDf$"INSTRUMENT_SETTING_ID")] %in% metaDataMassSpectrometryDf[row_INSTRUMENT_SETTING_ID, ])])
  
  if(massSpectrometryInformationCanByUsed){
    ## add columns to df
    msColumns <- rownames(metaDataMassSpectrometryDf)
    if(column_INSTRUMENT_SETTING_ID_isThere) msColumns <- msColumns[msColumns != "INSTRUMENT_SETTING_ID"]
    if(any(grepl(x = colnames(metaDataCompoundsDf), pattern = "^MS_SETTING_"))) stop("MS_SETTING column name pattern already there")
    
    ## TODO AC$MASS_SPECTROMETRY: ACQUISITION_MODE DDA
    ## TODO unify subtags: https://github.com/MassBank/MassBank-web/blob/master/Documentation/MassBankRecordFormat.md#2.4.5
    
    msColumns <- paste("MS_SETTING", msColumns, sep = "_")
    metaDataCompoundsDf[, msColumns] <- NA
    for(colIdx in seq_len(ncol(metaDataMassSpectrometryDf))){
      values <- metaDataMassSpectrometryDf[, colIdx]
      values[values==""] <- NA
      if(!is.na(row_INSTRUMENT_SETTING_ID))  values <- values[-row_INSTRUMENT_SETTING_ID]
      if(column_INSTRUMENT_SETTING_ID_isThere)  rowIndeces <- which(metaDataCompoundsDf$"INSTRUMENT_SETTING_ID" == metaDataMassSpectrometryDf[row_INSTRUMENT_SETTING_ID, colIdx])
      else                                      rowIndeces <- seq_len(nrow(metaDataCompoundsDf))
      if(length(rowIndeces) == 0) stop("No cpds for certain INSTRUMENT_SETTING_ID")
      
      for(rowIdx in rowIndeces) metaDataCompoundsDf[rowIdx, msColumns] <- values
    }
    
    metaDataCompoundsDf$"INSTRUMENT_SETTING_ID" <- NULL
  }
  
  ## collect compound information
  metaDataCompoundsDf <- preprocessData_compoundInformation(metaDataCompoundsDf)
  
  ## export
  cat(paste("\nWriting compound-meta-data file", fileMetaDataCompoundInformationProcessed))
  write.table(x = metaDataCompoundsDf, file = fileMetaDataCompoundInformationProcessed, sep = "\t")
  #if(chromatographyInformationCanByUsed){
  #  cat(paste("\nWriting chromatography-meta-data file", fileMetaDataChromatographyInformationProcessed))
  #  write.table(x = metaDataChromatographyDf, file = fileMetaDataChromatographyInformationProcessed, sep = "\t")
  #}
  #if(massSpectrometryInformationCanByUsed){
  #  cat(paste("\nWriting mass-spectrometry-meta-data file", fileMetaDataMassSpectrometryInformationProcessed))
  #  write.table(x = metaDataMassSpectrometryDf, file = fileMetaDataMassSpectrometryInformationProcessed, sep = "\t")
  #}
}
preprocessCompoundListFileList <- function(fileMetaDataCompoundInformationProcessed, fileCompoundListFileList){
  metaDataCompoundsDf <- read.table(file = fileMetaDataCompoundInformationProcessed, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  
  #########################################################################################################
  ## create CompoundList.csv / fileList.csv
  compoundListFileListDf <- metaDataCompoundsDf
  
  compoundListFileListDf$ID <- seq_len(nrow(compoundListFileListDf))
  compoundListFileListDf$Valid <- rep(x = FALSE, times = nrow(compoundListFileListDf))
  names(compoundListFileListDf)[names(compoundListFileListDf)=="RT (min)"] <- "RT"
  names(compoundListFileListDf)[names(compoundListFileListDf)=="File"] <- "Files"
  #compoundListFileListDf$RT <- trimsw(compoundListFileListDf$RT)
  for(cpdIdx in seq_len(nrow(compoundListFileListDf))){
    cat(c("\n", unname(unlist(compoundListFileListDf[cpdIdx, c("ID", "Name")]))))
    if(all(is.na(compoundListFileListDf$"InChI"[[cpdIdx]]), is.na(compoundListFileListDf$"SMILES"[[cpdIdx]]), is.na(compoundListFileListDf$"PubChem CID"[[cpdIdx]]))){
      stop(paste("\n### Error ### Compound with ID", cpdIdx, "has unknown structure"))
      #next
    }
    
    #if(is.na(compoundListFileListDf$SMILES[[cpdIdx]])) compoundListFileListDf$SMILES[[cpdIdx]] <- webchem::cs_inchi_smiles(inchi = compoundListFileListDf$InChI[[cpdIdx]], verbose = FALSE)
    #if(is.na(compoundListFileListDf$SMILES[[cpdIdx]])) cat(paste("\n### Warning ### Compound with ID", cpdIdx, "failed to convert to SMILES"))
    
    if(!grepl(x = compoundListFileListDf$RT[[cpdIdx]], pattern = "^\\d+([\\.,]\\d+)?")){
      cat(paste("\n### Warning ### Compound with ID", cpdIdx, "has invalid RT (", compoundListFileListDf$RT[[cpdIdx]], ")"))
      #next
    }
    
    #compoundListFileListDf$Valid[[cpdIdx]] <- TRUE
    compoundListFileListDf$Valid[[cpdIdx]] <- compoundListFileListDf$MspFileThere[[cpdIdx]]# & !is.na(compoundListFileListDf$InChI[[cpdIdx]])
    if(!compoundListFileListDf$Valid[[cpdIdx]])
      cat(paste("\n### Warning ### Compound with ID", cpdIdx, "has no corresponding spectrum file"))
  }
  
  compoundListFileListDf_valid <- compoundListFileListDf[compoundListFileListDf$Valid, ]
  
  cat(paste("\nWriting compound-list-file-list file", fileCompoundListFileList))
  write.table(x = compoundListFileListDf_valid, file = fileCompoundListFileList, sep = "\t")
}
removeItalicsTagsFromName <- function(name){
  #if(!grepl(x = name, pattern = "~\\{")){
  #  name <- NA
  #  next
  #  #stop(paste("invalid name", name))
  #}
  while(grepl(x = name, pattern = "~\\{")){
    list <- regexpr(text = name, pattern = "(?<italic>~\\{[^\\}]+\\})", perl = TRUE)
    name <- paste(
      substr(x = name, start = 1, stop = attr(list, "capture.start") - 1),
      substr(x = name, start = attr(list, "capture.start") + 2, stop = attr(list, "capture.start") + attr(list, "capture.length") - 2),
      substr(x = name, start = attr(list, "capture.start") + attr(list, "capture.length"), stop = nchar(name)),
      sep = ""
    )
  }
  return(name)
}
removeItalicsTagsFromNames <- function(names){
  if(any(grepl(x = names, pattern = "~"))){
    for(idx in which(grepl(x = names, pattern = "~")))
      names[[idx]] <- removeItalicsTagsFromName(names[[idx]])
    names <- names[!is.na(names)]
  }
  return(names)
}
preprocessInfoList <- function(fileCompoundListFileList, fileInfoList, takeRecordedNames){
  #########################################################################################################
  ## create infolist.csv
  
  ## columns:
  ## id	dbcas	dbname	dataused	COMMENT.CONFIDENCE	COMMENT.ID
  ## CH$NAME1	CH$NAME2	CH$NAME3	CH$COMPOUND_CLASS	CH$FORMULA	CH$EXACT_MASS	CH$SMILES	CH$IUPAC
  ## CH$LINK.CAS	CH$LINK.CHEBI	CH$LINK.HMDB	CH$LINK.KEGG	CH$LINK.LIPIDMAPS	CH$LINK.PUBCHEM	CH$LINK.INCHIKEY	CH$LINK.CHEMSPIDER
  
  compoundListFileListDf <- read.table(file = fileCompoundListFileList, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  
  infoListDf <- compoundListFileListDf
  
  names(infoListDf)[names(infoListDf)=="ID"  ] <- "id"
  infoListDf$"dbcas" <- rep(x = "", times = nrow(infoListDf))
  names(infoListDf)[names(infoListDf)=="Name"] <- "dbname"
  
  infoListDf$"dataused" <- rep(x = "TO_DO", times = nrow(infoListDf))
  names(infoListDf)[names(infoListDf)=="Confidence"] <- "COMMENT.CONFIDENCE"
  #infoListDf$COMMENT.ID
  infoListDf$"CH$NAME1" <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$NAME2" <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$NAME3" <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$NAME4" <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$NAME5" <- rep(x = "TO_DO", times = nrow(infoListDf))
  names(infoListDf)[names(infoListDf)=="Compound class"] <- "CH$COMPOUND_CLASS"
  infoListDf$"CH$FORMULA"    <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$EXACT_MASS" <- rep(x = "TO_DO", times = nrow(infoListDf))
  names(infoListDf)[names(infoListDf)=="SMILES"] <- "CH$SMILES"
  #infoListDf$"CH$SMILES"     <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$IUPAC"      <- rep(x = "TO_DO", times = nrow(infoListDf))
  
  infoListDf$"CH$LINK.CAS"        <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$LINK.CHEBI"      <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$LINK.HMDB"       <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$LINK.KEGG"       <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$LINK.LIPIDMAPS"  <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$LINK.PUBCHEM"    <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$LINK.INCHIKEY"   <- rep(x = "TO_DO", times = nrow(infoListDf))
  infoListDf$"CH$LINK.CHEMSPIDER" <- rep(x = "TO_DO", times = nrow(infoListDf))
  
  cat("\nProcessing InfoList cpd information using PubChem CID or InChI...")
  library("jsonlite")
  #library("Rdisop")
  library("rcdk")
  for(rowIdx in seq_len(nrow(compoundListFileListDf))){
    #for(rowIdx in 957:nrow(compoundListFileListDf)){
    pubchemCidThere <- !is.na(compoundListFileListDf$"PubChem CID"[[rowIdx]])
    inchiThere      <- !is.na(compoundListFileListDf$"InChI"      [[rowIdx]])
    smilesThere     <- !is.na(compoundListFileListDf$"SMILES"     [[rowIdx]])
    databaseLinksThere <- !is.na(compoundListFileListDf$"Database links"[[rowIdx]])
    
    smilesIsIsomeric <- ifelse(test = smilesThere, yes = any(c("/", "\\", "@") %in% strsplit(x = compoundListFileListDf$"SMILES"     [[rowIdx]], split = "")[[1]]), no = NA)
    
    dataused <- ifelse(
      test = pubchemCidThere, yes = "PubChem CID", no = ifelse(
      test = smilesThere    , yes = "SMILES",      no = ifelse(
      test = inchiThere     , yes = "InChI",       no = "NA"
    )))
    cat(paste("\nProcessing InfoList cpd", rowIdx, "/", nrow(compoundListFileListDf), "using", dataused))
    
    if(databaseLinksThere){
      databaseLinks <- compoundListFileListDf$"Database links"[[rowIdx]]
      databaseLinks <- strsplit(x = trimws(strsplit(x = databaseLinks, split = ";")[[1]]), split = "[ :]")
      databaseNames <- toupper(unlist(lapply(X = databaseLinks, FUN = "[", 1)))
      databaseIDs   <- unlist(lapply(X = databaseLinks, FUN = "[", 2))
    } else {
      databaseNames <- NA
      databaseIDs   <- NA
    }
    
    #givenName <- compoundListFileListDf$"Name"[[rowIdx]]
    if(takeRecordedNames)
      synonyms <- infoListDf$"dbname"[[rowIdx]]
    else
      synonyms <- NA
    if(all(!is.na(compoundListFileListDf$"Synonyms"[[rowIdx]]), compoundListFileListDf$"Synonyms"[[rowIdx]] != "") & takeRecordedNames)
      synonyms <- unique(trimws(c(synonyms, strsplit(x = compoundListFileListDf$"Synonyms"[[rowIdx]], split = "(; )|(, )")[[1]])))
    #else
    #  synonyms <- NA
    
    CH_NAME1 <- NA
    CH_NAME2 <- NA
    CH_NAME3 <- NA
    CH_NAME4 <- NA
    CH_NAME5 <- NA
    CH_NAME6 <- NA
    CH_FORMULA    <- "TO_DO"
    CH_EXACT_MASS <- "TO_DO"
    CH_SMILES     <- NA
    CH_IUPAC      <- NA
    
    CH_LINK.CAS        <- NA
    CH_LINK.CHEBI      <- NA
    CH_LINK.HMDB       <- NA
    CH_LINK.KEGG       <- NA
    CH_LINK.LIPIDMAPS  <- NA
    CH_LINK.PUBCHEM    <- NA
    CH_LINK.INCHIKEY   <- NA
    CH_LINK.CHEMSPIDER <- NA
    
    switch(dataused, 
           "PubChem CID"={
             ## PubChem CID available
             ## pubchemIds <- c(75215237, 46832173)
             ## pubchemIds <- c(9939941, 5280805)
             pubchemIds <- as.integer(strsplit(x = trimws(as.character(compoundListFileListDf$"PubChem CID"[[rowIdx]])), split = "; ")[[1]])
             
             ## xrefs
             ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON?inchi=InChI=1S/C3H8/c1-3-2/h3H2,1-2H3
             ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/AADVZSXPNRLYLV/JSON
             #url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/", inchi, "/JSON", sep = "")
             #url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON?inchi=", inchi, sep = "")
             #tmp <- fromJSON("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/AADVZSXPNRLYLV/JSON")[[1]]
             
             ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/CCCC/cids/TXT
             ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/CCCC/synonyms/XML
             ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/CCCC/xrefs/SBURL/XML
             ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/CCCC/xrefs/SBURL/JSON
             xrefs_m <- unique(unlist(sapply(X = pubchemIds, FUN = function(cid){
               ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/16219824/xrefs/SBURL/JSON
               url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid, "/xrefs/SBURL/JSON", sep = "")
               xrefs = tryCatch({    fromJSON(url)[[1]][[1]]$SBURL[[1]]  }, error = function(e) {  print(e); return(NA)  })
               if(all(!is.na(xrefs), xrefs==0)) xrefs <- NA
               return(xrefs)
             })))
             xrefs_m <- xrefs_m[!is.na(xrefs_m)]
             
             ## synonyms
             
             if(all(!is.na(synonyms), length(synonyms) > 0))
               synonyms_m <- getPubChemSynonyms(cids = pubchemIds, namesToRetain = synonyms)
             else
               synonyms_m <- getPubChemSynonyms(cids = pubchemIds)
             
             ## get records
             records <- sapply(X = pubchemIds, simplify = F, FUN = function(cid){
               recordHere <- NA
               while(length(recordHere)==1 && is.na(recordHere)){
                 recordHere <- tryCatch(expr = { webchem::pc_prop(cid = cid, verbose = F) },
                                        error = function(e){ if(grepl(x = e$message, pattern = "Too many requests or server too busy.")){print(e$message); Sys.sleep(time = 10); NA} else {stop(e$message)} }
                 )
                 return(recordHere)
               }
             })
             
             #record <- webchem::pc_prop(cid = pubchemIds[[1]], verbose = F)
             #smiles <- record$CanonicalSMILES
             #iupacName   <- record$IUPACName
             #monoMass <- record$MonoisotopicMass
             
             ## collect information
             CH_NAME1 <- ifelse(test = length(synonyms_m) > 0, yes = synonyms_m[[1]], no = NA)
             CH_NAME2 <- ifelse(test = length(synonyms_m) > 1, yes = synonyms_m[[2]], no = NA)
             CH_NAME3 <- ifelse(test = length(synonyms_m) > 2, yes = synonyms_m[[3]], no = NA)
             CH_NAME4 <- ifelse(test = length(synonyms_m) > 3, yes = synonyms_m[[4]], no = NA)
             CH_NAME5 <- ifelse(test = length(synonyms_m) > 4, yes = synonyms_m[[5]], no = NA)
             CH_NAME6 <- ifelse(test = length(synonyms_m) > 5, yes = synonyms_m[[6]], no = NA)
             CH_FORMULA    <- records[[1]]$MolecularFormula
             if(records[[1]]$Charge != 0){
               pattern <- ifelse(test = abs(records[[1]]$Charge) == 1, yes = ifelse(test = records[[1]]$Charge > 0, yes = "\\+", no = "\\-"), no = paste(abs(records[[1]]$Charge), ifelse(test = records[[1]]$Charge > 0, yes = "\\+", no = "\\-"), sep = ""))
               if(grepl(x = CH_FORMULA, pattern = paste(pattern, "$", sep = ""))){
                 CH_FORMULA <- gsub(x = CH_FORMULA, pattern = paste(pattern, "$", sep = ""), replacement = paste("]", pattern, sep = ""))
                 CH_FORMULA <- paste("[", CH_FORMULA, sep = "")
               }
             }
             CH_EXACT_MASS <- records[[1]]$MonoisotopicMass
             #CH_SMILES     <- records[[1]]$CanonicalSMILES
             CH_SMILES     <- records[[1]]$IsomericSMILES
             CH_IUPAC      <- ifelse(test = length(unlist(lapply(X = records, FUN = function(x){x$IUPACName})))>0, yes = lapply(X = records, FUN = function(x){x$IUPACName})[[1]], no = NA)
             
             CH_LINK.CAS.url        <- xrefs_m[grepl(x = xrefs_m, pattern = "^https://search\\.toxplanet\\.com/CategorySearch\\.aspx\\?cas_no=")]
             CH_LINK.CHEBI.url      <- xrefs_m[grepl(x = xrefs_m, pattern = "^http://www\\.ebi\\.ac\\.uk/chebi/searchId\\.do\\?chebiId=")]
             CH_LINK.HMDB.url       <- xrefs_m[grepl(x = xrefs_m, pattern = "^http://www\\.hmdb\\.ca/metabolites/")]
             CH_LINK.KEGG.url       <- xrefs_m[grepl(x = xrefs_m, pattern = "^http://www\\.genome\\.jp/dbget-bin/www_bget\\?cpd:")]
             CH_LINK.LIPIDMAPS.url  <- xrefs_m[grepl(x = xrefs_m, pattern = "^http://www\\.lipidmaps\\.org/data/LMSDRecord\\.php\\?LM_ID=")]
             CH_LINK.CHEMSPIDER.url <- xrefs_m[grepl(x = xrefs_m, pattern = "^http://www\\.chemspider\\.com/Chemical-Structure\\.")]
             
             CH_LINK.CAS.url        <- ifelse(test = length(CH_LINK.CAS.url)        > 0, yes = CH_LINK.CAS.url[[1]],        no = "")
             CH_LINK.CHEBI.url      <- ifelse(test = length(CH_LINK.CHEBI.url)      > 0, yes = CH_LINK.CHEBI.url[[1]],      no = "")
             CH_LINK.HMDB.url       <- ifelse(test = length(CH_LINK.HMDB.url)       > 0, yes = CH_LINK.HMDB.url[[1]],       no = "")
             CH_LINK.KEGG.url       <- ifelse(test = length(CH_LINK.KEGG.url)       > 0, yes = CH_LINK.KEGG.url[[1]],       no = "")
             CH_LINK.LIPIDMAPS.url  <- ifelse(test = length(CH_LINK.LIPIDMAPS.url)  > 0, yes = CH_LINK.LIPIDMAPS.url[[1]],  no = "")
             CH_LINK.CHEMSPIDER.url <- ifelse(test = length(CH_LINK.CHEMSPIDER.url) > 0, yes = CH_LINK.CHEMSPIDER.url[[1]], no = "")
             
             ## xrefs_m[grepl(x = xrefs_m, pattern = "153-18-4")]
             ## https://search.toxplanet.com/CategorySearch.aspx?cas_no=153-18-4
             ## "^\\d{1,10}-\\d{2,2}-\\d$" 153-18-4
             ## xrefs_m[grepl(x = xrefs_m, pattern = "^https://search\\.toxplanet\\.com/CategorySearch\\.aspx\\?cas_no=")]
             CH_LINK.CAS        <- unlist(lapply(X = strsplit(x = CH_LINK.CAS.url, split = "="), FUN = "[", 2))[[1]]
             ## http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:28527
             ## "^CHEBI:\\d{5,5}$" CHEBI:28527
             CH_LINK.CHEBI      <- unlist(lapply(X = strsplit(x = CH_LINK.CHEBI.url, split = "="), FUN = "[", 2))[[1]]
             ## http://www.hmdb.ca/metabolites/HMDB0003249
             ## "^HMDB\\d{7,7}$" C05625
             CH_LINK.HMDB       <- unlist(lapply(X = strsplit(x = CH_LINK.HMDB.url, split = "/"), FUN = "[", 5))[[1]]
             ## http://www.genome.jp/dbget-bin/www_bget?cpd:C05625
             ## "^C\\d{5,5}$" C05625
             CH_LINK.KEGG       <- unlist(lapply(X = strsplit(x = CH_LINK.KEGG.url, split = ":"), FUN = "[", 3))[[1]]
             ## http://www.lipidmaps.org/data/LMSDRecord.php?LM_ID=LMSP03010003
             ## "^LMSP\\d{8,8}$" LMSP03010003
             CH_LINK.LIPIDMAPS  <- unlist(lapply(X = strsplit(x = CH_LINK.LIPIDMAPS.url, split = "="), FUN = "[", 2))[[1]]
             ## trivial
             CH_LINK.PUBCHEM    <- pubchemIds[[1]]
             ## record
             CH_LINK.INCHIKEY   <- records[[1]]$InChIKey
             ## http://www.chemspider.com/Chemical-Structure.4444362.html
             CH_LINK.CHEMSPIDER <- unlist(lapply(X = strsplit(x = CH_LINK.CHEMSPIDER.url, split = "\\."), FUN = "[", 4))[[1]]
             
             CH_LINK.CAS        <- ifelse(test = !is.null(CH_LINK.CAS),        yes = CH_LINK.CAS,        no = "NA")
             CH_LINK.CHEBI      <- ifelse(test = !is.null(CH_LINK.CHEBI),      yes = CH_LINK.CHEBI,      no = "NA")
             CH_LINK.HMDB       <- ifelse(test = !is.null(CH_LINK.HMDB),       yes = CH_LINK.HMDB,       no = "NA")
             CH_LINK.KEGG       <- ifelse(test = !is.null(CH_LINK.KEGG),       yes = CH_LINK.KEGG,       no = "NA")
             CH_LINK.LIPIDMAPS  <- ifelse(test = !is.null(CH_LINK.LIPIDMAPS),  yes = CH_LINK.LIPIDMAPS,  no = "NA")
             #CH_LINK.PUBCHEM    <- 
             #CH_LINK.INCHIKEY   <- 
             CH_LINK.CHEMSPIDER <- ifelse(test = !is.null(CH_LINK.CHEMSPIDER), yes = CH_LINK.CHEMSPIDER, no = "NA")
           },
           "InChI"={
             ## no PubChem CID available
             ### get formula and exect mass from structure
             ### InChI=1S/C20H30O2/c1-12(2)14-11-13-7-8-15-19(3,4)9-6-10-20(15,5)16(13)18(22)17(14)21/h11-12,15,21-22H,6-10H2,1-5H3/t15-,20-/m0/s1
             #CH_FORMULA <- strsplit(x = strsplit(x = compoundListFileListDf$InChI[[rowIdx]], split = "=")[[1]][[2]], split = "/")[[1]][[2]]
             ### Rdisop
             #CH_EXACT_MASS <- Rdisop::getMass(molecule = Rdisop::getMolecule(formula = CH_FORMULA))
             
             molecule <- rinchi::parse.inchi(inchis = compoundListFileListDf$"InChI"[[rowIdx]])[[1]]
             mol <- rcdk::get.mol2formula(molecule = molecule)
             CH_FORMULA <- mol@string
             CH_EXACT_MASS <- mol@mass
             CH_SMILES <- ifelse(test = all(smilesThere, smilesIsIsomeric),
                                 yes = compoundListFileListDf$"SMILES"     [[rowIdx]],
                                 no = rcdk::get.smiles(molecule = molecule, flavor = smiles.flavors("Isomeric"))
             )
             
             if(all(!is.na(synonyms), length(synonyms) > 0)) CH_NAME1 <- synonyms[[1]]
             if(all(!is.na(synonyms), length(synonyms) > 1)) CH_NAME2 <- synonyms[[2]]
             if(all(!is.na(synonyms), length(synonyms) > 2)) CH_NAME3 <- synonyms[[3]]
             if(all(!is.na(synonyms), length(synonyms) > 3)) CH_NAME4 <- synonyms[[4]]
             if(all(!is.na(synonyms), length(synonyms) > 4)) CH_NAME5 <- synonyms[[5]]
             if(all(!is.na(synonyms), length(synonyms) > 5)) CH_NAME6 <- synonyms[[6]]
           },
           "SMILES"={
             ## no PubChem CID and no InChI available
             ## C=12[C@]3(C)[C@@]([H])(CCC1C=C(C(=C2O)O)C(C)C)C(CCC3)(C)C
             molecule <- rcdk::parse.smiles(smiles = compoundListFileListDf$"SMILES"[[rowIdx]])[[1]]
             mol <- rcdk::get.mol2formula(molecule = molecule)
             CH_FORMULA <- mol@string
             CH_EXACT_MASS <- mol@mass
             CH_SMILES <- compoundListFileListDf$"SMILES"[[rowIdx]]
             
             if(all(!is.na(synonyms), length(synonyms) > 0)) CH_NAME1 <- synonyms[[1]]
             if(all(!is.na(synonyms), length(synonyms) > 1)) CH_NAME2 <- synonyms[[2]]
             if(all(!is.na(synonyms), length(synonyms) > 2)) CH_NAME3 <- synonyms[[3]]
             if(all(!is.na(synonyms), length(synonyms) > 3)) CH_NAME4 <- synonyms[[4]]
             if(all(!is.na(synonyms), length(synonyms) > 4)) CH_NAME5 <- synonyms[[5]]
             if(all(!is.na(synonyms), length(synonyms) > 5)) CH_NAME6 <- synonyms[[6]]
           },
           "NA"={
             stop(paste("Unknown structure"))
           },
           {  stop(paste("Unknown dataused '", dataused, "'"))  }
    )
    
    ## take over present database links if there is nothing there
    if(is.na(CH_LINK.CAS)        & "CAS"        %in% databaseNames) CH_LINK.CAS        <- databaseIDs[databaseNames == "CAS"]       [[1]]
    if(is.na(CH_LINK.CHEBI)      & "CHEBI"      %in% databaseNames) CH_LINK.CHEBI      <- databaseIDs[databaseNames == "CHEBI"]     [[1]]
    if(is.na(CH_LINK.HMDB)       & "HMDB"       %in% databaseNames) CH_LINK.HMDB       <- databaseIDs[databaseNames == "HMDB"]      [[1]]
    if(is.na(CH_LINK.KEGG)       & "KEGG"       %in% databaseNames) CH_LINK.KEGG       <- databaseIDs[databaseNames == "KEGG"]      [[1]]
    if(is.na(CH_LINK.LIPIDMAPS)  & "LIPIDMAPS"  %in% databaseNames) CH_LINK.LIPIDMAPS  <- databaseIDs[databaseNames == "LIPIDMAPS"] [[1]]
    if(is.na(CH_LINK.PUBCHEM)    & "PUBCHEM"    %in% databaseNames) CH_LINK.PUBCHEM    <- databaseIDs[databaseNames == "PUBCHEM"]   [[1]]
    if(is.na(CH_LINK.INCHIKEY)   & "INCHIKEY"   %in% databaseNames) CH_LINK.INCHIKEY   <- databaseIDs[databaseNames == "INCHIKEY"]  [[1]]
    if(is.na(CH_LINK.CHEMSPIDER) & "CHEMSPIDER" %in% databaseNames) CH_LINK.CHEMSPIDER <- databaseIDs[databaseNames == "CHEMSPIDER"][[1]]
    
    CH_NAME1 <- gsub(x = CH_NAME1, pattern = "\"", replacement = "'")
    CH_NAME2 <- gsub(x = CH_NAME2, pattern = "\"", replacement = "'")
    CH_NAME3 <- gsub(x = CH_NAME3, pattern = "\"", replacement = "'")
    CH_NAME4 <- gsub(x = CH_NAME4, pattern = "\"", replacement = "'")
    CH_NAME5 <- gsub(x = CH_NAME5, pattern = "\"", replacement = "'")
    CH_NAME6 <- gsub(x = CH_NAME6, pattern = "\"", replacement = "'")
    
    if(!takeRecordedNames & !is.na(CH_NAME1)){
      infoListDf$"dbname"[[rowIdx]] <- CH_NAME1
      CH_NAME1 <- CH_NAME2
      CH_NAME2 <- CH_NAME3
      CH_NAME3 <- CH_NAME4
      CH_NAME4 <- CH_NAME5
      CH_NAME5 <- CH_NAME6
    }
    
    ## box
    infoListDf$dataused  [[rowIdx]] <- dataused
    infoListDf$"CH$NAME1"[[rowIdx]] <- CH_NAME1
    infoListDf$"CH$NAME2"[[rowIdx]] <- CH_NAME2
    infoListDf$"CH$NAME3"[[rowIdx]] <- CH_NAME3
    infoListDf$"CH$NAME4"[[rowIdx]] <- CH_NAME4
    infoListDf$"CH$NAME5"[[rowIdx]] <- CH_NAME5
    
    infoListDf$"CH$FORMULA"   [[rowIdx]] <- CH_FORMULA
    infoListDf$"CH$EXACT_MASS"[[rowIdx]] <- CH_EXACT_MASS
    infoListDf$"CH$SMILES"    [[rowIdx]] <- CH_SMILES
    infoListDf$"CH$IUPAC"     [[rowIdx]] <- CH_IUPAC
    
    infoListDf$"CH$LINK.CAS"       [[rowIdx]] <- CH_LINK.CAS
    infoListDf$"CH$LINK.CHEBI"     [[rowIdx]] <- CH_LINK.CHEBI
    infoListDf$"CH$LINK.HMDB"      [[rowIdx]] <- CH_LINK.HMDB
    infoListDf$"CH$LINK.KEGG"      [[rowIdx]] <- CH_LINK.KEGG
    infoListDf$"CH$LINK.LIPIDMAPS" [[rowIdx]] <- CH_LINK.LIPIDMAPS
    infoListDf$"CH$LINK.PUBCHEM"   [[rowIdx]] <- CH_LINK.PUBCHEM
    infoListDf$"CH$LINK.INCHIKEY"  [[rowIdx]] <- CH_LINK.INCHIKEY
    infoListDf$"CH$LINK.CHEMSPIDER"[[rowIdx]] <- CH_LINK.CHEMSPIDER
  }
  
  cat(paste("\nWriting info-list file", fileInfoList))
  write.table(x = infoListDf, file = fileInfoList, sep = "\t")
  if(any(apply(X = infoListDf, MARGIN = 2, FUN = function(x){sum(grepl(x = x, pattern = "\""))>0}))) stop("infoListDf with quotes")
  
  
  ## update SMILES
  compoundListFileListDf$SMILES <- infoListDf$`CH$SMILES`
  write.table(x = compoundListFileListDf, file = fileCompoundListFileList, sep = "\t")
}
adductToIonMode <- function(adduct){
  if(is.na(adduct)) return("")
  ## added:
  ## pK [M+K]+; 
  ## pACN_pH [M+ACN+H]+; 
  ## pACN_pNa [M+ACN+Na]+;  ## Acetonitrile C2H3N
  ## p2Na_mH [M+2Na-H]+; 
  
  ## pM_pH [2M+H]+; 
  ## pM_pK [2M+K]+; 
  ## pM_pNa [2M+Na]+; 
  ## pM_pNH4 [2M+NH4]+; 
  ## pM_pACN_pH [2M+ACN+H]+; 
  
  ## pACN_p2H [M+ACN+2H]2+; 
  ## p2H [M+2H]2+
  
  ## findMz.formula("C6H12O6", "")
  adducts <- getAdductInformation("")[, c("mode", "adductString")]
  if(!(adduct %in% adducts$adductString))
    stop(paste("Unknown adduct '", adduct, "'", sep = ""))
  ionMode <- adducts$mode[adducts$adductString==adduct]
  return(ionMode)
}
preprocessSettings <- function(fileCompoundListFileList_part, fileSettings_part, accessionPrefix, ms_type){
  #########################################################################################################
  ## create settings.ini
  
  metaDataCompoundsDf_part <- read.csv(file = fileCompoundListFileList_part, stringsAsFactors = FALSE, check.names = FALSE)
  #metaDataCompoundsDf_part <- read.table(file = fileCompoundListFileList_part, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  
  ## static: TODO store in extra settings file or discuss with Gerd as being fixed
  #ms_type <- "MS2"
  #accessionPrefix <- "IH"
  fragmentation_mode <- "CID"
  
  ## settings...
  openBabelPath <- paste(dirname(system(command = "which obabel", intern = TRUE)), "/", sep = "")
  authors <- as.character(metaDataCompoundsDf_part$Authors[[1]])
  copyright <- paste("Copyright (C)", authors)
  publication <- metaDataCompoundsDf_part$Publication[[1]]
  if(any(is.na(metaDataCompoundsDf_part$Publication[[1]]), is.null(metaDataCompoundsDf_part$Publication[[1]]), metaDataCompoundsDf_part$Publication[[1]] == ""))
    publication <- NULL
  else
    publication <- gsub(x = publication, pattern = ":", replacement = ";")
  instrument <- as.character(metaDataCompoundsDf_part$INSTRUMENT[[1]])
  instrument_type <- as.character(metaDataCompoundsDf_part$INSTRUMENT_TYPE[[1]])
  instrument_configuration <- NULL
  confidence_comment <- as.character(metaDataCompoundsDf_part$Confidence[[1]])
  comment <- NULL
  compound_class <- metaDataCompoundsDf_part$"Compound class"[[1]]
  internal_id_fieldname <- "INTERNAL_ID"
  lc_column <- as.character(metaDataCompoundsDf_part$CHROMATOGRAPHY[[1]])
  ionization <- as.character(metaDataCompoundsDf_part$IONIZATION[[1]])
  ce <- metaDataCompoundsDf_part$"Collision energy"[[1]]
  ces <- ce
  res <- NULL
  
  columnsMassSpectrometry <- colnames(metaDataCompoundsDf_part)[grepl(x = colnames(metaDataCompoundsDf_part), pattern = "^MS_SETTING_")]
  columnsChromatography   <- colnames(metaDataCompoundsDf_part)[grepl(x = colnames(metaDataCompoundsDf_part), pattern = "^CHROMATOGRAPHY_SETTING_")]
  tagMassSpectrometry     <- gsub(x = columnsMassSpectrometry, pattern = "^MS_SETTING_",             replacement = "AC$MASS_SPECTROMETRY_")
  tagChromatography       <- gsub(x = columnsChromatography,   pattern = "^CHROMATOGRAPHY_SETTING_", replacement = "AC$CHROMATOGRAPHY_")
  valueMassSpectrometry   <- unlist(lapply(X = columnsMassSpectrometry, FUN = function(column){metaDataCompoundsDf_part[[1,column]]}))
  valueChromatography     <- unlist(lapply(X = columnsChromatography,   FUN = function(column){metaDataCompoundsDf_part[[1,column]]}))
  
  valueMassSpectrometry   <- gsub(x = valueMassSpectrometry, pattern = " : ", replacement = ":")
  valueChromatography     <- gsub(x = valueChromatography, pattern = " : ", replacement = ":")
  valueMassSpectrometry   <- gsub(x = valueMassSpectrometry, pattern = ": ", replacement = ":")
  valueChromatography     <- gsub(x = valueChromatography, pattern = ": ", replacement = ":")
  
  ################################################################
  ## collect information
  listSettings <- list()
  # Deprofile input data?
  listSettings$"deprofile" <- ""
  
  # Deviation (in minutes) allowed the for retention time
  listSettings$"rtMargin" <- 0.2
  # Systematic retention time shift
  listSettings$"rtShift" <- 0
  
  # Points to the directory where babel.exe (or the Linux "babel" equivalent) lies.
  listSettings$"babeldir" <- openBabelPath
  
  # Which MassBank record version to use; version 2 is advised.
  listSettings$"use_version" <- 2
  
  # Include reanalyzed peaks? TODO necessary?
  listSettings$"use_rean_peaks" <- TRUE
  
  # annotate the spectra files with (putative) molecular formulas for fragments?
  listSettings$"add_annotation" <- TRUE
  
  # Annotations for the spectrum:
  listSettings$"annotations" <- ""
  # Author etc. annotation
  ## AUTHORS
  listSettings$"    authors" <- authors
  ## COPYRIGHT
  listSettings$"    copyright" <- copyright
  ## PUBLICATION
  if(!is.null(publication)) listSettings$"    publication" <- publication
  ## LICENSE
  listSettings$"    license" <- "CC BY"
  ## AC$INSTRUMENT
  listSettings$"    instrument" <- instrument
  ## AC$INSTRUMENT_TYPE
  listSettings$"    instrument_type" <- instrument_type
  ## ignored?
  if(!is.null(instrument_configuration)) listSettings$"    instrument_configuration" <- instrument_configuration
  ## COMMENT: CONFIDENCE
  listSettings$"    confidence_comment" <- confidence_comment
  ## ignored?
  if(!is.null(comment)) listSettings$"    comment" <- comment
  ## CH$COMPOUND_CLASS
  listSettings$"    compound_class" <- compound_class
  ## ignored?
  if(!is.null(internal_id_fieldname)) listSettings$"    internal_id_fieldname" <- internal_id_fieldname
  
  # HPLC annotations:
  #if(FALSE){
  ## example: lc_gradient: 90/10 at 0 min, 50/50 at 4 min, 5/95 at 17 min, 5/95 at 25 min, 90/10 at 25.1 min, 90/10 at 30 min
  ### AC$CHROMATOGRAPHY: FLOW_GRADIENT
  #lc_gradient: 98/2 at 0 min, 98/2 at 2 min, 64/36 at 18 min, 5/95 at 21 min, 5/95 at 22.5 min, 98/2 at 22.52 min, 98/2 at 24 min
  ## example: lc_flow: 200 uL/min
  ### AC$CHROMATOGRAPHY: FLOW_RATE
  #lc_flow: 0.4 mL/min
  ## example: lc_solvent_a: water with 0.1% formic acid
  ### AC$CHROMATOGRAPHY: SOLVENT A
  #lc_solvent_a: 0.33 mM ammonium formate with 0.66 mM formic acid in water
  ### AC$CHROMATOGRAPHY: SOLVENT B
  #lc_solvent_b: acetonitrile
  #}
  # example: lc_column: XBridge C18 3.5um, 2.1x50mm, Waters
  ## AC$CHROMATOGRAPHY: COLUMN_NAME
  listSettings$"    lc_column" <- lc_column
  
  for(idx in seq_along(tagMassSpectrometry))  listSettings[[paste("    ", tagMassSpectrometry[[idx]], sep = "")]] <- valueMassSpectrometry[[idx]]
  
  # Prefix for MassBank accession IDs
  listSettings$"    entry_prefix" <- accessionPrefix
  ## AC$MASS_SPECTROMETRY: MS_TYPE
  listSettings$"    ms_type" <- ms_type
  ## AC$MASS_SPECTROMETRY: IONIZATION
  listSettings$"    ionization" <- ionization
  #if(FALSE){
  ### MS$DATA_PROCESSING: RECALIBRATE
  #listSettings$"    ms_dataprocessing" <- ""
  #ms_dataprocessing:
  #  RECALIBRATE: APCI calibrant solution via a calibration delivery system (CDS)
  #listSettings$"        RECALIBRATE" <- ""
  #}
  
  for(idx in seq_along(tagChromatography  ))  listSettings[[paste("    ", tagChromatography  [[idx]], sep = "")]] <- valueChromatography  [[idx]]
  
  # Annotator:
  listSettings$"spectraList" <- ""
  ## AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE
  # mode: fragmentation mode, e.g. CID
  listSettings$"- mode" <- fragmentation_mode
  ## AC$MASS_SPECTROMETRY: COLLISION_ENERGY
  # ces: "short" format collision energy (for record title)
  # ce: "long" format collision energy (for annotation field)
  listSettings$"  ces" <- ces
  listSettings$"  ce" <- ce
  ## AC$MASS_SPECTROMETRY: RESOLUTION
  # res: FT resolution
  if(!is.null(res)) listSettings$"  res" <- res
  
  # Shifts of the starting points for RMassBank accession numbers.
  listSettings$"accessionNumberShifts" <- ""
  listSettings$"    pH"  <-  0  # [M+H]+: Accession numbers 1-14
  listSettings$"    pM"  <- 16  # [M]+: 17-30
  listSettings$"    pNa" <- 32  # [M+Na]+: 33-46
  listSettings$"    mH"  <- 50  # [M-H]-: 51-64
  listSettings$"    mFA" <- 66  # [M+FA]-: 67-80
  
  ionModes <- sapply(X = unique(metaDataCompoundsDf_part$Adduct), FUN = adductToIonMode)
  presentIonModes <- c("pH", "pM", "pNa", "mH", "mFA")
  if(!all(ionModes %in% presentIonModes)){
    missingIonModes <- ionModes[!(ionModes %in% presentIonModes)]
    missingIonModes <- missingIonModes[missingIonModes!=""]
    cat(paste("\n### Warning ### Adding unusual adducts", paste(missingIonModes, collapse = "; ")))
    for(ionMode in missingIonModes)
      listSettings[[paste("    ", ionMode, sep = "")]]  <-  0
  }
  
  ## TODO discuss with Gerd
  ## A list of known electronic noise peaks 
  listSettings$"electronicNoise" <- paste(
    "",
    #"- 189.825",
    #"- 201.725",
    #"- 196.875",
    sep = "\n"
  )
  ## Exclusion width of electronic noise peaks (from unmatched peaks, prior to reanalysis)
  listSettings$"electronicNoiseWidth" <- 0.3
  
  # recalibration settings:
  listSettings$"recalibrateBy" <- "dppm"
  
  # recalibrate MS1:
  listSettings$"recalibrateMS1" <- "common"
  # Window width to look for MS1 peaks to recalibrate (in ppm)
  listSettings$"recalibrateMS1Window" <- 15
  
  # Custom recalibration function: You can overwrite the recal function by
  listSettings$"recalibrator" <- ""
  listSettings$"    MS1" <- "recalibrate.identity"
  listSettings$"    MS2" <- "recalibrate.identity"
  
  # Define the multiplicity filtering level
  listSettings$"multiplicityFilter" <- 1
  
  # Define the title format.
  listSettings$"titleFormat" <- paste(
    "",
    "- \"{CH$NAME}\"",
    "- \"{AC$INSTRUMENT_TYPE}\"",
    "- \"{AC$MASS_SPECTROMETRY: MS_TYPE}\"",
    "- \"CE: {RECORD_TITLE_CE}\"",
    #"- \"R={AC$MASS_SPECTROMETRY: RESOLUTION}\"", ## TODO
    "- \"{MS$FOCUSED_ION: PRECURSOR_TYPE}\"",
    sep = "\n"
  )
  
  # Define filter settings.
  ## TODO discuss with Gerd
  listSettings$"filterSettings" <- ""
  listSettings$"    ppmHighMass"       <-  60
  listSettings$"    ppmLowMass"        <- 120
  listSettings$"    massRangeDivision" <- 150
  listSettings$"    ppmFine"           <-  60
  listSettings$"    prelimCut"         <- 000
  listSettings$"    prelimCutRatio"    <-   0
  listSettings$"    fineCut"           <-   0
  listSettings$"    fineCutRatio"      <-   0
  listSettings$"    specOkLimit"       <- 000
  listSettings$"    dbeMinLimit"       <-  -0.5
  listSettings$"    satelliteMzLimit"  <-   0.
  listSettings$"    satelliteIntLimit" <-   0.0
  
  # Define raw MS retrieval settings.
  listSettings$"findMsMsRawSettings" <- ""
  listSettings$"    ppmFine" <- 10
  listSettings$"    mzCoarse" <- 0.5
  
  # fillPrecursorScan is FALSE for "good" mzML files which have all the info needed.
  # However, for example AB Sciex files will have missing precursor scan information,
  # in which case fillPrecursorScan = TRUE is needed. Try it out.
  listSettings$"    fillPrecursorScan" <- FALSE ## TODO adapt if needed
  
  
  ################################################################
  ## box and write information
  fileLines <- paste(names(listSettings), unname(listSettings), sep = ": ")
  cat(paste("\nWriting settings file", fileSettings_part))
  writeLines(text = fileLines, con = fileSettings_part)
}

preprocessContributorToMassBankWorkflow <- function(folder, accessionPrefix, xlsxFile, ms_type, takeRecordedNames, reprocessMetaDataCompoundInformation = FALSE, reprocessCompoundList_FileList = FALSE, reprocessInfoList = FALSE, reprocessSubFolders = FALSE){
  if(reprocessMetaDataCompoundInformation){
    reprocessCompoundList_FileList <- TRUE
    reprocessInfoList <- TRUE
    reprocessSubFolders <- TRUE
  }
  if(reprocessCompoundList_FileList){
    reprocessInfoList <- TRUE
    reprocessSubFolders <- TRUE
  }
  if(reprocessInfoList){
    reprocessSubFolders <- TRUE
  }
  
  if(accessionPrefix %in% knownAccessionPrefixes) stop(paste("Accession prefix", accessionPrefix, "is already known."))
  
  #folder <- "/mnt/data/IPB/Projects/2017_005_MS-databases/mFam contributions/Ulschan Bathe di and sesquiterpenes"
  #folderClassifier <- paste(folder, "classificator", sep = "/")
  folderMsp        <- paste(folder, "converted to msp", sep = "/")
  folderMassBank   <- paste(folder, "converted to MassBank", sep = "/")
  folderMetaData   <- paste(folder, "meta data", sep = "/")
  folderRawData    <- paste(folder, "raw data", sep = "/")
  folderRMassBankData <- paste(folder, "RMassBank", sep = "/")
  if(!file.exists(folderRMassBankData)) if(!dir.create(folderRMassBankData))  stop(paste("Could not create folder", folderRMassBankData))
  
  if(!all(file.exists(folderMsp, folderMassBank, folderMetaData, folderRawData)))
    stop(paste("Folder structure is not conventional. Missing folder", c(folderMsp, folderMassBank, folderMetaData, folderRawData)[file.exists(folderMsp, folderMassBank, folderMetaData, folderRawData)]))
  cat("\nFolder structure is conventional")
  
  ##################################################################################
  ## handle files
  fileMetaData <- list.files(path = folderMetaData, pattern = "^[^~]*\\.xlsx$", full.names = TRUE, recursive = TRUE)
  #fileMetaData <- file.choose()
  #if(!length(fileMetaData) == 1)  stop(paste("Not exactly one meta data file (", length(fileMetaData), "): ", paste(basename(fileMetaData), collapse = "; "), sep = ""))
#------  if(!(xlsxFile %in% basename(fileMetaData)))  stop(paste("Xlsx file not there (", length(fileMetaData), "): ", paste(basename(fileMetaData), collapse = "; "), sep = ""))
  fileMetaData <- fileMetaData[grepl(x = fileMetaData, pattern = paste(xlsxFile, "$", sep = ""))]
  
  fileMetaDataCompoundInformationProcessed         <- paste(folderRMassBankData, "MetaDataCompoundInformationProcessed.csv", sep = "/")
  #fileMetaDataChromatographyInformationProcessed   <- paste(folderRMassBankData, "MetaDataChromatographyInformationProcessed.csv", sep = "/")
  #fileMetaDataMassSpectrometryInformationProcessed <- paste(folderRMassBankData, "MetaDataMassSpectrometryInformationProcessed.csv", sep = "/")
  
  ## process given meta data
  cat(paste("\n", ifelse(test = file.exists(fileMetaDataCompoundInformationProcessed), yes = "exists: ", no = "is being created: "), fileMetaDataCompoundInformationProcessed, sep = ""))
  if(!file.exists(fileMetaDataCompoundInformationProcessed) | reprocessMetaDataCompoundInformation)
    preprocessData(folder, fileMetaData, fileMetaDataCompoundInformationProcessed, ms_type)#, fileMetaDataChromatographyInformationProcessed, fileMetaDataMassSpectrometryInformationProcessed)
  #metaDataCompoundsDf <- read.table(file = fileMetaDataCompoundInformationProcessed, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  
  #chromatographyInformationCanByUsed   <- file.exists(fileMetaDataChromatographyInformationProcessed)
  #massSpectrometryInformationCanByUsed <- file.exists(fileMetaDataMassSpectrometryInformationProcessed)
  #if(chromatographyInformationCanByUsed)   metaDataChromatographyDf   <- read.table(file = fileMetaDataCompoundInformationProcessed,         header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  #if(massSpectrometryInformationCanByUsed) metaDataMassSpectrometryDf <- read.table(file = fileMetaDataMassSpectrometryInformationProcessed, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  
  fileCompoundListFileList <- paste(folderRMassBankData, "CompoundList_FileList.csv", sep = "/")
  fileInfoList             <- paste(folderRMassBankData, "InfoList.csv", sep = "/")
  #fileSettings             <- paste(folderRMassBankData, "Settings.ini", sep = "/")
  
  ## compound list
  cat(paste("\n", ifelse(test = file.exists(fileCompoundListFileList), yes = "exists: ", no = "is being created :"), fileCompoundListFileList, sep = ""))
  if(!file.exists(fileCompoundListFileList) | reprocessCompoundList_FileList)
    preprocessCompoundListFileList(fileMetaDataCompoundInformationProcessed, fileCompoundListFileList)
  #compoundListFileListDf <- read.table(file = fileCompoundListFileList, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  
  ## info list
  cat(paste("\n", ifelse(test = file.exists(fileInfoList), yes = "exists: ", no = "is being created: "), fileInfoList, sep = ""))
  if(!file.exists(fileInfoList) | reprocessInfoList)
    preprocessInfoList(fileCompoundListFileList, fileInfoList, takeRecordedNames)
  ##infoListDf <- read.table(file = fileInfoList, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  
  ####################################################################################################
  ## split into setting-wise parts
  compoundListFileListDf <- read.table(file = fileCompoundListFileList, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  infoListDf             <- read.table(file = fileInfoList,             header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  
  compoundListFileListDf[is.na(compoundListFileListDf)] <- "NA" ## for splitting
  
  cat("\n\nSplitting compound list")
  columnsMassSpectrometry <- colnames(compoundListFileListDf)[grepl(x = colnames(compoundListFileListDf), pattern = "^MS_SETTING_")]
  columnsChromatography   <- colnames(compoundListFileListDf)[grepl(x = colnames(compoundListFileListDf), pattern = "^CHROMATOGRAPHY_SETTING_")]
  splitColumns <- c(
    "Authors",
    "Publication",
    "INSTRUMENT",
    "INSTRUMENT_TYPE",
    "Confidence",
    "Compound class",
    "CHROMATOGRAPHY",
    "IONIZATION",
    "Adduct",
    "Ionization mode",
    "Collision energy",
    columnsMassSpectrometry,
    columnsChromatography
  )
  
  if(FALSE){
    splitList <- list()
    for(splitColumn in splitColumns){
      cat(paste("\nColumn '", splitColumn, "': ", length(unique(compoundListFileListDf[, splitColumn])), " values", sep = ""))
      splitList[[splitColumn]] <- compoundListFileListDf[, splitColumn]
    }
    
    #split(x = compoundListFileListDf, f = list(compoundListFileListDf$"Ionization mode", compoundListFileListDf$"Compound class"))
    compoundListFileListDf_splitted <- split(x = compoundListFileListDf, f = splitList)
    compoundListFileListDf_splitted <- compoundListFileListDf_splitted[unlist(lapply(X = compoundListFileListDf_splitted, FUN = nrow))>0]
  }
  
  for(splitColumn in splitColumns) {
    values <- unique(compoundListFileListDf[, splitColumn])
    if(length(values) == 1) next
    cat(paste("\nColumn '", splitColumn, "': ", length(values), " values:\n", sep = ""))
    cat(paste(paste("   ", values), collapse = "\n"))
  }
  
  splitDf <- unique(compoundListFileListDf[, splitColumns])
  compoundListFileListDf_splitted <- apply(X = splitDf, MARGIN = 1, FUN = function(row){
    rowSelection <- apply(X = compoundListFileListDf[, splitColumns], MARGIN = 1, FUN = function(row2){all(row==row2)})
    compoundListFileListDf[rowSelection, , drop=FALSE]
  })
  
  maximumNumberOfSpectraPerBatch <- 50
  compoundListFileListDf_splitted <- unlist(lapply(X = compoundListFileListDf_splitted, FUN = function(compoundListFileListDf_splittedHere){
    numberOfBatches <- ceiling(nrow(compoundListFileListDf_splittedHere) / maximumNumberOfSpectraPerBatch)
    lapply(X = seq_len(numberOfBatches), FUN = function(batchIdx){
      from = (batchIdx - 1) * maximumNumberOfSpectraPerBatch + 1
      to   = from + maximumNumberOfSpectraPerBatch - 1
      to   = ifelse(test = to > nrow(compoundListFileListDf_splittedHere), yes = nrow(compoundListFileListDf_splittedHere), no = to)
      return(compoundListFileListDf_splittedHere[from:to,])
    })
  }), recursive = FALSE)
  #compoundListFileListDf_splitted <- unlist(lapply())
  
  
  cat(paste("\n", "Splitted ", nrow(compoundListFileListDf), " compounds into ", length(compoundListFileListDf_splitted), " settings", sep = ""))
  
  ## create setting folders and files
  compoundListPathSuffixes <- vector(mode = "character", length = length(compoundListFileListDf_splitted))
  for(partIdx in seq_along(compoundListFileListDf_splitted))
    compoundListPathSuffixes[[partIdx]] <- paste(partIdx, "_", length(compoundListFileListDf_splitted), "__", length(unique(compoundListFileListDf_splitted[[partIdx]]$MspFile)), "_files__", nrow(compoundListFileListDf_splitted[[partIdx]]), "_cpds", sep = "")
  
  folder_parts <- paste(folderRMassBankData, "/", "Part", "_", compoundListPathSuffixes, sep = "")
  
  fileCompoundListFileList_parts <- paste(folder_parts, "/", "CompoundList", "_", compoundListPathSuffixes, ".csv", sep = "")
  fileInfoList_parts             <- paste(folder_parts, "/", "InfoList",     "_", compoundListPathSuffixes, ".csv", sep = "")
  fileSettings_parts             <- paste(folder_parts, "/", "Settings",     "_", compoundListPathSuffixes, ".ini", sep = "")
  
  cat(paste("\n\nCreating ", length(folder_parts), " compound list part folders", sep = ""))
  #partIdx <- 1
  for(partIdx in seq_along(folder_parts)){
    cat(paste("\n", ifelse(test = file.exists(folder_parts[[partIdx]]), yes = "exists: ", no = "is being created :"), folder_parts[[partIdx]], sep = ""))
    cat(paste(ifelse(test = file.exists(folder_parts[[partIdx]]), yes = "exists:", no = "is being created:"), folder_parts[[partIdx]]))
    if(file.exists(folder_parts[[partIdx]]))  next
    
    if(!dir.create(folder_parts[[partIdx]]))  stop(paste("Could not create folder", folder_parts[[partIdx]]))
    
    infoListDf_split <- infoListDf[infoListDf$id %in% compoundListFileListDf_splitted[[partIdx]]$ID, ]
    compoundListFileListDf_splitted[[partIdx]]$"CAS" <- infoListDf_split$"CH$LINK.CAS"[match(x = compoundListFileListDf_splitted[[partIdx]]$"ID", table = infoListDf_split$id)]
    
    cat(paste("\nWriting compound-list-file-list part file", fileCompoundListFileList_parts[[partIdx]]))
    write.csv(x = compoundListFileListDf_splitted[[partIdx]], file = fileCompoundListFileList_parts[[partIdx]])
    
    cat(paste("\nWriting info-list part file", fileInfoList_parts[[partIdx]]))
    write.csv(x = infoListDf_split, file = fileInfoList_parts[[partIdx]])
    
    ## settings
    preprocessSettings(fileCompoundListFileList_part = fileCompoundListFileList_parts[[partIdx]], fileSettings_part = fileSettings_parts[[partIdx]], accessionPrefix = accessionPrefix, ms_type = ms_type)
  }
}

compileSummary_precursor <- function(mb, folder){
  #########################################################################################################################
  ## gather summary of precursors
  
  ## gather values
  attrPrecursors <- lapply(X = mb@spectra, FUN = function(spec){
    attr1 <- attributes(spec)
    attr1 <- attr1[-which(names(attr1) %in% c("children", "parent", "class"))]
    attr2 <- attributes(spec@parent)
    attr2 <- attr2[-which(names(attr2) %in% c(".__classVersion__", "class"))]
    attrPrecursor <- c(attr1, attr2)
    attrPrecursor <- unlist(attrPrecursor)
    attrPrecursor <- attrPrecursor[!duplicated(names(attrPrecursor))]
    if(all(!(c("mz", "polarity", "rt", "mz", "intensity") %in% names(attrPrecursor)))) attrPrecursor[c("polarity", "rt", "mz", "intensity")] <- NA
    return(attrPrecursor)
  })
  
  attrPrecursors <- lapply(X = attrPrecursors, FUN = function(x){x[names(attrPrecursors[[1]])]})
  
  ## box values
  dfAttrPrecursors <- as.data.frame(do.call(rbind, attrPrecursors), stringsAsFactors = FALSE)
  
  ## reorder columns
  dfAttrPrecursors <- cbind(dfAttrPrecursors[, "id", drop=FALSE], dfAttrPrecursors[, -which(colnames(dfAttrPrecursors) %in% c("id"))])
  
  ## write values
  fileName_precursor <- paste("Summary_", basename(folder), "_precursor.tsv", sep = "")
  filePath_precursor <- file.path(getOption("RMassBank")$annotations$entry_prefix, fileName_precursor)
  write.table(x = dfAttrPrecursors, file = filePath_precursor, sep = "\t", row.names = FALSE)
}
compileSummary_children <- function(mb, folder){
  #########################################################################################################################
  ## gather summary of children
  
  ## gather values
  children <- unlist(lapply(X = mb@spectra, FUN = function(spec){
    childList <- as.list(attributes(spec)[["children"]])
    for(childIdx in seq_along(childList)){
      childList[[childIdx]]@fromFile  <- as.integer(spec@id)
      childList[[childIdx]]@scanIndex <- childIdx
    }
    return(childList)
  }))
  attrChildren <- lapply(X = children, FUN = function(child){
    attr3 <- attributes(child)
    attr3 <- attr3[-which(names(attr3) %in% c("satellite", "low", "rawOK", "good", "mzCalc", "formula", "dbe", "formulaCount", "dppm", "dppmBest", "mz", "intensity", "class", ".__classVersion__"))]
    attr3$info <- paste(unlist(attr3$info), collapse = "; ")
    names(attr3)[names(attr3)=="fromFile"]  <- "cpdID"
    names(attr3)[names(attr3)=="scanIndex"] <- "childIdx"
    return(unlist(attr3))
  })
  
  ## box values
  dfAttrChildren   <- as.data.frame(do.call(rbind, attrChildren), stringsAsFactors = FALSE)
  
  ## reorder columns
  dfAttrChildren   <- cbind(dfAttrChildren[, "cpdID", drop=FALSE], dfAttrChildren[, "childIdx", drop=FALSE], dfAttrChildren[, -which(colnames(dfAttrChildren) %in% c("cpdID", "childIdx"))])
  
  ## write values
  fileName_children  <- paste("Summary_", basename(folder), "_spectra.tsv",   sep = "")
  filePath_children  <- file.path(getOption("RMassBank")$annotations$entry_prefix, fileName_children)
  write.table(x = dfAttrChildren,   file = filePath_children,  sep = "\t", row.names = FALSE)
}
compileSummary_peaks <- function(mb, folder){
  #########################################################################################################################
  ## gather summary of peaks
  
  ## gather values
  children <- unlist(lapply(X = mb@spectra, FUN = function(spec){
    childList <- as.list(attributes(spec)[["children"]])
    for(childIdx in seq_along(childList)){
      childList[[childIdx]]@fromFile  <- as.integer(spec@id)
      childList[[childIdx]]@scanIndex <- childIdx
    }
    return(childList)
  }))
  attrPeaks <- lapply(X = children, FUN = function(child){
    attr3 <- attributes(child)
    attr3 <- attr3[which(names(attr3) %in% c("satellite", "low", "rawOK", "good", "mzCalc", "formula", "dbe", "formulaCount", "dppm", "dppmBest", "mz", "intensity"))]
    attr3$"cpdID"    <- child@fromFile
    attr3$"childIdx" <- child@scanIndex
    return(as.data.frame(attr3))
  })
  
  ## box values
  dfAttrPeaks      <- as.data.frame(do.call(rbind, attrPeaks), stringsAsFactors = FALSE)
  
  ## merge with mb@aggregated
  dfAttrPeaks <- cbind(dfAttrPeaks[, c("childIdx", "satellite", "low", "rawOK")], mb@aggregated)
  
  ## reorder columns
  dfAttrPeaks      <- cbind(dfAttrPeaks[, "cpdID", drop=FALSE], dfAttrPeaks[, "childIdx", drop=FALSE], dfAttrPeaks[, -which(colnames(dfAttrPeaks) %in% c("cpdID", "childIdx"))])
  
  ## write values
  fileName_peaks     <- paste("Summary_", basename(folder), "_peaks.tsv",     sep = "")
  filePath_peaks     <- file.path(getOption("RMassBank")$annotations$entry_prefix, fileName_peaks)
  write.table(x = dfAttrPeaks,      file = filePath_peaks,     sep = "\t", row.names = FALSE)
}
compileSummary_records <- function(folder, babel_dir){
  #########################################################################################################################
  ## gather summary of records
  filePath_recData_valid   <- file.path(getOption("RMassBank")$annotations$entry_prefix, "recdata")
  filePath_recData_invalid <- file.path(getOption("RMassBank")$annotations$entry_prefix, "recdata_invalid")
  
  fileName_records_valid   <- paste("Summary_", basename(folder), "_records_valid.tsv",   sep = "")
  filePath_records_valid   <- file.path(getOption("RMassBank")$annotations$entry_prefix, fileName_records_valid)
  fileName_records_invalid <- paste("Summary_", basename(folder), "_records_invalid.tsv", sep = "")
  filePath_records_invalid <- file.path(getOption("RMassBank")$annotations$entry_prefix, fileName_records_invalid)
  
  getInfoFixKey(Directory = filePath_recData_valid,   csvname = filePath_records_valid,   babel_dir = babel_dir, verbose.output = FALSE)
  getInfoFixKey(Directory = filePath_recData_invalid, csvname = filePath_records_invalid, babel_dir = babel_dir, verbose.output = FALSE)
  
  if(file.exists(filePath_records_valid)){
    recordDf_valid   <- read.csv(file = filePath_records_valid,   stringsAsFactors = FALSE, comment.char = "")
    write.table(x = recordDf_valid,   file = filePath_records_valid,   sep = "\t", row.names = FALSE)
  }
  
  if(file.exists(filePath_records_invalid)){
    recordDf_invalid <- read.csv(file = filePath_records_invalid, stringsAsFactors = FALSE, comment.char = "")
    write.table(x = recordDf_invalid, file = filePath_records_invalid, sep = "\t", row.names = FALSE)
  }
}
compileSummary <- function(mb, folder, babel_dir, compilePrecursor=TRUE, compileChildren=TRUE, compilePeaks=TRUE, compileRecord=TRUE){
  if(compilePrecursor) compileSummary_precursor(mb, folder)
  if(compileChildren)  compileSummary_children(mb, folder)
  if(compilePeaks)     compileSummary_peaks(mb, folder)
  if(compileRecord)    compileSummary_records(folder, babel_dir)
}
removePseudoStructures <- function(recordFolder){
  recordFiles <- list.files(path = recordFolder, full.names = TRUE)
  for(recordFile in recordFiles){
    lines <- readLines(recordFile)
    if(any(grepl(x = lines, pattern = "^CH\\$SMILES: (\\[[A-Z][a-z]?[a-z]?\\])+$"))){
      lines[[which(grepl(x = lines, pattern = "^CH\\$SMILES: (\\[[A-Z][a-z]?[a-z]?\\])+$"))]] <- "CH$SMILES: N/A"
      lines[[which(grepl(x = lines, pattern = "^CH\\$IUPAC: "))]] <- "CH$IUPAC: N/A"
      lines <- lines[-which(grepl(x = lines, pattern = "^CH\\$LINK: "))]
      writeLines(text = lines, con = recordFile)
    }
  }
}
recordToMsp <- function(recordFolder, filePath_msp){
  # + NAME: acetyl-CoA
  # - SCANNUMBER: 888
  # + RETENTIONTIME: 17.07485
  # + PRECURSORMZ: 808.1318
  # + PRECURSORTYPE: [M-H]-
  # + IONMODE: Negative
  # + INCHIKEY: ZSLZBFCDCINBPY-ZSJPKINUSA-N
  # + INCHI: InChI=1S/C23H38N7O17P3S/c1-12(31)51-7-6-25-14(32)4-5-26-21(35)18(34)23(2,3)9-44-50(41,42)47-49(39,40)43-8-13-17(46-48(36,37)38)16(33)22(45-13)30-11-29-15-19(24)27-10-28-20(15)30/h10-11,13,16-18,22,33-34H,4-9H2,1-3H3,(H,25,32)(H,26,35)(H,39,40)(H,41,42)(H2,24,27,28)(H2,36,37,38)/t13-,16-,17-,18+,22-/m1/s1
  # + SMILES: CC(=O)SCCNC(=O)CCNC(=O)[C@@H](C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@@H]1[C@H]([C@H]([C@@H](O1)N2C=NC3=C2N=CN=C3N)O)OP(=O)(O)O)O
  # + FORMULA: C23H38N7O17P3S
  # ? INTENSITY: 706897.6
  # - ISOTOPE: M + 0
  # + AUTHORS: Gerd Balcke and Xuemei Lin
  # + LICENSE: ???
  # + COLLISIONENERGY: CES -80, -45, -10 V
  # + INSTRUMENTTYPE: QTOF
  # + INSTRUMENT: Sciex 5600 TripleToF
  # + COMMENT: pure standard
  # - COMMENT2: NRG
  # - COMMENT3: SWATH
  # + Num Peaks: 73
  # + 78.95848    8025
  #  96.96798    1216
  #  107.03656    40
  #  134.04816    6877
  #  139.49448    22
  
  recordFiles <- list.files(path = recordFolder, full.names = TRUE)
  if(length(recordFiles) == 0){
    cat(paste("### Warning ### No msp generated: No record files in folder", recordFolder, "\n"))
    return()
  }
  
  ## iterate MassBank records and convert to msp format
  mspFileLines <- vector(mode = "character")
  
  getRecordEntry <- function(record, tag){
    substring(text = grep(pattern = tag, x = record, value = TRUE, fixed = TRUE), first = nchar(tag) + 1)
  }
  setMspEntry <- function(mspFileLines, tag, value){
    mspFileLines[length(mspFileLines)+1] <- paste(tag, value, sep = ": ")
    return(mspFileLines)
  }
  
  for(recordFile in recordFiles){
    record <- readLines(con = recordFile)
    
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "NAME",            value =          getRecordEntry(record = record, tag = "CH$NAME: ")[[1]])
    #mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "SCANNUMBER",      value = NA)
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "RETENTIONTIME",   value = gsub(x = getRecordEntry(record = record, tag = "AC$CHROMATOGRAPHY: RETENTION_TIME "), pattern = " min$", replacement = ""))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "PRECURSORMZ",     value =          getRecordEntry(record = record, tag = "MS$FOCUSED_ION: PRECURSOR_M/Z "))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "PRECURSORTYPE",   value =          getRecordEntry(record = record, tag = "MS$FOCUSED_ION: PRECURSOR_TYPE "))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "IONMODE",         value =          getRecordEntry(record = record, tag = "AC$MASS_SPECTROMETRY: ION_MODE "))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "INCHIKEY",        value =          getRecordEntry(record = record, tag = "CH$LINK: INCHIKEY "))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "INCHI",           value =          getRecordEntry(record = record, tag = "CH$IUPAC: "))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "SMILES",          value =          getRecordEntry(record = record, tag = "CH$SMILES: "))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "FORMULA",         value =          getRecordEntry(record = record, tag = "CH$FORMULA: "))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "INTENSITY",       value =          getRecordEntry(record = record, tag = "MS$FOCUSED_ION: PRECURSOR_INTENSITY "))
    #mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "ISOTOPE",         value = NA)
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "AUTHORS",         value =          getRecordEntry(record = record, tag = "AUTHORS: "))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "LICENSE",         value =          getRecordEntry(record = record, tag = "LICENSE: "))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "COLLISIONENERGY", value =          getRecordEntry(record = record, tag = "AC$MASS_SPECTROMETRY: COLLISION_ENERGY "))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "INSTRUMENTTYPE",  value =          getRecordEntry(record = record, tag = "AC$INSTRUMENT_TYPE: "))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "INSTRUMENT",      value =          getRecordEntry(record = record, tag = "AC$INSTRUMENT: "))
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "COMMENT",         value =          getRecordEntry(record = record, tag = "COMMENT: CONFIDENCE "))
    #mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "COMMENT2",        value = NA)
    #mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "COMMENT3",        value = NA)
    mspFileLines <- setMspEntry(mspFileLines = mspFileLines, tag = "Num Peaks",       value =          getRecordEntry(record = record, tag = "PK$NUM_PEAK: "))
    
    peaksStart <- which(grepl(x = record, pattern = "^PK\\$PEAK:")) + 1
    peaksEnd   <- which(grepl(x = record, pattern = "^//$"))      - 1
    peaks  <- strsplit(x = trimws(record[peaksStart: peaksEnd]), split = " ")
    peakMz  <- unlist(lapply(X = peaks, FUN = "[", 1))
    peakInt <- unlist(lapply(X = peaks, FUN = "[", 2))
    mspFileLines[(length(mspFileLines)+1):(length(mspFileLines)+0 + length(peaks))] <- paste(peakMz, peakInt, sep = "\t")
    
    mspFileLines[length(mspFileLines)+1] <- ""
  }
  
  ## write msp
  writeLines(text = mspFileLines, con = filePath_msp)
}
aggregateMsp <- function(folder_parts, folderMsp){
  #########################################################################################################################
  ## split by polarity and copy msp to right destination
  #folderMsp_pos <- paste(folderMsp, "pos", sep = "/")
  #folderMsp_neg <- paste(folderMsp, "neg", sep = "/")
  #if(!file.exists(folderMsp_pos)) if(!dir.create(folderMsp_pos))  stop(paste("Could not create folder", folderMsp_pos))
  #if(!file.exists(folderMsp_neg)) if(!dir.create(folderMsp_neg))  stop(paste("Could not create folder", folderMsp_neg))
  
  filePaths_msp_valid <- paste(
    folder_parts, "/",
    getOption("RMassBank")$annotations$entry_prefix, "/",
    "Summary_", basename(folder_parts), "_records_valid.msp", 
    sep = ""
  )
  cat(paste("\nAggregate msp from ", length(folder_parts), " folders (", sum(file.exists(filePaths_msp_valid)), ") msp files", sep = ""))
  
  ## collect
  fileLinesMspPos <- vector(mode = "character")
  fileLinesMspNeg <- vector(mode = "character")
  ## partIdx <- 1
  for(partIdx in seq_along(folder_parts)){
    #setwd(folder_parts[[partIdx]])
    #fileName_msp_valid   <- paste("Summary_", basename(folder_parts[[partIdx]]), "_records_valid.msp",   sep = "")
    #filePath_msp_valid   <- file.path(getOption("RMassBank")$annotations$entry_prefix, fileName_msp_valid)
    filePath_msp_valid   <- filePaths_msp_valid[[partIdx]]
    
    if(!file.exists(filePath_msp_valid))  next
    
    fileLinesMsp_part <- readLines(con = filePath_msp_valid)
    ionMode <- strsplit(x = unique(grep(x = fileLinesMsp_part, pattern = "^IONMODE: ", value = TRUE)), split = ": ")[[1]][[2]]
    
    switch(ionMode, 
           "POSITIVE"={
             fileLinesMspPos <- c(fileLinesMspPos, fileLinesMsp_part)
           },
           "NEGATIVE"={
             fileLinesMspNeg <- c(fileLinesMspNeg, fileLinesMsp_part)
           },
           {  stop(paste("Unknown ionMode '", ionMode, "'"))  }
    )
  }
  
  ## write
  fileNameMspPos <- paste(getOption("RMassBank")$annotations$entry_prefix, "_pos.msp", sep="")
  fileNameMspNeg <- paste(getOption("RMassBank")$annotations$entry_prefix, "_neg.msp", sep="")
  #filePathMspPos <- paste(folderMsp_pos, fileNameMspPos, sep="")
  #filePathMspNeg <- paste(folderMsp_neg, fileNameMspNeg, sep="")
  filePathMspPos <- paste(folderMsp, "/", fileNameMspPos, sep="")
  filePathMspNeg <- paste(folderMsp, "/", fileNameMspNeg, sep="")
  if(length(fileLinesMspPos) > 0)
    writeLines(text = fileLinesMspPos, con = filePathMspPos)
  if(length(fileLinesMspNeg) > 0)
    writeLines(text = fileLinesMspNeg, con = filePathMspNeg)
}
aggregateMassBankRecords <- function(folder_parts, folderMassBank){
  #########################################################################################################################
  ## collect all MassBank records and check for duplicated accessions
  
  cat(paste("\nAggregate MassBank records from", length(folder_parts), "folders"))
  
  ## collect
  library("tools")
  ## partIdx <- 1
  filePaths <- NULL
  for(partIdx in seq_along(folder_parts)){
    setwd(folder_parts[[partIdx]])
    filePath_recData_valid   <- file.path(getOption("RMassBank")$annotations$entry_prefix, "recdata")
    files <- list.files(path = filePath_recData_valid, full.names = TRUE)
    if(length(files) == 0) next
    files <- sapply(X = files, FUN = file_path_as_absolute)
    filePaths <- c(filePaths, files)
  }
  
  ## check
  if(sum(duplicated(basename(filePaths))) > 0)
    stop("Duplicated accessions")
  
  ## copy
  destinationPaths <- paste(folderMassBank, basename(filePaths), sep = "/")
  if(any(file.exists(destinationPaths))) if(!all(file.remove(destinationPaths[file.exists(destinationPaths)]))) stop("Not all record files successfully removed before copying")
  
  if(!all(file.copy(from = filePaths, to = destinationPaths))) stop("Not all record files successfully copied")
  
  ## remove: COMMENT: INTERNAL_ID
  for(massbankRecordFilePath in list.files(path = folderMassBank, full.names = TRUE)){
    lines <- readLines(con = massbankRecordFilePath)
    lines <- lines[!grepl(x = lines, pattern = "^COMMENT: INTERNAL_ID")]
    writeLines(text = lines, con = massbankRecordFilePath)
  }
}
aggregateSummary <- function(folder, folder_parts){
  parentFolder <- folder
  
  cat(paste("\nAggregate summary records from", length(folder_parts), "folders"))
  cat("\nPart X / XX Prec. Child Peaks")
  
  aggregatedPrecursorDfList <- NULL
  aggregatedChildrenDfList <- NULL
  aggregatedPeaksDfList <- NULL
  for(partIdx in seq_along(folder_parts)){
    setwd(folder_parts[[partIdx]])
    folderHere <- folder_parts[[partIdx]]
    fileName_precursor <- paste("Summary_", basename(folderHere), "_precursor.tsv", sep = "")
    filePath_precursor <- file.path(getOption("RMassBank")$annotations$entry_prefix, fileName_precursor)
    fileName_children  <- paste("Summary_", basename(folderHere), "_spectra.tsv",   sep = "")
    filePath_children  <- file.path(getOption("RMassBank")$annotations$entry_prefix, fileName_children)
    fileName_peaks     <- paste("Summary_", basename(folderHere), "_peaks.tsv",     sep = "")
    filePath_peaks     <- file.path(getOption("RMassBank")$annotations$entry_prefix, fileName_peaks)
    
    cat(paste("\nPart", partIdx, "/", length(folder_parts), paste(c("pr", "ch", "pk"), file.exists(c(filePath_precursor, filePath_children, filePath_peaks)), collapse = "; ")))
    
    #aggregatedPrecursorDf  <- rbind(aggregatedPrecursorDf, read.table(file = filePath_precursor, sep = "\t", header = TRUE))
    aggregatedPrecursorDfList[[partIdx]]  <- read.table(file = filePath_precursor, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
    aggregatedPrecursorDfList[[partIdx]] <- aggregatedPrecursorDfList[[partIdx]][, c("id", "found", "complete", "empty", "formula", "mz", "name", "mode", "polarity", "msLevel", "peaksCount", "rt", "acquisitionNum", "tic", "intensity", "centroided", "smoothed")]
    
    if(file.exists(filePath_children))
      #aggregatedChildrenDf <- rbind(aggregatedChildrenDf,  read.table(file = filePath_children,  sep = "\t", header = TRUE))
      aggregatedChildrenDfList[[partIdx]] <- read.table(file = filePath_children,  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
    if(file.exists(filePath_peaks))
      #aggregatedPeaksDf    <- rbind(aggregatedPeaksDf,     read.table(file = filePath_peaks,     sep = "\t", header = TRUE))
      aggregatedPeaksDfList[[partIdx]]               <- read.table(file = filePath_peaks,     sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, comment.char = "")
    
    #if("Petunidin-3-Caffeoyl rutinoside-5-glucoside" %in% precursorDf$name)
    #  print(folder_parts[[partIdx]])
  }
  aggregatedPrecursorDf <- data.table::rbindlist(aggregatedPrecursorDfList)
  aggregatedChildrenDf <- data.table::rbindlist(aggregatedChildrenDfList)
  aggregatedPeaksDf <- data.table::rbindlist(aggregatedPeaksDfList)
  
  
  cat("\nPostprocess summary")
  aggregatedPrecursorDf <- aggregatedPrecursorDf[order(aggregatedPrecursorDf$id   ), ]
  aggregatedChildrenDf  <- aggregatedChildrenDf[ order(aggregatedChildrenDf $cpdID), ]
  aggregatedPeaksDf     <- aggregatedPeaksDf[    order(aggregatedPeaksDf    $cpdID), ]
  
  folderRMassBankData <- paste(parentFolder, "RMassBank", sep = "/")
  
  fileName_summaryPrecursor <- paste("Summary_precursor.tsv", sep = "")
  filePath_summaryPrecursor <- paste(folderRMassBankData, fileName_summaryPrecursor, sep = "/")
  fileName_summaryChildren  <- paste("Summary_spectra.tsv",   sep = "")
  filePath_summaryChildren  <- paste(folderRMassBankData, fileName_summaryChildren, sep = "/")
  fileName_summaryPeaks     <- paste("Summary_peaks.tsv",     sep = "")
  filePath_summaryPeaks     <- paste(folderRMassBankData, fileName_summaryPeaks, sep = "/")
  fileName_summaryPeaksSlim <- paste("Summary_precursors_and_peaks.tsv",     sep = "")
  filePath_summaryPeaksSlim <- paste(folderRMassBankData, fileName_summaryPeaksSlim, sep = "/")
  
  polarities <- substring(first = 1, last = 1, text = aggregatedPrecursorDf$mode)
  polarities <- ifelse(test = polarities == "p", yes = "Positive", no = ifelse(test = polarities == "m", yes = "Negative", no = NA))
  aggregatedPrecursorDf$polarity <- polarities
  colnames(aggregatedPrecursorDf)[[which(colnames(aggregatedPrecursorDf) == "id"       )]] <- "cpdID"
  
  aggregatedPrecursorDf2 <- aggregatedPrecursorDf[, c("cpdID", "name", "mode", "polarity", "rt", "mz", "intensity", "found", "empty")]
  niceColumnNames <- c("cpdID", "precursor_name", "precursor_mode", "precursor_polarity", "precursor_rt", "precursor_mz", "precursor_intensity", "precursor_found", "spectrum_empty")
  colnames(aggregatedPrecursorDf2) <- niceColumnNames
  #colnames(aggregatedPrecursorDf2)[[which(colnames(aggregatedPrecursorDf2) == "intensity")]] <- "precursor_intensity"
  colnames(aggregatedPeaksDf    )[[which(colnames(aggregatedPeaksDf)     == "intensity")]] <- "fragment_intensity"
  aggregatedChildrenDf <- merge(x = aggregatedChildrenDf, y = aggregatedPrecursorDf2, by = "cpdID")
  aggregatedPeaksDf    <- merge(x = aggregatedPeaksDf,    y = aggregatedPrecursorDf2, by = "cpdID")
  
  write.table(x = aggregatedPrecursorDf, file = filePath_summaryPrecursor, sep = "\t", row.names = FALSE)
  write.table(x = aggregatedChildrenDf,  file = filePath_summaryChildren,  sep = "\t", row.names = FALSE)
  write.table(x = aggregatedPeaksDf,     file = filePath_summaryPeaks,     sep = "\t", row.names = FALSE)
  
  #aggregatedPeaksDf <- read.table(file = filePath_summaryPeaks, header = TRUE, sep = "\t", stringsAsFactors = F, check.names = FALSE, comment.char = "")
  
  aggregatedPeaksDf      <- as.data.frame(aggregatedPeaksDf)
  aggregatedPrecursorDf2 <- as.data.frame(aggregatedPrecursorDf2)
  
  
  cat("\nCollect nice summary")
  #aggregatedPeaksDf$"formula" <- as.character(aggregatedPeaksDf$"formula")
  #aggregatedPeaksDf$"reanalyzed.formula" <- as.character(aggregatedPeaksDf$"reanalyzed.formula")
  #aggregatedPeaksDf$"precursor_name" <- as.character(aggregatedPeaksDf$"precursor_name")
  #aggregatedPeaksDf$"precursor_mode" <- as.character(aggregatedPeaksDf$"precursor_mode")
  #
  #aggregatedPrecursorDf2$"precursor_name" <- as.character(aggregatedPrecursorDf2$"precursor_name")
  #aggregatedPrecursorDf2$"precursor_mode" <- as.character(aggregatedPrecursorDf2$"precursor_mode")
  
  ## generate decent peaks list, i.e. one row per input peak
  #aggregatedPeaksDf <- read.table(file = "/mnt/data/IPB/Projects/2017_005_MS-databases/mFam contributions/Jean-Luc Wolfender University of Geneva/RMassBank/Summary_peaks.tsv", sep = "", header = T, check.names = F, comment.char = "", stringsAsFactors = FALSE)
  #cpdID	childIdx	satellite	low	rawOK	mzFound	fragment_intensity	good	mzCalc	formula	dbe	formulaCount	dppm	dppmBest	scan	parentScan	dppmRc	index	noise	reanalyzed.formula	reanalyzed.mzCalc	reanalyzed.dppm	reanalyzed.formulaCount	reanalyzed.dbe	matchedReanalysis	formulaMultiplicity	filterOK	problematicPeak	precursor_name	precursor_mode	precursor_polarity	precursor_rt	precursor_mz	precursor_intensity
  aggregatedPeaksDf_slim <- data.frame(
    "precursor_name"  = character(),
    "precursor_mode"  = character(),
    "precursor_polarity"  = character(),
    "precursor_mz"    = numeric(),
    "precursor_rt"    = numeric(),
    "precursor_intensity" = numeric(),
    "cpdID"     = integer(),
    "precursor_found" = logical(),
    "spectrum_empty"  = logical(),
    #"childIdx"  = integer(),
    "mzFound"   = numeric(),
    "fragment_intensity"  = integer(),
    "filterOK"  = logical(),
    "rawOK"     = logical(),
    "noise"     = logical(),
    "problematicPeak" = logical(),
    "good"      = logical(),
    "satellite" = logical(),
    "low"       = logical(),
    "mzCalc"    = numeric(),
    "formula"   = character(),
    "dbe"       = numeric(),
    "formulaCount"    = integer(),
    "dppm"      = numeric(),
    "formulaMultiplicity" = integer(),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  colNamesHere <- colnames(aggregatedPeaksDf_slim)
  #aggregatedPrecursorDf2 <- data.frame(aggregatedPrecursorDf2)
  #for(cpdID in unique(aggregatedPeaksDf$cpdID)){
  for(cpdID in as.integer(unique(aggregatedPrecursorDf2$cpdID))){
    if(cpdID %in% aggregatedPeaksDf$cpdID){
      aggregatedPeaksDf_part <- aggregatedPeaksDf[aggregatedPeaksDf$cpdID==cpdID, , drop=FALSE]
      for(childIdx in unique(aggregatedPeaksDf_part$childIdx)){
        aggregatedPeaksDf_part2 <- aggregatedPeaksDf_part[aggregatedPeaksDf_part$childIdx==childIdx, , drop=FALSE]
        for(mzFound in unique(aggregatedPeaksDf_part2$mzFound)){
          aggregatedPeaksDf_part3 <- aggregatedPeaksDf_part2[aggregatedPeaksDf_part2$mzFound==mzFound, , drop=FALSE]
          if(nrow(aggregatedPeaksDf_part3) == 1){
            #print("Case 1")
            aggregatedPeaksDf_slim[nrow(aggregatedPeaksDf_slim)+1, ] <- aggregatedPeaksDf_part3[, colNamesHere]
            next
          } 
          if(sum(aggregatedPeaksDf_part3$filterOK) == 1){
            #print("Case 2")
            aggregatedPeaksDf_slim[nrow(aggregatedPeaksDf_slim)+1, ] <- aggregatedPeaksDf_part3[aggregatedPeaksDf_part3$filterOK, colNamesHere]
            next
          }
          #print("Case 3")
          aggregatedPeaksDf_slim[nrow(aggregatedPeaksDf_slim)+1, ] <- aggregatedPeaksDf_part3[which.min(aggregatedPeaksDf_part3$dppm), colNamesHere]
        }
      }
    } else {
      #print("Case 0")
      aggregatedPeaksDf_slim[nrow(aggregatedPeaksDf_slim)+1, ] <- c(
        aggregatedPrecursorDf2[aggregatedPrecursorDf2$cpdID==cpdID, colNamesHere[colNamesHere %in% niceColumnNames]], rep(x = NA, times = ncol(aggregatedPeaksDf_slim) - ncol(aggregatedPrecursorDf2))
      )
    }
  }
  
  colnames(aggregatedPeaksDf_slim) <- c(colNamesHere[colNamesHere %in% niceColumnNames], c(
    "fragment_mzFound",
    "fragment_intensity",
    "fragment_filterOK",
    "fragment_rawOK",
    "fragment_noise",
    "fragment_problematicPeak",
    "fragment_good",
    "fragment_satellite",
    "fragment_low",
    "fragment_mzCalc",
    "fragment_formula",
    "fragment_dbe",
    "fragment_formulaCount",
    "fragment_dppm",
    "fragment_formulaMultiplicity"
  ))
  
  cat("\nPostprocess nice summary")
  ## replace precursor mode by readyble adduct
  adducts <- getAdductInformation("")[, c("mode", "adductString")]
  aggregatedPeaksDf_slim$precursor_mode <- unlist(sapply(X = aggregatedPeaksDf_slim$precursor_mode, FUN = function(precursor_mode){
    if(precursor_mode %in% adducts$mode)
      return(adducts$adductString[adducts$mode==precursor_mode])
    else return(NA)
  }))
  
  ## add meta data
  fileMetaDataCompoundInformationProcessed         <- paste(folderRMassBankData, "MetaDataCompoundInformationProcessed.csv", sep = "/")
  metaDataCompoundsDf <- read.table(file = fileMetaDataCompoundInformationProcessed, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  
  targetColIdx <- which(names(aggregatedPeaksDf_slim) == 'precursor_mz')[[1]]
  aggregatedPeaksDf_slim <- cbind(aggregatedPeaksDf_slim[,1:targetColIdx,drop=F], "Collision energy" = "", aggregatedPeaksDf_slim[,(targetColIdx+1):ncol(aggregatedPeaksDf_slim),drop=F], stringsAsFactors = FALSE)
  for(cpdID in unique(aggregatedPeaksDf_slim$cpdID))
    if(cpdID %in% metaDataCompoundsDf$ID)
      aggregatedPeaksDf_slim[aggregatedPeaksDf_slim$cpdID==cpdID, c("Collision energy")] <- metaDataCompoundsDf[metaDataCompoundsDf$ID==cpdID, c("Collision energy")]
  
  write.table(x = aggregatedPeaksDf_slim,     file = filePath_summaryPeaksSlim,     sep = "\t", row.names = FALSE)
  #aggregatedPeaksDf_slim <- read.table(file = filePath_summaryPeaksSlim,     sep = "\t", header = T)
}

runContributorToMassBankWorkflow <- function(folder, applyIntensityThreshold = FALSE, reprocess = FALSE, fileAnnoTable = NULL){
  folderMsp        <- paste(folder, "converted to msp", sep = "/")
  folderMassBank   <- paste(folder, "converted to MassBank", sep = "/")
  
  ## check folder structure
  folderRMassBankData <- paste(folder, "RMassBank", sep = "/")
  folder_parts <- list.files(path = folderRMassBankData, full.names = TRUE, pattern = "^Part_\\d+_\\d+__\\d+_files__\\d+_cpds$")
  
  if(length(folder_parts) == 0) stop("No part folders")
  
  compoundListPathSuffixes <- gsub(x = basename(path = folder_parts), pattern = "^Part_", replacement = "")
  fileCompoundListFileList_parts <- paste(folder_parts, "/", "CompoundList", "_", compoundListPathSuffixes, ".csv", sep = "")
  fileInfoList_parts             <- paste(folder_parts, "/", "InfoList",     "_", compoundListPathSuffixes, ".csv", sep = "")
  fileSettings_parts             <- paste(folder_parts, "/", "Settings",     "_", compoundListPathSuffixes, ".ini", sep = "")
  
  if(!all(file.exists(fileCompoundListFileList_parts))) stop(paste("Missing file", fileCompoundListFileList_parts))
  if(!all(file.exists(fileInfoList_parts)))             stop(paste("Missing file", fileInfoList_parts))
  if(!all(file.exists(fileSettings_parts)))             stop(paste("Missing file", fileSettings_parts))
  
  #########################################################################################################
  ## run RMassBank
  library("RMassBank")
  
  ## adapt RMassBank settings
  RMassBank.env$verbose.output <- TRUE
  RMassBank.env$export.invalid <- TRUE
  RMassBank.env$export.molfiles <- FALSE
  RMassBank.env$strictMsMsSpectraSelection <- TRUE
  
  #loadRmbSettings(file_or_list = fileSettings_parts[[1]])
  
  #for(sourceFile in list.files(path = "/home/htreutle/Code/Java/RMassBank/R", full.names = TRUE))
  #  source(sourceFile)
  partIdxStart <- ifelse(test = reprocess, yes = 1, no = min(which(sapply(X = folder_parts, FUN = function(folder_part){length(list.dirs(path = folder_part, recursive = FALSE))==0}))))
  # folder_parts <- "/mnt/data/IPB/Projects/2017_005_MS-databases/mFam contributions/Ulschan Bathe di and sesquiterpenes/RMassBank/Part_1_1__4_files"
  # partIdx <- 1
  # partIdx <- partIdxStart
  for(partIdx in partIdxStart:length(folder_parts)){
    #for(partIdx in seq_along(folder_parts)){
    cat("\n################################################################################################")
    cat(paste("\n### Processing part", partIdx, "/", length(folder_parts), folder_parts[[partIdx]]))
    setwd(folder_parts[[partIdx]])
    
    compoundListFileListDf_split <- read.csv(file = fileCompoundListFileList_parts[[partIdx]], stringsAsFactors = F)
    
    listSettingsLines   <- readLines(con = fileSettings_parts[[partIdx]])
    listSettingsNames   <- unlist(lapply(X = strsplit(x = listSettingsLines, split = ": "), FUN = function(x){x[[1]]}))
    listSettings        <- unlist(lapply(X = strsplit(x = listSettingsLines, split = ": "), FUN = function(x){ifelse(test = length(x)> 1, yes = x[[2]], no = "")}))
    names(listSettings) <- listSettingsNames
    
    #switch(compoundListFileListDf_split$Ionization.mode[[1]], 
    #       "Positive"={  ionMode <- "pH" },
    #       "Negative"={  ionMode <- "mH" },
    #       {  stop(paste("Unknown Ionization mode '", compoundListFileListDf_split$Ionization.mode[[1]], "'"))  }
    #)
    ## "pH", "pNa", "pM", "pNH4", "mH", "mM", "mFA" for different ions ([M+H]+, [M+Na]+, [M]+, [M+NH4]+, [M-H]-, [M]-, [M+FA]-).
    ionMode <- adductToIonMode(compoundListFileListDf_split$Adduct[[1]])
    runId <- listSettings[["    entry_prefix"]]
    fileMspData <- unique(compoundListFileListDf_split$Files)
    
    cat(paste("\nProcessing spectra ", ionMode, ": ", nrow(compoundListFileListDf_split), " compounds\n", sep = ""))
    
    
    ## settings and compound list
    loadRmbSettings(file_or_list = fileSettings_parts[[partIdx]])
    #options("RMassBank")
    ## do not fix the multiplicityFilter
    #settings <- getOption("RMassBank")
    #settings$multiplicityFilter <- NULL
    #options(RMassBank = settings)
    
    
    #if(length(list.files(file.path(getOption("RMassBank")$annotations$entry_prefix, "recdata"))) > 0)      {print("done");next;}
    #else                                                                                                   {print("TODO"; break;}
    
    
    #loadList(path = "/home/htreutle/Code/Java/IPB_library/MM48/compoundlist/IPBMM48.csv")
    loadList(path = fileCompoundListFileList_parts[[partIdx]])
    
    msmsList <- newMsmsWorkspace()
    msmsList <- msmsRead(w = msmsList, filetable = fileCompoundListFileList_parts[[partIdx]], files = fileMspData, readMethod = "msp", mode = ionMode)
    ## --> leMsmsRow.R --> findMsMsHRperMsp.direct : compoundListFileListDf_split[, c("ID", "Files")]
    #mode <- ionMode
    #fileName <- "/mnt/data/IPB/Projects/2017_005_MS-databases/mFam contributions/Univ Athens Maria Halabalaki 2nd/raw data/exported as raw msp/neg/loganin_neg.msp"
    #cpdIDs <- 19
    
    
    numberOfPeaksThere <- sum(unlist(lapply(X = msmsList@spectra, FUN = function(spec){ sum(unlist(lapply(X = spec@children, FUN = function(child){ child@peaksCount }))) })))
    if(numberOfPeaksThere == 0){
      cat(paste("\n### Warning ### No spectra found or spectra contain no peaks."))
      filePath_recData_prefix   <- file.path(getOption("RMassBank")$annotations$entry_prefix)
      if(!file.exists(filePath_recData_prefix)) if(!dir.create(filePath_recData_prefix,recursive=TRUE))  stop(paste("Could not create folder", filePath_recData_prefix))
      compileSummary(mb = msmsList, folder = folder_parts[[partIdx]], babel_dir = NULL, compilePrecursor=TRUE, compileChildren=FALSE, compilePeaks=FALSE, compileRecord=FALSE)
      next
    }
    
    ## rt: min --> secs
    for(idx1 in seq_along(msmsList@spectra@listData)){
      msmsList@spectra@listData[[idx1]]@parent@rt <- msmsList@spectra@listData[[idx1]]@parent@rt * 60
      for(idx2 in seq_along(msmsList@spectra@listData[[idx1]]@children))
        msmsList@spectra@listData[[idx1]]@children[[idx2]]@rt <- msmsList@spectra@listData[[idx1]]@children[[idx2]]@rt * 60
    }
    
    #msmsList <- msmsWorkflow(w = msmsList, mode = ionMode, steps = 2:3, archivename = runId)
    ##peaksMatched(msmsList)$formulaCount
    #msmsList <- msmsWorkflow(w = msmsList, mode = ionMode, steps = 4:8, archivename = runId)
    #msmsList <- msmsWorkflow(w = msmsList, mode = ionMode, steps = 2:8, archivename = runId, newRecalibration = FALSE)
    #msmsList <- msmsWorkflow(w = msmsList, mode = ionMode, steps = 2:7, archivename = runId)
    msmsList <- msmsWorkflow(w = msmsList, mode = ionMode, steps = c(2,3,6,7,8), archivename = runId, newRecalibration = FALSE)
    
    if(all(!msmsList@aggregated$filterOK)){
      cat(paste("\n### Warning ### All peaks have been filtered."))
      filePath_recData_prefix   <- file.path(getOption("RMassBank")$annotations$entry_prefix)
      if(!file.exists(filePath_recData_prefix)) if(!dir.create(filePath_recData_prefix,recursive=TRUE))  stop(paste("Could not create folder", filePath_recData_prefix))
      compileSummary(mb = msmsList, folder = folder_parts[[partIdx]], babel_dir = NULL, compilePrecursor=TRUE, compileChildren=TRUE, compilePeaks=TRUE, compileRecord=FALSE)
      next
    }
    
    mb <- newMbWorkspace(w = msmsList)
    
    newInfoListFile    <- paste(folder_parts[[partIdx]], "/infolist.csv",        sep = "")
    mergedInfoListFile <- paste(folder_parts[[partIdx]], "/infolist_merged.csv", sep = "")
    
    if(!file.exists(newInfoListFile))
      mb <- mbWorkflow(mb = mb, steps = 1:2, gatherData = "online")
    
    if(!file.exists(mergedInfoListFile)){
      ## merge both infoLists
      infoListNew     <- read.csv(file = newInfoListFile,               check.names = FALSE, stringsAsFactors = FALSE)
      infoListPresent <- read.csv(file = fileInfoList_parts[[partIdx]], check.names = FALSE, stringsAsFactors = FALSE)
      infoListMerged  <- infoListPresent
      
      #cleanedNames <- data.frame("CH$NAME1" = infoListNew$`CH$NAME1`, "CH$NAME2" = infoListNew$`CH$NAME2`, "CH$NAME3" = infoListNew$`CH$NAME3`)
      cleanedNames <- list("CH$NAME1" = infoListNew$`CH$NAME1`, "CH$NAME2" = infoListNew$`CH$NAME2`, "CH$NAME3" = infoListNew$`CH$NAME3`)
      cleanedNames <- lapply(X = cleanedNames, FUN = removeItalicsTagsFromNames)
      cleanedNames <- lapply(X = cleanedNames, FUN = function(cleanedNamesHere){
        if(all(is.na(cleanedNamesHere))) return(NA)
        unlist(lapply(X = strsplit(x = cleanedNamesHere, split = "; "), FUN = function(tokens){
        ifelse(test = length(tokens)==0, yes = "", no = head(x = tokens, n=1))
      }
      ))})
      cleanedNames <- lapply(X = cleanedNames, FUN = function(cleanedNamesHere){gsub(x = cleanedNamesHere, pattern = ", [<>=]*\\d+%$", replacement = "")})
      cleanedNames <- lapply(X = cleanedNames, FUN = trimws)
      cleanedNames <- lapply(X = cleanedNames, FUN = function(cleanedNamesHere){gsub(x = cleanedNamesHere, pattern = "’", replacement = "'")})
      cleanedNames <- lapply(X = cleanedNames, FUN = function(cleanedNamesHere){gsub(x = cleanedNamesHere, pattern = "→", replacement = "->")})
      
      infoListNew$`CH$NAME1` <- cleanedNames[[1]]
      infoListNew$`CH$NAME2` <- cleanedNames[[2]]
      infoListNew$`CH$NAME3` <- cleanedNames[[3]]
      
      infoListNew <- infoListNew[order(infoListNew$id),]
      
      ## remove compounds which have not been found
      infoListMerged <- infoListMerged[infoListMerged$id %in% infoListNew$id, ]
      
      columnNames <- c("dbcas",	"dbname",	"dataused",	#"COMMENT.CONFIDENCE",	"COMMENT.ID",
                       "CH$NAME1", "CH$NAME2", "CH$NAME3", #CH$COMPOUND_CLASS
                       "CH$FORMULA", "CH$EXACT_MASS", #"CH$SMILES", "CH$IUPAC",
                       "CH$LINK.CAS", "CH$LINK.CHEBI", "CH$LINK.HMDB", "CH$LINK.KEGG", "CH$LINK.LIPIDMAPS", "CH$LINK.PUBCHEM", #"CH$LINK.INCHIKEY", 
                       "CH$LINK.CHEMSPIDER")
      for(columnName in columnNames)
        infoListMerged[is.na(infoListMerged[, columnName]), columnName] <- infoListNew[is.na(infoListMerged[, columnName]), columnName]
      infoListMerged$"CH$SMILES"        <- infoListNew$"CH$SMILES"
      infoListMerged$"CH$IUPAC"         <- infoListNew$"CH$IUPAC"
      infoListMerged$"CH$LINK.INCHIKEY" <- infoListNew$"CH$LINK.INCHIKEY"
      infoListMerged$"COMMENT.ID"       <- infoListNew$"COMMENT.ID"
      infoListMerged[is.na(infoListMerged)] <- ""
      
      ##  add SP$SAMPLE
      compoundListFileListDf_split_part <- compoundListFileListDf_split[compoundListFileListDf_split$ID %in% infoListMerged$id, ]
      if(nrow(infoListMerged) != nrow(compoundListFileListDf_split_part))  stop("Cannot flll in SP$SAMPLE")
      infoListMerged[["SP.SAMPLE"]] <- compoundListFileListDf_split_part$Sample[match(x = compoundListFileListDf_split_part$ID, table = infoListMerged$COMMENT.ID)]
      #infoListMerged[["SP.SAMPLE"]] <- compoundListFileListDf_split_part$Sample
      
      write.csv(x = infoListMerged, file = mergedInfoListFile)
    }
    
    mb <- loadInfolist(mb = mb, fileName = mergedInfoListFile)
    #mb <- mbWorkflow(mb = mb, steps = 3)
    mb <- mbWorkflow(mb = mb, steps = c(3,4,5,7))
    #mb <- mbWorkflow(mb = mb, steps = 3:8)
    
    babel_dir <- listSettings["babeldir"]
    compileSummary(mb = mb, folder = folder_parts[[partIdx]], babel_dir=babel_dir)
    
    #########################################################################################################################
    ## gather and export msp
    filePath_recData_valid   <- file.path(getOption("RMassBank")$annotations$entry_prefix, "recdata")
    filePath_recData_invalid <- file.path(getOption("RMassBank")$annotations$entry_prefix, "recdata_invalid")
    
    fileName_msp_valid   <- paste("Summary_", basename(folder_parts[[partIdx]]), "_records_valid.msp",   sep = "")
    filePath_msp_valid   <- file.path(getOption("RMassBank")$annotations$entry_prefix, fileName_msp_valid)
    fileName_msp_invalid <- paste("Summary_", basename(folder_parts[[partIdx]]), "_records_invalid.msp", sep = "")
    filePath_msp_invalid <- file.path(getOption("RMassBank")$annotations$entry_prefix, fileName_msp_invalid)
    
    removePseudoStructures(recordFolder = filePath_recData_valid  )
    removePseudoStructures(recordFolder = filePath_recData_invalid)
    correctRecords(folder = filePath_recData_valid,   applyIntensityThreshold = applyIntensityThreshold)
    correctRecords(folder = filePath_recData_invalid, applyIntensityThreshold = applyIntensityThreshold)
    recordToMsp(recordFolder = filePath_recData_valid,   filePath_msp = filePath_msp_valid)
    recordToMsp(recordFolder = filePath_recData_invalid, filePath_msp = filePath_msp_invalid)
  }
  
  aggregateMsp(folder_parts, folderMsp)
  aggregateMassBankRecords(folder_parts, folderMassBank)
  aggregateSummary(folder, folder_parts)
  createSunBurstPlot(folder, fileAnnoTable)
  resultsWithError <- validateRecords(folder)
  if(length(resultsWithError) > 0){
    print(paste("Error in", length(resultsWithError), "files"))
    print(paste(basename(names(resultsWithError))))
    #print(resultsWithError)
    print(lapply(X = resultsWithError, FUN = function(lines){lines[grepl(x = lines, pattern = "\\[main\\] ERROR")]}))
  }
}

cleanContribtorDirectoryForReprocessing <- function(folder, onlyUpdatePeaks = FALSE){
  folderMsp           <- paste(folder, "converted to msp", sep = "/")
  folderMassBank      <- paste(folder, "converted to MassBank", sep = "/")
  folderRMassBankData <- paste(folder, "RMassBank", sep = "/")
  
  filesMsp           <- list.files(path = folderMsp,           full.names = TRUE)
  filesMassBank      <- list.files(path = folderMassBank,      full.names = TRUE)
  filesRMassBankData <- list.files(path = folderRMassBankData, full.names = TRUE)
  filesRMassBankDataSubFolder <- grepl(x = filesRMassBankData, pattern = "Part_\\d+_\\d+__\\d+_files__\\d+_cpds$")
  filesRMassBankDataSummary   <- grepl(x = filesRMassBankData, pattern = "Summary_.+\\.tsv$")
  filesRMassBankDataMetaData  <- !(filesRMassBankDataSubFolder | filesRMassBankDataSummary)
  
  if(!onlyUpdatePeaks) if(unlink(filesRMassBankData[filesRMassBankDataMetaData]))           stop(paste("Could not remove files from folder", folderRMassBankData))
  if(unlink(filesRMassBankData[filesRMassBankDataSummary]))           stop(paste("Could not remove files from folder", folderRMassBankData))
  if(!onlyUpdatePeaks){
    if(unlink(filesRMassBankData[filesRMassBankDataSubFolder], recursive = T))           stop(paste("Could not remove files from folder", folderRMassBankData))
  } else {
    sapply(X = filesRMassBankData[filesRMassBankDataSubFolder], FUN = function(subFolder){
      filesInSubFolder <- list.files(path = subFolder, full.names = TRUE)
      
      isMetaDataRegEx <- paste(paste(c(
        "(CompoundList_\\d+_\\d+__\\d+_files__\\d+_cpds.csv)",
        "(InfoList_\\d+_\\d+__\\d+_files__\\d+_cpds\\.csv)",
        "(infolist\\.csv)",
        "(infolist_merged\\.csv)",
        "(Settings_\\d+_\\d+__\\d+_files__\\d+_cpds\\.ini)"), 
        collapse="|"), "$", sep="")
      
      isMetaData <- grepl(x = filesInSubFolder, pattern = isMetaDataRegEx)
      cat(paste("\nRemoving", sum(!isMetaData), "files from folder", subFolder))
      if(unlink(x = filesInSubFolder[!isMetaData], recursive = T))           stop(paste("Could not remove files from folder", subFolder))
    })
  }
  
  cat(paste("\nRemoving", length(filesMsp), "files from folder", folderMsp))
  if(unlink(filesMsp))           stop(paste("Could not remove files from folder", folderMsp))
  cat(paste("\nRemoving", length(filesMassBank), "files from folder", folderMassBank))
  if(unlink(filesMassBank))      stop(paste("Could not remove files from folder", folderMassBank))
  
  #cat(paste("\nRemoving", length(filesRMassBankData), "files/folders from folder", folderRMassBankData))
  #if(unlink(filesRMassBankData, recursive = TRUE)) stop(paste("Could not remove files from folder", folderRMassBankData))
  cat("\nCleaning ready")
}
removeEmptySpectra <- function(folder){
  folderRMassBankData <- paste(folder, "RMassBank", sep = "/")
  
  fileCompoundListFileList <- paste(folderRMassBankData, "CompoundList_FileList.csv", sep = "/")
  compoundListFileListDf <- read.table(file = fileCompoundListFileList, header = TRUE, sep = "\t", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  existingFiles <- compoundListFileListDf$Files
  
  cat(paste("\nRemoving empty spectra", length(existingFiles), "files"))
  
  originalDirectories <- paste(unique(dirname(existingFiles)), "original", sep = "/")
  sapply(X = originalDirectories, FUN = function(origDir){ if(!file.exists(origDir)) if(!dir.create(origDir)) stop(paste("Could not create directory", origDir))})
  
  newFiles <- paste(dirname(existingFiles), "/", "original", "/", basename(existingFiles), "_", sep = "")
  for(idx in seq_along(existingFiles)){
    if(file.exists(newFiles[[idx]])) next
    cat(paste("\nRemoving empty spectra from file", idx, "/", length(existingFiles), basename(existingFiles[[idx]])))
    if(!file.copy(from = existingFiles[[idx]], to = newFiles[[idx]])) stop(paste("Failed to copy file", existingFiles[[idx]], "to", newFiles[[idx]]))
    fileLines <- readLines(con = existingFiles[[idx]])
    ends <- which(grepl(x = fileLines, pattern = "^Num Peaks: 0$"))
    starts <- which(grepl(x = fileLines, pattern = "^NAME"))
    starts <- sapply(X = ends, FUN = function(end){
      max(starts[starts < end])
    })
    ends <- ends + 1
    indeces <- unlist(sapply(X = seq_along(starts), FUN = function(idx2){seq(from = starts[[idx2]], to = ends[[idx2]])}))
    fileLines <- fileLines[-indeces]
    writeLines(text = fileLines, con = existingFiles[[idx]])
  }
}
correctRecords <- function(folder, applyIntensityThreshold = FALSE){
  options(warn = 1)
  
  recordFolder <- paste(folder, "converted to MassBank", sep = "/")
  if(!file.exists(recordFolder)) recordFolder <- folder
  filePaths <- list.files(path = recordFolder, pattern = "[A-Z0-9]{8,8}\\.txt$", full.names=TRUE)
  print(paste("Correcting", length(filePaths), "files in folder", folder))
  for(filePath in filePaths){
    lines <- readLines(con = filePath)
    lines2 <- lines
    
    ## switch two lines
    if(FALSE){
      if(!any(grepl(x = lines2, pattern = "^SP\\$SAMPLE")))
        next
      instrumentIdx <- which(grepl(x = lines2, pattern = "^AC\\$INSTRUMENT: "))
      spSampleIdx   <- which(grepl(x = lines2, pattern = "^SP\\$SAMPLE"))
      
      if(spSampleIdx < instrumentIdx) next
      
      lines2 <- c(
        lines2[1:(instrumentIdx-1)],
        lines2[[spSampleIdx]],
        lines2[instrumentIdx:(spSampleIdx-1)],
        lines2[(spSampleIdx + 1):length(lines2)]
      )
      if(!all(all(lines2 %in% lines), all(lines %in% lines2), length(lines) == length(lines2))) stop("oh oh")
    }
    ## switch two lines
    if(TRUE){
      acms_ionization <- which(grepl(x = lines2, pattern = "^AC\\$MASS_SPECTROMETRY: IONIZATION"))
      acms_ionmode    <- which(grepl(x = lines2, pattern = "^AC\\$MASS_SPECTROMETRY: ION_MODE"))
      if(acms_ionmode > acms_ionization){
        if(acms_ionmode - acms_ionization != 1) stop("hae?")
        print(paste("IONIZATION vs ION_MODE", basename(filePath), filePath, sep = "\t"))
        lines2 <- c(
          lines2[1:(acms_ionization-1)],
          lines2[[acms_ionmode]],
          lines2[[acms_ionization]],
          lines2[(acms_ionmode + 1):length(lines2)]
        )
        if(!all(all(lines2 %in% lines), all(lines %in% lines2), length(lines) == length(lines2))) stop("oh oh")
      }
    }
    ## correct lines
    if(TRUE){
      lines2 <- gsub(x = lines2, pattern = "LC-ESI-qTOF", replacement = "LC-ESI-QTOF")
      lines2 <- gsub(x = lines2, pattern = "^\\^PK\\$NUM_PEAK:", replacement = "PK$NUM_PEAK:")
      lines2 <- gsub(x = lines2, pattern = "^PK\\\\\\$NUM_PEAK:", replacement = "PK$NUM_PEAK:")
    }
    ## duplicated CH$NAME
    if(TRUE){
      indeces <- which(grepl(x = lines2, pattern = "^CH\\$NAME"))
      if(any(duplicated(lines2[indeces]))){
        lines2 <- lines2[-indeces[duplicated(lines2[indeces])]]
        print(paste("CH$NAME duplicated", basename(filePath), filePath, sep = "\t"))
      }
    }
    ## spaces in CH$NAME and RECORD_TITLE
    if(TRUE){
      indeces <- which(grepl(x = lines2, pattern = "^CH\\$NAME"))
      names <- gsub(x = lines2[indeces], pattern = "^CH\\$NAME: ", replacement = "")
      if(!all(names == trimws(names))){
        names <- trimws(names)
        lines2[indeces] <- paste("CH$NAME: ", names, sep = "")
        
        index <- which(grepl(x = lines2, pattern = "^RECORD_TITLE"))
        lines2[index] <- gsub(x = lines2[index], pattern = "  ", replacement = " ")
        lines2[index] <- gsub(x = lines2[index], pattern = " ;", replacement = ";")
        
        print(paste("CH$NAME trimws", basename(filePath), filePath, sep = "\t"))
      }
    }
    ## correct italics and some chars in CH$NAME and TITLE
    if(TRUE){
      lines3 <- lines2
      indeces <- c(which(grepl(x = lines2, pattern = "^CH\\$NAME")), which(grepl(x = lines2, pattern = "^RECORD_TITLE")))
      for(index in indeces){
        lines2[[index]] <- removeItalicsTagsFromName(name = lines2[[index]])
        lines2[[index]] <- gsub(x = lines2[[index]], pattern = "’", replacement = "'")
        lines2[[index]] <- gsub(x = lines2[[index]], pattern = "→", replacement = "->")
      }
      if(!all(lines3==lines2))
        print(paste("CH$NAME/TITLE italics and characters ’→", basename(filePath), filePath, sep = "\t"))
    }
    ## remove character :
    if(FALSE){
      if(!any(grepl(x = lines2, pattern = "^AC\\$CHROMATOGRAPHY: COLUMN_TEMPERATURE:")))
        next
      idx <- which(grepl(x = lines2, pattern = "^AC\\$CHROMATOGRAPHY: COLUMN_TEMPERATURE:"))
      lines2[[idx]] <- gsub(x = lines2[[idx]], pattern = "^AC\\$CHROMATOGRAPHY: COLUMN_TEMPERATURE:", replacement = "AC$CHROMATOGRAPHY: COLUMN_TEMPERATURE")
    }
    ## duplicated peak list
    if(TRUE){
      isPeakLine <- grepl(x = lines2, pattern = "^  \\d+.*")
      peakLines <- lines2[isPeakLine]
      
      PK_NUM_PEAK_idx <- which(grepl(x = lines2, pattern = "^PK\\$NUM_PEAK"))
      PK_ANNOTATION_idx <- which(grepl(x = lines2, pattern = "^PK\\$ANNOTATION"))
      if(length(PK_ANNOTATION_idx) == 0){
        peakMzs <- as.numeric(unlist(lapply(X = strsplit(x = trimws(peakLines), split = " "), FUN = "[", 1)))
        if(any(duplicated(peakMzs))){
          peakInts <- as.numeric(unlist(lapply(X = strsplit(x = trimws(peakLines), split = " "), FUN = "[", 2)))
          peakIndecesToRemove <- NULL
          for(idx in which(duplicated(peakMzs))){
            indeces <- which(peakMzs==peakMzs[[idx]])
            peakIndecesToRemove <- c(peakIndecesToRemove, indeces[[which.min(peakInts[indeces])]])
          }
          peakIndecesToRemove <- which(isPeakLine)[peakIndecesToRemove]
          lines2 <- lines2[-peakIndecesToRemove]
          print(paste("Peak duplicates", basename(filePath), filePath, sep = "\t"))
        }
      } else {
        warning(paste("### Peak duplicates: anno not implemented yet", paste(PK_ANNOTATION_idx, collapse = "; "), filePath))
      }
    }
    ## sort peak list
    if(TRUE){
      isPeakLine <- grepl(x = lines2, pattern = "^  \\d+.*")
      peakLines <- lines2[isPeakLine]
      
      PK_NUM_PEAK_idx <- which(grepl(x = lines2, pattern = "^PK\\$NUM_PEAK"))
      PK_ANNOTATION_idx <- which(grepl(x = lines2, pattern = "^PK\\$ANNOTATION"))
      if(length(PK_ANNOTATION_idx) == 0){
        peakMzs <- as.numeric(unlist(lapply(X = strsplit(x = trimws(peakLines), split = " "), FUN = "[", 1)))
        order <- order(peakMzs)
        if(!all(order == seq_along(peakLines))){
          peakLines <- peakLines[order]
          lines2[isPeakLine] <- peakLines
          print(paste("Peak order", basename(filePath), filePath, sep = "\t"))
        }
      } else {
        warning(paste("### Peak order: anno not implemented yet", paste(PK_ANNOTATION_idx, collapse = "; "), filePath))
      }
    }
    ## correct num peaks
    if(TRUE){
      isPeakLine <- grepl(x = lines2, pattern = "^  \\d+.*")
      peakLines <- lines2[isPeakLine]
      
      PK_NUM_PEAK_idx <- which(grepl(x = lines2, pattern = "^PK\\$NUM_PEAK"))
      PK_ANNOTATION_idx <- which(grepl(x = lines2, pattern = "^PK\\$ANNOTATION"))
      if(length(PK_ANNOTATION_idx) == 0){
        peakMzs <- as.numeric(unlist(lapply(X = strsplit(x = trimws(peakLines), split = " "), FUN = "[", 1)))
        numPeaks <- as.numeric(gsub(x = lines2[[PK_NUM_PEAK_idx]], pattern = "^PK\\$NUM_PEAK: ", replacement = ""))
        if(numPeaks != length(peakMzs)){
          lines2[[PK_NUM_PEAK_idx]] <- paste("PK$NUM_PEAK: ", length(peakMzs), sep = "")
          print(paste("Num Peak", basename(filePath), filePath, sep = "\t"))
        }
      } else {
        warning(paste("### Num Peak: anno not implemented yet", paste(PK_ANNOTATION_idx, collapse = "; "), filePath))
      }
    }
    ## correct ion charge
    if(TRUE){
      if(any(grepl(x = lines2, pattern = "^CH\\$FORMULA: [A-Z0-9]+[\\+\\-]$"))){
        idx <- which(grepl(x = lines2, pattern = "^CH\\$FORMULA: [A-Z0-9]+[\\+\\-]$"))
        formula <- substr(x = lines2[[idx]], start = nchar("CH$FORMULA: ") + 1, stop = nchar(lines2[[idx]]))
        formula <- paste("[", substr(x = formula, start = 1, stop = nchar(formula) - 1), "]", substr(x = formula, start = nchar(formula), stop = nchar(formula)), sep = "")
        lines2[[idx]] <- paste("CH$FORMULA: ", formula, sep = "")
        print(paste("Charged formula", basename(filePath), filePath, sep = "\t"))
      }
    }
    ## correct PK$NUM_PEAK
    if(FALSE){
      PK_ANNOTATION_idx <- which(grepl(x = lines2, pattern = "^PK\\$ANNOTATION"))
      NUM_PEAK_idx <- which(grepl(x = lines2, pattern = "^PK\\$NUM_PEAK"))
      
      isPeakishLine <- grepl(x = lines2, pattern = "^  \\d+.*")
      isPeakLine <- isPeakishLine & c(rep(x = FALSE, times = NUM_PEAK_idx), rep(x = TRUE , times = length(lines2) - NUM_PEAK_idx))
      isAnnoLine <- isPeakishLine & c(rep(x = TRUE , times = NUM_PEAK_idx), rep(x = FALSE, times = length(lines2) - NUM_PEAK_idx))
      peakLines <- lines2[isPeakLine]
      annoLines <- lines2[isAnnoLine]
      
      lines2[[NUM_PEAK_idx]] <- paste("PK$NUM_PEAK: ", length(peakLines), sep = "")
    }
    ## intensity threshold
    if(applyIntensityThreshold){
      intensityProportionThreshold <- 0.005
      
      PK_ANNOTATION_idx <- which(grepl(x = lines2, pattern = "^PK\\$ANNOTATION"))
      NUM_PEAK_idx <- which(grepl(x = lines2, pattern = "^PK\\$NUM_PEAK"))
      
      isPeakishLine <- grepl(x = lines2, pattern = "^  \\d+.*")
      isPeakLine <- isPeakishLine & c(rep(x = FALSE, times = NUM_PEAK_idx), rep(x = TRUE , times = length(lines2) - NUM_PEAK_idx))
      isAnnoLine <- isPeakishLine & c(rep(x = TRUE , times = NUM_PEAK_idx), rep(x = FALSE, times = length(lines2) - NUM_PEAK_idx))
      peakLines <- lines2[isPeakLine]
      annoLines <- lines2[isAnnoLine]
      
      peakInts <- as.numeric(unlist(lapply(X = strsplit(x = trimws(peakLines), split = " "), FUN = "[", 2)))
      removedPeaks <- (peakInts < intensityProportionThreshold * max(peakInts)) | (peakInts == 1)
      #print(max(peakInts))
      
      if(TRUE & sum(removedPeaks) > 0){
        peakLines <- peakLines[!removedPeaks]
        annoLines <- annoLines[!removedPeaks]
        #numPeaks <- as.numeric(gsub(x = lines2[[NUM_PEAK_idx]], pattern = "^PK\\$NUM_PEAK: ", replacement = ""))
        
        lines2[[NUM_PEAK_idx]] <- paste("PK$NUM_PEAK: ", length(peakLines), sep = "")
        lines2 <- c(
          lines2[seq(from=1, to=min(which(isAnnoLine)) - 1)],
          annoLines,
          lines2[seq(from=max(which(isAnnoLine)) + 1, to=min(which(isPeakLine)) - 1)],
          peakLines,
          lines2[seq(from=max(which(isPeakLine)) + 1, to=length(lines2))]
        )
        
        print(paste("IntThreshold", sum(isPeakLine), "-->", length(peakLines), basename(filePath), filePath, sep = "\t"))
      }
    }
    ## WeizMass library
    if(FALSE){
      lines2[lines2 == "AC$INSTRUMENT: TODO"]                   <- "AC$INSTRUMENT: HDMS Synapt, Waters"
      lines2[lines2 == "AC$INSTRUMENT_TYPE: TODO"]              <- "AC$INSTRUMENT_TYPE: LC-ESI-QTOF"
      lines2[lines2 == "AC$MASS_SPECTROMETRY: IONIZATION TODO"] <- "AC$MASS_SPECTROMETRY: IONIZATION ESI"
      lines2[lines2 == "AC$CHROMATOGRAPHY: COLUMN_NAME TODO"]   <- "AC$CHROMATOGRAPHY: COLUMN_NAME Waters Acquity UPLC system"
      lines2[[2]]   <- gsub(x = lines2[[2]], pattern = "TODO; MS2", replacement = "LC-ESI-QTOF; MS2")
      print(paste("WeizMass", basename(filePath), filePath, sep = "\t"))
    }
    ## correct splash
    if(TRUE){
      PK_NUM_PEAK_idx <- which(grepl(x = lines2, pattern = "^PK\\$NUM_PEAK"))
      PK_ANNOTATION_idx <- which(grepl(x = lines2, pattern = "^PK\\$ANNOTATION"))
      PK_SPLASH_idx <- which(grepl(x = lines2, pattern = "^PK\\$SPLASH"))
      
      splashThere <- gsub(x = lines2[[PK_SPLASH_idx]], pattern = "^PK\\$SPLASH: ", replacement = "")
      
      isPeakishLine <- grepl(x = lines2, pattern = "^  \\d+.*")
      isPeakLine <- isPeakishLine & c(rep(x = FALSE, times = PK_NUM_PEAK_idx), rep(x = TRUE , times = length(lines2) - PK_NUM_PEAK_idx))
      isAnnoLine <- isPeakishLine & c(rep(x = TRUE , times = PK_NUM_PEAK_idx), rep(x = FALSE, times = length(lines2) - PK_NUM_PEAK_idx))
      peakLines <- lines2[isPeakLine]
      annoLines <- lines2[isAnnoLine]
      
      peakMzs <- as.numeric(unlist(lapply(X = strsplit(x = trimws(peakLines), split = " "), FUN = "[", 1)))
      peakInts <- as.numeric(unlist(lapply(X = strsplit(x = trimws(peakLines), split = " "), FUN = "[", 2)))
      splash <- RMassBank:::getSplash(peaks = data.frame("m/z" = peakMzs, "int." = peakInts))
      if(splashThere != splash){
        lines2[[PK_SPLASH_idx]] <- paste("PK$SPLASH: ", splash, sep = "")
        print(paste("Peak splash", basename(filePath), filePath, sep = "\t"))
      }
    }
    
    if(!isTRUE(all.equal(lines, lines2))){
      writeLines(con = filePath, text = lines2)
      print(paste("Corrected", filePath))
    }
  }
  
  options(warn = 2)
}
renameRecord <- function(oldRecordPath, newRecordPath, newAccession){
  lines <- readLines(con = oldRecordPath)
  lines[[1]] <- paste("ACCESSION", newAccession, sep = ": ")
  writeLines(text = lines, con = newRecordPath)
  #cmd <- paste("mv '", oldRecordPath, "' '", newRecordPath, "'", sep = "")
  #print(cmd)
  #file.rename(from = oldRecordPath, to = newRecordPath)
  ##out <- system(command = cmd, intern = TRUE)
}
getIDsList <- function(allNames){
  idLogicalList <- list(
    "isNumberRich"= ((nchar(allNames) - nchar(gsub(x = allNames, pattern = "[0-9]", replacement = ""))) / nchar(allNames)) > 0.2,
    "isInChIKey" = grepl(x = allNames, pattern = "^[A-Z]{14,14}-[A-Z]{10,10}-[A-Z]$"),
    "isUNII"     = grepl(x = allNames, pattern = "^UNII-[A-Z0-9]{10,10}$"),
    "isNSC"      = grepl(x = allNames, pattern = "^NSC ?\\d+$"),
    "isNCGC"     = grepl(x = allNames, pattern = "^NCGC\\d+-\\d+$"),
    #isFDA      = grepl(x = allNames, pattern = "^[A-Z0-9]{10,10}$"),
    "isEINECS"   = grepl(x = allNames, pattern = "^EINECS \\d+-\\d+-\\d+$"),
    "isPubChem"  = grepl(x = allNames, pattern = "^((CID )|(PubChem))\\d+$"),
    "isCHEMBL"   = grepl(x = allNames, pattern = "^S?CHEMBL\\d+$"),
    "isChEBI"    = grepl(x = allNames, pattern = "^CHEBI:\\d+$"),
    "isKEGG"     = grepl(x = allNames, pattern = "^C\\d{5,5}$"),
    "isCAS"      = grepl(x = allNames, pattern = "^(CAS-)?\\d+-\\d+-\\d+$"),
    "isSLING"    = grepl(x = allNames, pattern = "^\\d+-EP\\d+[A-Z]\\d$"),
    "isABI_Chem" = grepl(x = allNames, pattern = "^AC[A-Z0-9]{6,6}$"),
    "isAchemica" = grepl(x = allNames, pattern = "^ACMC-[a-zA-Z0-9]{5,6}$"),
    "isACN"      = grepl(x = allNames, pattern = "^ACN-S\\d{6,6}$"),
    "isAC"       = grepl(x = allNames, pattern = "^AC-\\d+$"),
    "isAB"       = grepl(x = allNames, pattern = "^AB\\d+([_-]\\d+)?$"),
    "isA"        = grepl(x = allNames, pattern = "^A[- ]?\\d+$"),
    "is4CN"      = grepl(x = allNames, pattern = "^4CN-\\d+$"),
    "isAcon1"    = grepl(x = allNames, pattern = "^ACon1_\\d{6,6}$"),
    "isAI3"      = grepl(x = allNames, pattern = "^AI3-\\d{5,5}$"),
    "isAJ"       = grepl(x = allNames, pattern = "^AJ-\\d+$"),
    "isAK"       = grepl(x = allNames, pattern = "^(AK\\d+)|(AK-\\d+)$"),
    "isAKOS"     = grepl(x = allNames, pattern = "^AKOS\\d{9,9}$"),
    "isALBB"     = grepl(x = allNames, pattern = "^ALBB-\\d{6,6}$"),
    "isAE"       = grepl(x = allNames, pattern = "^AE-\\d+/\\d+$"),
    "isAKOS"     = grepl(x = allNames, pattern = "^AKOS [A-Z0-9]+-[A-Z0-9]+$"),
    "isDSSTox"   = grepl(x = allNames, pattern = "^((DSSTox_((GSID)|(CID)|(RID))_)|(DTXSID))\\d+$"),
    "isDB"       = grepl(x = allNames, pattern = "^DB-?\\d+$"),
    "isD"        = grepl(x = allNames, pattern = "^D\\d{5,5}$"),
    "isCTK"      = grepl(x = allNames, pattern = "^CTK[A-Z0-9]{6,6}$"),
    "isCS"       = grepl(x = allNames, pattern = "^CS-[A-Z0-9]+$"),
    "isCM"       = grepl(x = allNames, pattern = "^CM\\d{5,5}$"),
    "isCJ"       = grepl(x = allNames, pattern = "^CJ-\\d{5,5}$"),
    "isCID"      = grepl(x = allNames, pattern = "^cid_\\d+$"),
    "isCCRIS"    = grepl(x = allNames, pattern = "^CCRIS \\d+$"),
    "isCCG"      = grepl(x = allNames, pattern = "^CCG-\\d+$"),
    "isCC"       = grepl(x = allNames, pattern = "^CC-\\d+$"),
    "isC"        = grepl(x = allNames, pattern = "^C-\\d{5,5}$"),
    "isBSPBio"   = grepl(x = allNames, pattern = "^BSPBio_\\d{6,6}$"),
    "isBRN"      = grepl(x = allNames, pattern = "^BRN \\d{7,7}$"),
    "isBRD"      = grepl(x = allNames, pattern = "^BRD-[A-Z]\\d{8,8}-\\d{3,3}-\\d{2,2}-\\d$"),
    "isBR"       = grepl(x = allNames, pattern = "^BR-\\d{5,5}$"),
    "isBPBio1"   = grepl(x = allNames, pattern = "^BPBio1_\\d{6,6}$"),
    "isBP"       = grepl(x = allNames, pattern = "^BP-\\d{5,5}$"),
    "isbmse"     = grepl(x = allNames, pattern = "^bmse\\d{6,6}$"),
    "isBio1"     = grepl(x = allNames, pattern = "^Bio1_\\d{6,6}$"),
    "isBIDD"     = grepl(x = allNames, pattern = "^BIDD:[A-Z]{2,2}\\d{4,4}$"),
    "isBDBM"     = grepl(x = allNames, pattern = "^BDBM\\d+$"),
    "isBCP"      = grepl(x = allNames, pattern = "^BCP\\d{5,5}$"),
    "isBC"       = grepl(x = allNames, pattern = "^BC\\d{6,6}$"),
    "isBBL"      = grepl(x = allNames, pattern = "^BBL\\d{6,6}$"),
    "isBB"       = grepl(x = allNames, pattern = "^BB \\d{7,7}$"),
    "isAX"       = grepl(x = allNames, pattern = "^AX\\d{7,7}$"),
    "isAS"       = grepl(x = allNames, pattern = "^AS-\\d{5,5}$"),
    "isANW"      = grepl(x = allNames, pattern = "^ANW-\\d{5,5}$"),
    "isAN"       = grepl(x = allNames, pattern = "^AN-\\d+$"),
    "isAM"       = grepl(x = allNames, pattern = "^AM\\d+$"),
    "isARONIS"   = grepl(x = allNames, pattern = "^ARONIS\\d+$"),
    "isAS"       = grepl(x = allNames, pattern = "^AS\\d+$"),
    "isMBZ"      = grepl(x = allNames, pattern = "^AMBZ\\d+$"),
    "isBC"       = grepl(x = allNames, pattern = "^BC\\d+$"),
    "isBCP"      = grepl(x = allNames, pattern = "^BCP\\d+$"),
    "isBI"       = grepl(x = allNames, pattern = "^BI[A-Z]+\\d+$"),
    "isC"        = grepl(x = allNames, pattern = "^C-\\d+$"),
    "isBiomol"   = grepl(x = allNames, pattern = "^BiomolK\\d+_\\d+$"),
    "isCaswell"  = grepl(x = allNames, pattern = "^Caswell No. [A-Z0-9]+$"),
    "isCMC"      = grepl(x = allNames, pattern = "^CMC_\\d+$"),
    "isDicK1c"   = grepl(x = allNames, pattern = "^DiKk1c_\\d{6,6}$"),
    "isDS"       = grepl(x = allNames, pattern = "^DS-\\d+$"),
    "isEBD"      = grepl(x = allNames, pattern = "^EBD\\d+$"),
    "isEC"       = grepl(x = allNames, pattern = "^EC \\d{3,3}-\\d{3,3}-\\d$"),
    "isCP"       = grepl(x = allNames, pattern = "^CP[A-Z]\\d+$"),
    "isEpitope"  = grepl(x = allNames, pattern = "^Epitope ID:\\d+$"),
    "isF"        = grepl(x = allNames, pattern = "^F\\d{4,4}-\\d{4,4}$"),
    "isFCH"      = grepl(x = allNames, pattern = "^FCH\\d{6,6}$"),
    "isFEMA"     = grepl(x = allNames, pattern = "^FEMA (((No.)|(Number)) )?\\d{4,4}$"),
    "isFT"       = grepl(x = allNames, pattern = "^FT-\\d{7,7}$"),
    "isGTPL"     = grepl(x = allNames, pattern = "^GTPL\\d{4,4}$"),
    "isHMS"      = grepl(x = allNames, pattern = "^HMS\\d+[A-Z]\\d{2,2}$"),
    "isHSDB"     = grepl(x = allNames, pattern = "^HSDB \\d+$"),
    "isHTS"      = grepl(x = allNames, pattern = "^HTS\\d{6,6}$"),
    "isJ"        = grepl(x = allNames, pattern = "^((J-\\d{6,6})|(J\\d{5,5}))$"),
    "isInChI"    = grepl(x = allNames, pattern = "^InChI=.+$"),
    "isIDI1"     = grepl(x = allNames, pattern = "^IDI1_\\d{6,6}$"),
    "isI"        = grepl(x = allNames, pattern = "^I\\d+-\\d+$"),
    "isI2"       = grepl(x = allNames, pattern = "^I-?\\d+$"),
    "isHY"       = grepl(x = allNames, pattern = "^HY-[A-Z]\\d+$"),
    "isKSC"      = grepl(x = allNames, pattern = "^KSC[A-Z0-9]{6,6}$"),
    "isKS"       = grepl(x = allNames, pattern = "^KS-[A-Z0-9]{8,8}$"),
    "isKBioGR"   = grepl(x = allNames, pattern = "^KBio((GR)|(SS)|(\\d+))_\\d{6,6}$"),
    "isK"        = grepl(x = allNames, pattern = "^K[ -]?\\d+$"),
    "isKS"       = grepl(x = allNames, pattern = "^KS-(([A-Z0-9]{8,8})|(\\d{4,4}))$"),
    "isKSC"      = grepl(x = allNames, pattern = "^KSC[A-Z0-9]{6,6}$"),
    "isL"        = grepl(x = allNames, pattern = "^L\\d+$"),
    "isQ"        = grepl(x = allNames, pattern = "^Q-\\d{6,6}$"),
    "isPS"       = grepl(x = allNames, pattern = "^PS-\\d{4,4}$"),
    "isPrifac"   = grepl(x = allNames, pattern = "^Prifr?ac \\d+$"),
    "isPrestwick"= grepl(x = allNames, pattern = "^Prestwick\\d?_\\d{6,6}$"),
    "isPharmakon"= grepl(x = allNames, pattern = "^Pharmakon\\d+-\\d{8,8}$"),
    "isPelargidanon"= grepl(x = allNames, pattern = "^pelargidanon \\d{4,4}$"),
    "isPDSP"     = grepl(x = allNames, pattern = "^PDSP\\d_\\d{6,6}$"),
    "isP"        = grepl(x = allNames, pattern = "^P\\d{4,4}$"),
    "isOprea"    = grepl(x = allNames, pattern = "^Oprea\\d_\\d{6,6}$"),
    "isOpera"    = grepl(x = allNames, pattern = "^Opera_ID_\\d+$"),
    "isNSC"      = grepl(x = allNames, pattern = "^NSC-\\d+$"),
    "isNINDS"    = grepl(x = allNames, pattern = "^NINDS_\\d{6,6}$"),
    "isNE"       = grepl(x = allNames, pattern = "^NE\\d{5,5}$"),
    "isNCI"      = grepl(x = allNames, pattern = "^NCI\\d{4,4}$"),
    "isNCGC"     = grepl(x = allNames, pattern = "^NCGC\\d+$"),
    "isNC"       = grepl(x = allNames, pattern = "^NC\\d{5,5}$"),
    "isN"        = grepl(x = allNames, pattern = "^N\\d{4,4}$"),
    "isMP"       = grepl(x = allNames, pattern = "^MP-\\d{4,4}$"),
    "isMolPort"  = grepl(x = allNames, pattern = "^MolPort-\\d{3,3}-\\d{3,3}-\\d{3,3}$"),
    "isMLS"      = grepl(x = allNames, pattern = "^MLS\\d{9,9}$"),
    "isMFCD"     = grepl(x = allNames, pattern = "^MFCD\\d{8,8}( \\(\\d+\\+?%\\))?$"),
    "isMEG"      = grepl(x = allNames, pattern = "^MEG[a-z0-9]{3,3}_\\d{6,6}$"),
    "isMCULE"    = grepl(x = allNames, pattern = "^MCULE-\\d{10,10}$"),
    "isM"        = grepl(x = allNames, pattern = "^M-\\d{4,4}$"),
    "isM2"       = grepl(x = allNames, pattern = "^M\\d+$"),
    "isLS"       = grepl(x = allNames, pattern = "^LS-?\\d+$"),
    "isLM"       = grepl(x = allNames, pattern = "^LM[A-Z]{2,2}\\d+$"),
    "isLM2"      = grepl(x = allNames, pattern = "^LM[ -]\\d+$"),
    "isLP"       = grepl(x = allNames, pattern = "^LP\\d+$"),
    "isSMR"      = grepl(x = allNames, pattern = "^SMR\\d{9,9}$"),
    "isSMP"      = grepl(x = allNames, pattern = "^SMP_\\d{6,6}$"),
    "isSDCCGMLS" = grepl(x = allNames, pattern = "^SDCCGMLS-\\d{7,7}.P\\d{3,3}$"),
    "isSC"       = grepl(x = allNames, pattern = "^SC-\\d{5,5}$"),
    "isSBI"      = grepl(x = allNames, pattern = "^SBI-\\d{7,7}.P\\d{3,3}$"),
    "isSBB"      = grepl(x = allNames, pattern = "^SBB\\d{6,6}$"),
    "isSAM"      = grepl(x = allNames, pattern = "^SAM\\d{9,9}$"),
    "isS"        = grepl(x = allNames, pattern = "^s\\d{4,4}$"),
    "isS2"       = grepl(x = allNames, pattern = "^S\\d+$"),
    "isRTR"      = grepl(x = allNames, pattern = "^RTR-\\d{6,6}$"),
    "isRTC"      = grepl(x = allNames, pattern = "^RTC-\\d{6,6}$"),
    "isRP"       = grepl(x = allNames, pattern = "^RP\\d{5,5}$"),
    "isRL"       = grepl(x = allNames, pattern = "^RL\\d{5,5}$"),
    "isNCI"      = grepl(x = allNames, pattern = "^NCI60_\\d{6,6}$"),
    "isRW"       = grepl(x = allNames, pattern = "^RW\\d{4,4}$"),
    "isSEL"      = grepl(x = allNames, pattern = "^SEL\\d{8,8}$"),
    "isSMP"      = grepl(x = allNames, pattern = "^SMP\\d?_\\d{6,6}$"),
    "isSPBio"    = grepl(x = allNames, pattern = "^SPBio_\\d{6,6}$"),
    "isSpectrum" = grepl(x = allNames, pattern = "^((Spectrum)|(SPECTRUM))\\d?_?\\d{6,6}$"),
    "isSR"       = grepl(x = allNames, pattern = "^SR-\\d{11,11}(-\\d+)?$"),
    "isST"       = grepl(x = allNames, pattern = "^ST\\d+$"),
    "isSTK"      = grepl(x = allNames, pattern = "^STK\\d{6,6}$"),
    "isSTL"      = grepl(x = allNames, pattern = "^STL\\d{6,6}$"),
    "isSTR"      = grepl(x = allNames, pattern = "^STR\\d{5,5}$"),
    "isJsp"      = grepl(x = allNames, pattern = "^Jsp\\d{6,6}$")
  )
  return(idLogicalList)
}
getPubChemSynonyms <- function(cids, namesSuggested = NULL, namesToRetain = NULL){
  ## get synonyms
  synonyms_m <- unique(unlist(sapply(X = cids, FUN = function(cid){
    ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/16219824/synonyms/JSON
    url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid, "/synonyms/JSON", sep = "")
    synonyms = tryCatch({    fromJSON(url)[[1]][[1]]$Synonym[[1]]  }, error = function(e) {  print(e); return(NA)  })
    if(all(!is.na(synonyms), synonyms==0)) synonyms <- NA
    return(synonyms)
  })))
  ## normalize
  synonyms_m <- synonyms_m[!is.na(synonyms_m)]
  synonyms_m <- trimws(unique(c(namesSuggested, synonyms_m)))
  synonyms_m <- unlist(strsplit(x = synonyms_m, split = "; "))
  synonyms_m <- removeItalicsTagsFromNames(synonyms_m)
  synonyms_m <- gsub(x = synonyms_m, pattern = ", [<>=]*\\d+%$", replacement = "")
  synonyms_m <- gsub(x = synonyms_m, pattern = "’", replacement = "'")
  synonyms_m <- gsub(x = synonyms_m, pattern = "→", replacement = "->")
  
  ## duplicated
  synonyms_m <- synonyms_m[!duplicated(toupper(synonyms_m))]
  ## prefer Ajmalicine over Ajmalicin and 16-hydroxy-tabersonine over 16-hydroxytabersonine
  synonyms_m <- sort(synonyms_m, decreasing = T)
  synonyms_m <- synonyms_m[order(nchar(synonyms_m))]
  synonyms_m <- synonyms_m[!duplicated(gsub(x = toupper(synonyms_m), pattern = "E$", replacement = ""))]
  synonyms_m <- synonyms_m[!duplicated(gsub(x = toupper(synonyms_m), pattern = "\\-", replacement = ""))]
  
  ## remove IDs
  idLogicalList <- getIDsList(synonyms_m)
  isID <- logical(length = length(synonyms_m))
  for(idLogical in idLogicalList) isID <- isID | idLogical
  #sum(isID)
  #synonyms_mSave <- synonyms_m
  synonyms_m <- synonyms_m[!isID]
  
  ## select names
  synonyms_m <- synonyms_m[order(nchar(synonyms_m))]
  if(length(synonyms_m) > 5)
    synonyms_m <- c(synonyms_m[1:2], synonyms_m[seq(from = 3, to = length(synonyms_m), length.out = 3)])
  ## retain names and remove duplicated
  synonyms_m <- c(namesToRetain, synonyms_m)
  synonyms_m <- synonyms_m[!duplicated(toupper(synonyms_m))]
  synonyms_m <- synonyms_m[!duplicated(gsub(x = toupper(synonyms_m), pattern = "E$", replacement = ""))]
  synonyms_m <- synonyms_m[!duplicated(gsub(x = toupper(synonyms_m), pattern = "\\-", replacement = ""))]
  return(synonyms_m)
}
enrichRecordWithNames <- function(folder){
  #folder <- "/mnt/data/IPB/Projects/2017_005_MS-databases/mFam contributions/Stefanie Doell IPB Halle_2nd/"
  recordFolder <- paste(folder, "converted to MassBank", sep = "/")
  if(!file.exists(recordFolder)) recordFolder <- folder
  filePaths <- list.files(path = recordFolder, pattern = "[A-Z0-9]{8,8}\\.txt$", full.names=TRUE)
  print(paste("Enriching", length(filePaths), "files in folder", folder))
  
  library("jsonlite")
  for(filePath in filePaths){
    lines <- readLines(con = filePath)
    indecesNames <- which(grepl(x = lines, pattern = "^CH\\$NAME"))
    names <- gsub(x = lines[indecesNames], pattern = "^CH\\$NAME: ", replacement = "")
    
    indecesCIDs <- which(grepl(x = lines, pattern = "^CH\\$LINK: PUBCHEM"))
    cids <- gsub(x = lines[indecesCIDs], pattern = "^CH\\$LINK: PUBCHEM ", replacement = "")
    if(length(cids)==0){
      print(paste("no cid", filePath))
      next
    }
    
    synonyms <- getPubChemSynonyms(cids = cids, namesToRetain = names)
    print(paste(basename(filePath), ": ", paste(synonyms, collapse = "; "), sep = ""))
    linesNew <- paste("CH$NAME: ", synonyms, sep = "")
    lines2 <- c(
      lines[seq_len(min(indecesNames) - 1)], 
      linesNew, 
      lines[seq(from = max(indecesNames) + 1, to = length(lines))]
    )
    writeLines(con = filePath, text = lines2)
  }
}
enrichRecordWithNames_save <- function(folder){
  #folder <- "/mnt/data/IPB/Projects/2017_005_MS-databases/mFam contributions/Stefanie Doell IPB Halle_2nd/"
  recordFolder <- paste(folder, "converted to MassBank", sep = "/")
  if(!file.exists(recordFolder)) recordFolder <- folder
  filePaths <- list.files(path = recordFolder, pattern = "[A-Z0-9]{8,8}\\.txt$", full.names=TRUE)
  print(paste("Correcting", length(filePaths), "files in folder", folder))
  allNames <- NULL
  allCIDs <- NULL
  for(filePath in filePaths){
    lines <- readLines(con = filePath)
    indeces <- which(grepl(x = lines, pattern = "^CH\\$NAME"))
    names <- gsub(x = lines[indeces], pattern = "^CH\\$NAME: ", replacement = "")
    allNames <- c(allNames, names)
    
    indeces <- which(grepl(x = lines, pattern = "^CH\\$LINK: PUBCHEM"))
    cids <- gsub(x = lines[indeces], pattern = "^CH\\$LINK: PUBCHEM ", replacement = "")
    allCIDs <- c(allCIDs, cids)
  }
  allNames <- sort(unique(allNames))
  idLogicalList <- list(
    "isInChIKey" = grepl(x = allNames, pattern = "^[A-Z]{14,14}-[A-Z]{10,10}-[A-Z]$"),
    "isUNII"     = grepl(x = allNames, pattern = "^UNII-[A-Z0-9]{10,10}$"),
    "isNSC"      = grepl(x = allNames, pattern = "^NSC ?\\d+$"),
    "isNCGC"     = grepl(x = allNames, pattern = "^NCGC\\d+-\\d+$"),
    #isFDA      = grepl(x = allNames, pattern = "^[A-Z0-9]{10,10}$"),
    "isEINECS"   = grepl(x = allNames, pattern = "^EINECS \\d+-\\d+-\\d+$"),
    "isPubChem"  = grepl(x = allNames, pattern = "^((CID )|(PubChem))\\d+$"),
    "isCHEMBL"   = grepl(x = allNames, pattern = "^S?CHEMBL\\d+$"),
    "isChEBI"    = grepl(x = allNames, pattern = "^CHEBI:\\d+$"),
    "isKEGG"     = grepl(x = allNames, pattern = "^C\\d{5,5}$"),
    "isCAS"      = grepl(x = allNames, pattern = "^(CAS-)?\\d+-\\d+-\\d+$"),
    "isSLING"    = grepl(x = allNames, pattern = "^\\d+-EP\\d+[A-Z]\\d$"),
    "isABI_Chem" = grepl(x = allNames, pattern = "^AC[A-Z0-9]{6,6}$"),
    "isAchemica" = grepl(x = allNames, pattern = "^ACMC-[a-zA-Z0-9]{5,6}$"),
    "isACN"      = grepl(x = allNames, pattern = "^ACN-S\\d{6,6}$"),
    "isAC"       = grepl(x = allNames, pattern = "^AC-\\d+$"),
    "isAB"       = grepl(x = allNames, pattern = "^AB\\d+([_-]\\d+)?$"),
    "isA"        = grepl(x = allNames, pattern = "^A[- ]?\\d+$"),
    "is4CN"      = grepl(x = allNames, pattern = "^4CN-\\d+$"),
    "isAcon1"    = grepl(x = allNames, pattern = "^ACon1_\\d{6,6}$"),
    "isAI3"      = grepl(x = allNames, pattern = "^AI3-\\d{5,5}$"),
    "isAJ"       = grepl(x = allNames, pattern = "^AJ-\\d+$"),
    "isAK"       = grepl(x = allNames, pattern = "^(AK\\d+)|(AK-\\d+)$"),
    "isAKOS"     = grepl(x = allNames, pattern = "^AKOS\\d{9,9}$"),
    "isALBB"     = grepl(x = allNames, pattern = "^ALBB-\\d{6,6}$"),
    "isAE"       = grepl(x = allNames, pattern = "^AE-\\d+/\\d+$"),
    "isAKOS"     = grepl(x = allNames, pattern = "^AKOS [A-Z0-9]+-[A-Z0-9]+$"),
    "isDSSTox"   = grepl(x = allNames, pattern = "^((DSSTox_((GSID)|(CID)|(RID))_)|(DTXSID))\\d+$"),
    "isDB"       = grepl(x = allNames, pattern = "^DB-?\\d+$"),
    "isD"        = grepl(x = allNames, pattern = "^D\\d{5,5}$"),
    "isCTK"      = grepl(x = allNames, pattern = "^CTK[A-Z0-9]{6,6}$"),
    "isCS"       = grepl(x = allNames, pattern = "^CS-[A-Z0-9]+$"),
    "isCM"       = grepl(x = allNames, pattern = "^CM\\d{5,5}$"),
    "isCJ"       = grepl(x = allNames, pattern = "^CJ-\\d{5,5}$"),
    "isCID"      = grepl(x = allNames, pattern = "^cid_\\d+$"),
    "isCCRIS"    = grepl(x = allNames, pattern = "^CCRIS \\d+$"),
    "isCCG"      = grepl(x = allNames, pattern = "^CCG-\\d+$"),
    "isCC"       = grepl(x = allNames, pattern = "^CC-\\d+$"),
    "isC"        = grepl(x = allNames, pattern = "^C-\\d{5,5}$"),
    "isBSPBio"   = grepl(x = allNames, pattern = "^BSPBio_\\d{6,6}$"),
    "isBRN"      = grepl(x = allNames, pattern = "^BRN \\d{7,7}$"),
    "isBRD"      = grepl(x = allNames, pattern = "^BRD-[A-Z]\\d{8,8}-\\d{3,3}-\\d{2,2}-\\d$"),
    "isBR"       = grepl(x = allNames, pattern = "^BR-\\d{5,5}$"),
    "isBPBio1"   = grepl(x = allNames, pattern = "^BPBio1_\\d{6,6}$"),
    "isBP"       = grepl(x = allNames, pattern = "^BP-\\d{5,5}$"),
    "isbmse"     = grepl(x = allNames, pattern = "^bmse\\d{6,6}$"),
    "isBio1"     = grepl(x = allNames, pattern = "^Bio1_\\d{6,6}$"),
    "isBIDD"     = grepl(x = allNames, pattern = "^BIDD:[A-Z]{2,2}\\d{4,4}$"),
    "isBDBM"     = grepl(x = allNames, pattern = "^BDBM\\d+$"),
    "isBCP"      = grepl(x = allNames, pattern = "^BCP\\d{5,5}$"),
    "isBC"       = grepl(x = allNames, pattern = "^BC\\d{6,6}$"),
    "isBBL"      = grepl(x = allNames, pattern = "^BBL\\d{6,6}$"),
    "isBB"       = grepl(x = allNames, pattern = "^BB \\d{7,7}$"),
    "isAX"       = grepl(x = allNames, pattern = "^AX\\d{7,7}$"),
    "isAS"       = grepl(x = allNames, pattern = "^AS-\\d{5,5}$"),
    "isANW"      = grepl(x = allNames, pattern = "^ANW-\\d{5,5}$"),
    "isAN"       = grepl(x = allNames, pattern = "^AN-\\d+$"),
    "isAM"       = grepl(x = allNames, pattern = "^AM\\d+$"),
    "isARONIS"   = grepl(x = allNames, pattern = "^ARONIS\\d+$"),
    "isAS"       = grepl(x = allNames, pattern = "^AS\\d+$"),
    "isMBZ"      = grepl(x = allNames, pattern = "^AMBZ\\d+$"),
    "isBC"       = grepl(x = allNames, pattern = "^BC\\d+$"),
    "isBCP"      = grepl(x = allNames, pattern = "^BCP\\d+$"),
    "isBI"       = grepl(x = allNames, pattern = "^BI[A-Z]+\\d+$"),
    "isC"        = grepl(x = allNames, pattern = "^C-\\d+$"),
    "isBiomol"   = grepl(x = allNames, pattern = "^BiomolK\\d+_\\d+$"),
    "isCaswell"  = grepl(x = allNames, pattern = "^Caswell No. [A-Z0-9]+$"),
    "isCMC"      = grepl(x = allNames, pattern = "^CMC_\\d+$"),
    "isDicK1c"   = grepl(x = allNames, pattern = "^DiKk1c_\\d{6,6}$"),
    "isDS"       = grepl(x = allNames, pattern = "^DS-\\d+$"),
    "isEBD"      = grepl(x = allNames, pattern = "^EBD\\d+$"),
    "isEC"       = grepl(x = allNames, pattern = "^EC \\d{3,3}-\\d{3,3}-\\d$"),
    "isCP"       = grepl(x = allNames, pattern = "^CP[A-Z]\\d+$"),
    "isEpitope"  = grepl(x = allNames, pattern = "^Epitope ID:\\d+$"),
    "isF"        = grepl(x = allNames, pattern = "^F\\d{4,4}-\\d{4,4}$"),
    "isFCH"      = grepl(x = allNames, pattern = "^FCH\\d{6,6}$"),
    "isFEMA"     = grepl(x = allNames, pattern = "^FEMA (((No.)|(Number)) )?\\d{4,4}$"),
    "isFT"       = grepl(x = allNames, pattern = "^FT-\\d{7,7}$"),
    "isGTPL"     = grepl(x = allNames, pattern = "^GTPL\\d{4,4}$"),
    "isHMS"      = grepl(x = allNames, pattern = "^HMS\\d+[A-Z]\\d{2,2}$"),
    "isHSDB"     = grepl(x = allNames, pattern = "^HSDB \\d+$"),
    "isHTS"      = grepl(x = allNames, pattern = "^HTS\\d{6,6}$"),
    "isJ"        = grepl(x = allNames, pattern = "^((J-\\d{6,6})|(J\\d{5,5}))$"),
    "isInChI"    = grepl(x = allNames, pattern = "^InChI=.+$"),
    "isIDI1"     = grepl(x = allNames, pattern = "^IDI1_\\d{6,6}$"),
    "isI"        = grepl(x = allNames, pattern = "^I\\d+-\\d+$"),
    "isI2"       = grepl(x = allNames, pattern = "^I-?\\d+$"),
    "isHY"       = grepl(x = allNames, pattern = "^HY-[A-Z]\\d+$"),
    "isKSC"      = grepl(x = allNames, pattern = "^KSC[A-Z0-9]{6,6}$"),
    "isKS"       = grepl(x = allNames, pattern = "^KS-[A-Z0-9]{8,8}$"),
    "isKBioGR"   = grepl(x = allNames, pattern = "^KBio((GR)|(SS)|(\\d+))_\\d{6,6}$"),
    "isK"        = grepl(x = allNames, pattern = "^K[ -]?\\d+$"),
    "isKS"       = grepl(x = allNames, pattern = "^KS-(([A-Z0-9]{8,8})|(\\d{4,4}))$"),
    "isKSC"      = grepl(x = allNames, pattern = "^KSC[A-Z0-9]{6,6}$"),
    "isL"        = grepl(x = allNames, pattern = "^L\\d+$"),
    "isQ"        = grepl(x = allNames, pattern = "^Q-\\d{6,6}$"),
    "isPS"       = grepl(x = allNames, pattern = "^PS-\\d{4,4}$"),
    "isPrifac"   = grepl(x = allNames, pattern = "^Prifr?ac \\d+$"),
    "isPrestwick"= grepl(x = allNames, pattern = "^Prestwick\\d?_\\d{6,6}$"),
    "isPharmakon"= grepl(x = allNames, pattern = "^Pharmakon\\d+-\\d{8,8}$"),
    "isPelargidanon"= grepl(x = allNames, pattern = "^pelargidanon \\d{4,4}$"),
    "isPDSP"     = grepl(x = allNames, pattern = "^PDSP\\d_\\d{6,6}$"),
    "isP"        = grepl(x = allNames, pattern = "^P\\d{4,4}$"),
    "isOprea"    = grepl(x = allNames, pattern = "^Oprea\\d_\\d{6,6}$"),
    "isOpera"    = grepl(x = allNames, pattern = "^Opera_ID_\\d+$"),
    "isNSC"      = grepl(x = allNames, pattern = "^NSC-\\d+$"),
    "isNINDS"    = grepl(x = allNames, pattern = "^NINDS_\\d{6,6}$"),
    "isNE"       = grepl(x = allNames, pattern = "^NE\\d{5,5}$"),
    "isNCI"      = grepl(x = allNames, pattern = "^NCI\\d{4,4}$"),
    "isNCGC"     = grepl(x = allNames, pattern = "^NCGC\\d+$"),
    "isNC"       = grepl(x = allNames, pattern = "^NC\\d{5,5}$"),
    "isN"        = grepl(x = allNames, pattern = "^N\\d{4,4}$"),
    "isMP"       = grepl(x = allNames, pattern = "^MP-\\d{4,4}$"),
    "isMolPort"  = grepl(x = allNames, pattern = "^MolPort-\\d{3,3}-\\d{3,3}-\\d{3,3}$"),
    "isMLS"      = grepl(x = allNames, pattern = "^MLS\\d{9,9}$"),
    "isMFCD"     = grepl(x = allNames, pattern = "^MFCD\\d{8,8}( \\(\\d+\\+?%\\))?$"),
    "isMEG"      = grepl(x = allNames, pattern = "^MEG[a-z0-9]{3,3}_\\d{6,6}$"),
    "isMCULE"    = grepl(x = allNames, pattern = "^MCULE-\\d{10,10}$"),
    "isM"        = grepl(x = allNames, pattern = "^M-\\d{4,4}$"),
    "isM2"       = grepl(x = allNames, pattern = "^M\\d+$"),
    "isLS"       = grepl(x = allNames, pattern = "^LS-?\\d+$"),
    "isLM"       = grepl(x = allNames, pattern = "^LM[A-Z]{2,2}\\d+$"),
    "isLM2"      = grepl(x = allNames, pattern = "^LM[ -]\\d+$"),
    "isLP"       = grepl(x = allNames, pattern = "^LP\\d+$"),
    "isSMR"      = grepl(x = allNames, pattern = "^SMR\\d{9,9}$"),
    "isSMP"      = grepl(x = allNames, pattern = "^SMP_\\d{6,6}$"),
    "isSDCCGMLS" = grepl(x = allNames, pattern = "^SDCCGMLS-\\d{7,7}.P\\d{3,3}$"),
    "isSC"       = grepl(x = allNames, pattern = "^SC-\\d{5,5}$"),
    "isSBI"      = grepl(x = allNames, pattern = "^SBI-\\d{7,7}.P\\d{3,3}$"),
    "isSBB"      = grepl(x = allNames, pattern = "^SBB\\d{6,6}$"),
    "isSAM"      = grepl(x = allNames, pattern = "^SAM\\d{9,9}$"),
    "isS"        = grepl(x = allNames, pattern = "^s\\d{4,4}$"),
    "isS2"       = grepl(x = allNames, pattern = "^S\\d+$"),
    "isRTR"      = grepl(x = allNames, pattern = "^RTR-\\d{6,6}$"),
    "isRTC"      = grepl(x = allNames, pattern = "^RTC-\\d{6,6}$"),
    "isRP"       = grepl(x = allNames, pattern = "^RP\\d{5,5}$"),
    "isRL"       = grepl(x = allNames, pattern = "^RL\\d{5,5}$"),
    "isNCI"      = grepl(x = allNames, pattern = "^NCI60_\\d{6,6}$"),
    "isRW"       = grepl(x = allNames, pattern = "^RW\\d{4,4}$"),
    "isSEL"      = grepl(x = allNames, pattern = "^SEL\\d{8,8}$"),
    "isSMP"      = grepl(x = allNames, pattern = "^SMP\\d?_\\d{6,6}$"),
    "isSPBio"    = grepl(x = allNames, pattern = "^SPBio_\\d{6,6}$"),
    "isSpectrum" = grepl(x = allNames, pattern = "^((Spectrum)|(SPECTRUM))\\d?_?\\d{6,6}$"),
    "isSR"       = grepl(x = allNames, pattern = "^SR-\\d{11,11}(-\\d+)?$"),
    "isST"       = grepl(x = allNames, pattern = "^ST\\d+$"),
    "isSTK"      = grepl(x = allNames, pattern = "^STK\\d{6,6}$"),
    "isSTL"      = grepl(x = allNames, pattern = "^STL\\d{6,6}$"),
    "isSTR"      = grepl(x = allNames, pattern = "^STR\\d{5,5}$"),
    "isJsp"      = grepl(x = allNames, pattern = "^Jsp\\d{6,6}$")
  )
  if(FALSE){
    print(unlist(lapply(idLogicalList, sum)))
    print(names(idLogicalList)[unlist(lapply(idLogicalList, sum))==0])
    #print(allNames[idLogicalList$"isHSDB"])
    
    isID <- logical(length = length(allNames))
    for(idLogical in idLogicalList) isID <- isID | idLogical
    sum(isID)
    allNames[!isID][2001:3000]
    
    allCIDs <- unique(allCIDs)
    library("jsonlite")
    synonyms_m <- unique(unlist(sapply(X = allCIDs, FUN = function(cid){
      ## https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/16219824/synonyms/JSON
      url <- paste("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid, "/synonyms/JSON", sep = "")
      synonyms = tryCatch({    fromJSON(url)[[1]][[1]]$Synonym[[1]]  }, error = function(e) {  print(e); return(NA)  })
      if(all(!is.na(synonyms), synonyms==0)) synonyms <- NA
      return(synonyms)
    })))
    synonyms_m <- synonyms_m[!is.na(synonyms_m)]
    allNames <- synonyms_m
  }
  
}
aggregateSpectra_all <- function(parentFolder, aggregationFolder){
  contributorFolders <- list.files(path = parentFolder, full.names = TRUE)
  aggregateSpectra_ready(contributorFolders, aggregationFolder, tag = "all")
}
aggregateSpectra_ready <- function(contributorFolders, aggregationFolder, tag){
  library("stringr")
  
  ## msps
  mspFolders <- paste(contributorFolders, "converted to msp", sep = "/")
  mspFolders <- mspFolders[file.exists(mspFolders)]
  mspsPos <- list.files(path = mspFolders, pattern = "_pos.msp$", full.names = TRUE)
  mspsNeg <- list.files(path = mspFolders, pattern = "_neg.msp$", full.names = TRUE)
  
  ## records
  recordFolders  <- paste(contributorFolders, "converted to MassBank", sep = "/")
  recordFolders  <- recordFolders[file.exists(recordFolders)]
  recordFiles    <- lapply(X = recordFolders, FUN = list.files, full.names = TRUE, pattern = "^[A-Z]{2,3}\\d{5,6}\\.txt")
  recordCount    <- unlist(lapply(X = recordFiles, FUN = length))
  recordFolders  <- recordFolders[recordCount > 0]
  recordFiles    <- recordFiles[recordCount > 0]
  accessionCodes <- unlist(lapply(X = recordFiles, FUN = function(records){
    unique(str_extract(string = basename(records), pattern = "^[A-Z]{2,3}"))
  }))
  
  allRecordPaths <- unlist(recordFiles)
  #isDuplicated <- duplicated(basename(allRecordPaths))
  #if(any(isDuplicated)){
  #  print(paste("There are", sum(isDuplicated), "accession conflicts"))
  #  duplicatedRecordPaths <- allRecordPaths[isDuplicated]
  #  if(any(table(basename(duplicatedRecordPaths)) > 1)) stop("Accession conflicts need to be solved in a clever way / iteratively")
  #  for(duplicatedRecordPath in duplicatedRecordPaths){
  #    accessionCode <- str_extract(string = basename(duplicatedRecordPath), pattern = "^[A-Z]{2,3}")
  #    suggestedAccession <- gsub(x = basename(duplicatedRecordPath), pattern = "\\.txt$", replacement = "")
  #    while(paste(suggestedAccession, ".txt", sep = "") %in% basename(allRecordPaths)){
  #      accNumber <- as.integer(gsub(x = suggestedAccession, pattern = paste("^", accessionCode, sep = ""), replacement = ""))
  #      accNumber <- accNumber + 1
  #      suggestedAccession <- sprintf("%s%06d", accessionCode, accNumber)
  #    }
  #    newFileName <- paste(suggestedAccession, ".txt", sep = "")
  #    newRecordPath <- paste(dirname(duplicatedRecordPath), newFileName, sep = "/")
  #    renameRecord(oldRecordPath = duplicatedRecordPath, newRecordPath = newRecordPath, newAccession = suggestedAccession)
  #  }
  #  
  #  ## update
  #  recordFiles    <- lapply(X = recordFolders, FUN = list.files, full.names = TRUE, pattern = "^[A-Z]{2,3}\\d{5,6}\\.txt")
  #  recordCount    <- unlist(lapply(X = recordFiles, FUN = length))
  #  recordFolders  <- recordFolders[recordCount > 0]
  #  recordFiles    <- recordFiles[recordCount > 0]
  #  accessionCodes <- unlist(lapply(X = recordFiles, FUN = function(records){
  #    unique(str_extract(string = basename(records), pattern = "^[A-Z]{2,3}"))
  #  }))
  #  
  #  allRecordPaths <- unlist(recordFiles)
  #  isDuplicated <- duplicated(basename(allRecordPaths))
  #  if(sum(isDuplicated) > 0) stop("Accession conflicts after accession conflict resolve")
  #} else {
  #  print("There are no accession conflicts")
  #}
  
  ## correct spectra
  #for(contributorFolder in contributorFolders)    correctRecords(folder = contributorFolder)
  #for(contributorFolder in contributorFolders)    validateRecords(folder = contributorFolder)
  
  ## aggregate
  currentFolder <- paste(aggregationFolder, "/", Sys.Date(), "_", tag, sep = "")
  if(file.exists(currentFolder)) if(unlink(currentFolder, recursive = TRUE)) stop(paste("Could not remove folder", currentFolder))
  if(!dir.create(currentFolder)) stop(paste("Could not create folder", currentFolder))
  print(paste("New aggregation folder:", currentFolder))
  
  mspFolder    <- paste(currentFolder, "converted to msp", sep = "/")
  recordFolder <- paste(currentFolder, "converted to MassBank", sep = "/")
  if(!dir.create(mspFolder))    stop(paste("Could not create folder", mspFolder))
  if(!dir.create(recordFolder)) stop(paste("Could not create folder", recordFolder))
  
  accessionCodes <- gsub(x = basename(allRecordPaths), pattern = "\\d+\\.txt", replacement = "")
  newAccessions <- vector(mode = "character", length = length(allRecordPaths))
  newRecordPaths <- vector(mode = "character", length = length(allRecordPaths))
  for(accessionCode in unique(accessionCodes)){
    newAccessionsHere <- sprintf("%s%06d", accessionCode, seq(from = 1, to = sum(accessionCodes == accessionCode)))
    newAccessions[accessionCodes == accessionCode] <- newAccessionsHere
    accessionCodeFolder <- paste(recordFolder, "/", accessionCode, sep = "")
    if(!file.exists(accessionCodeFolder)) if(!dir.create(accessionCodeFolder))    stop(paste("Could not create folder", accessionCodeFolder))
    newRecordPaths[accessionCodes == accessionCode] <- paste(accessionCodeFolder, "/", newAccessionsHere, ".txt", sep = "")
  }
  print(paste("Copying", length(allRecordPaths), "records..."))
  for(idx in seq_along(allRecordPaths))
    renameRecord(oldRecordPath = allRecordPaths[[idx]], newRecordPath = newRecordPaths[[idx]], newAccession = newAccessions[[idx]])
  
  ## aggregate records
  #for(idx in seq_along(recordFiles)){
  #  recordFilesHere   <- recordFiles[[idx]]
  #  accessionCodeHere <- accessionCodes[[idx]]
  #  recordSubFolder <- paste(recordFolder, accessionCodeHere, sep = "/")
  #  print(paste("Copying to", recordSubFolder))
  #  if(!file.exists(recordSubFolder)) if(!dir.create(recordSubFolder))    stop(paste("Could not create folder", recordSubFolder))
  #  newRecordPaths <- paste(recordSubFolder, basename(recordFilesHere), sep = "/")
  #  file.copy(from = recordFilesHere, to = newRecordPaths)
  #}
  
  ## aggregate msps
  print("Aggregating msp files...")
  mspPos_all <- paste(mspFolder, "mFam_pos.msp", sep = "/")
  mspNeg_all <- paste(mspFolder, "mFam_neg.msp", sep = "/")
  msp_all    <- paste(mspFolder, "mFam.msp", sep = "/")
  mspLinesPos <- NULL
  mspLinesNeg <- NULL
  for(mspPos in mspsPos) mspLinesPos <- c(mspLinesPos, readLines(con = mspPos))
  for(mspNeg in mspsNeg) mspLinesNeg <- c(mspLinesNeg, readLines(con = mspNeg))
  print(paste("Write file", mspPos_all))
  writeLines(text = mspLinesPos, con = mspPos_all)
  print(paste("Write file", mspNeg_all))
  writeLines(text = mspLinesNeg, con = mspNeg_all)
  print(paste("Write file", msp_all))
  writeLines(text = c(mspLinesPos, mspLinesNeg), con = msp_all)
  
  ## finally validate again
  contributorFoldersHere <- paste(recordFolder, "/", unique(accessionCodes), sep = "")
  print(paste("Validating", length(contributorFoldersHere), "folders with", length(newRecordPaths), "files"))
  for(contributorFolderHere in contributorFoldersHere){
    resultsWithError <- validateRecords(contributorFolderHere)
    if(length(resultsWithError) > 0){
      print(paste("Error in", length(resultsWithError), "files"))
      print(paste(basename(names(resultsWithError))))
      #print(resultsWithError)
      print(lapply(X = resultsWithError, FUN = function(lines){lines[grepl(x = lines, pattern = "\\[main\\] ERROR")]}))
    }
  }
}
getMonoisotopicMassFromFormula <- function(formula){
  m <- rcdk::parse.smiles(molecularFormulaToSMILES(formula))[[1]]
  do.aromaticity(m)
  do.typing(m)
  do.isotopes(m)
  return(get.exact.mass(m))
}
escapeCommandLineChars <- function(str){
  chars <- c(" ", "|", "<", ">", "\"", "\'", "\'", "\`", "$", ",", ";", "*", "#")
  #str <- " "
  for(char in chars)
    str <- gsub(x = str, pattern = char, replacement = paste("\\", char, sep = ""), fixed=TRUE)
  str <- gsub(x = str, pattern = "\\\\", replacement = "\\", fixed=TRUE)
  return(str)
}
validateRecords <- function(folder){
  folderMassBank   <- paste(folder, "converted to MassBank", sep = "/")
  if(!file.exists(folderMassBank)) folderMassBank <- folder
  #cmd <- paste("/bin/bash /Massbank/MassBank-data/.scripts/validate.sh '", folderMassBank, "'", sep = "")
  #cmd <- paste("/bin/bash /Massbank/MassBank-data/.scripts/validate.sh \"", folderMassBank, "\"", sep = "")
  #cmd <- paste("/bin/bash /Massbank/MassBank-data/.scripts/validate.sh ", escapeCommandLineChars(folderMassBank), "", sep = "")
  #cmd <- paste("bash /Massbank/MassBank-data/.scripts/validate.sh '", escapeCommandLineChars(folderMassBank), "'", sep = "")
  #out <- system(command = cmd, intern = TRUE)
  #out <- system2(command = "/bin/bash /Massbank/MassBank-data/.scripts/validate.sh", args = folderMassBank, stdout = TRUE)
  
  files <- list.files(path = folderMassBank, full.names = TRUE)
  if(length(files)==0) return(character())
  cat(paste("\nValidating", length(files), "records in folder", folder, "\n"))
  results <- lapply(X = files, FUN = function(file){
    cmd <- paste("/bin/bash /Massbank/MassBank-web/MassBank-Project/MassBank-lib/target/MassBank-lib-0.0.1-default/MassBank-lib-0.0.1/bin/Validator ", gsub(x = file, pattern = " ", replacement = "\\\\ "), "", sep = "")
    out <- suppressWarnings(expr = {
      system(command = cmd, intern = TRUE)
    })
    return(out)
  })
  names(results) <- files
  errorHappened <- unlist(lapply(X = results, FUN = function(result){any(grepl(x = result, pattern = "ERROR"))}))
  errorHappenedSomeWhere <- any(errorHappened)
  errorFiles <- which(errorHappened)
  resultsWithError <- results[errorFiles]
  
  
  #aFile <- "/mnt/data/IPB/Projects/2017_005_MS-databases/mFam contributions/Corey Broeckling Colorado State University/converted to MassBank/CL000001.txt"
  #cmd <- paste("/bin/bash /Massbank/MassBank-data/.scripts/validate.sh '", aFile, "'", sep = "")
  #cmd <- paste("/bin/bash /Massbank/MassBank-web/MassBank-Project/MassBank-lib/target/MassBank-lib-0.0.1-default/MassBank-lib-0.0.1/bin/Validator '", aFile, "'", sep = "")
  
  #files <- list.files(path = folderMassBank, full.names = TRUE)
  #cmd <- paste("/bin/bash /Massbank/MassBank-web/MassBank-Project/MassBank-lib/target/MassBank-lib-0.0.1-default/MassBank-lib-0.0.1/bin/Validator ", paste("'", files, "'", sep = "", collapse = " "), "", sep = "")
  return(resultsWithError)
}
inchiKeysToInchiKeysBlock1 <- function(inchiKeys){
  inchiKeysSplitted <- strsplit(x = inchiKeys, split = "-")
  inchiKeysBlock1 <- unlist(lapply(X = inchiKeysSplitted, FUN = function(x){
    if(length(x) == 0){
      return("")
    } else {
      return(x[[1]])
    }
  }))
  return(inchiKeysBlock1)
}
createSunBurstPlot <- function(folder, fileAnnoTable){
  if(is.null(fileAnnoTable)) return()
  
  folderMassBank   <- paste(folder, "converted to MassBank", sep = "/")
  if(!file.exists(folderMassBank)) folderMassBank <- folder
  
  #folders <- c(
  #  "/mnt/data/IPB/Projects/2017_005_MS-databases/mFam contributions/Stefanie Doell IPB Halle/converted to MassBank/",
  #  "/mnt/data/IPB/Projects/2017_005_MS-databases/mFam contributions/Stefanie Doell IPB Halle_2nd//converted to MassBank/"
  #)
  inchiKeys <- sapply(X = list.files(folderMassBank, full.names = TRUE), FUN = function(file){
    lines <- readLines(file)
    line  <- lines[grepl(x = lines, pattern = "^CH\\$LINK: INCHIKEY ")]
    inchiKey <- gsub(x = line, pattern = "^CH\\$LINK: INCHIKEY ", replacement = "")
    return(inchiKey)
  })
  
  #source("/home/htreutle/Code/Java/MetFam_util/Classifier/SubstanceClassClassifier.R")
  #fileAnnoTable <- "/home/htreutle/Downloads/MetSWATH/MONA/181019_MSMS_merge_HR_scaffolds.tsv"
  #fileAnnoTable <- "/home/htreutle/Downloads/MetSWATH/MONA/190523_MSMS_HR_someScaffolds.tsv"
  annoTable <- read.table(file = fileAnnoTable, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
  structuresPrefixes <- inchiKeysToInchiKeysBlock1(inchiKeys)
  structuresPrefixesThere <- inchiKeysToInchiKeysBlock1(annoTable[, "InChIKey"])
  mapping <- match(x = structuresPrefixes, table = structuresPrefixesThere)
  structureWithEntry <- !is.na(mapping)
  
  print(paste(sum(!structureWithEntry), "/", length(inchiKeys), "spectra have no annotation entry,", sum(structureWithEntry), "remain"))
  inchiKeys <- inchiKeys[structureWithEntry]
  mapping   <- mapping  [structureWithEntry]
  annoTable <- annoTable[mapping, ]
  
  structureWithAnno <- !(annoTable[, "CHEMONT_name"] == "NA" | annoTable[, "CHEMONT_name"] == "")
  print(paste(sum(!structureWithAnno, na.rm = TRUE) + sum(is.na(structureWithAnno)), "/", length(structureWithAnno), "spectra have no annotation,", sum(structureWithAnno), "remain"))
  
  inchiKeys       <- inchiKeys[structureWithAnno]
  annoTable       <- annoTable[structureWithAnno, ]
  
  substanceclasses   <- annoTable[, "CHEMONT_name"]
  tab <- table(substanceclasses)
  
  classifierClasses <- names(tab)
  numberOfSpectra   <- tab
  
  #source("/home/htreutle/Code/Java/MetFam_util/Annotation/Sunburst_plots.R")
  
  #setwd("/home/htreutle/Downloads/tmp/")
  timeStamp <- gsub(pattern = "[ ]", x = Sys.time(), replacement = "_")
  widthToHeightRatio = 10/10;
  size = 15
  resolution = 300
  width = size * widthToHeightRatio
  height = size
  #fileName = paste("Sunburst_", basename(classifierFile), ".png", sep = "")
  fileName = paste(folder, "/", timeStamp, "_SubstanceClasses_", length(classifierClasses), ".png", sep = "")
  print(fileName)
  png(filename = fileName, res = resolution, width = width * resolution, height = height * resolution)
  sunBurstPlotFromSubstanceClasses(classifierClasses, numberOfSpectra, pValues = NULL)
  dev.off()
}

sunBurstPlotFromSubstanceClasses <- function(classifierClasses, numberOfSpectra, pValues = NULL, qualities = NULL, qualitiesMinMax = NULL, qualitiesName = NULL,
                                             ## colors
                                             colorStart = 0, colorAlpha = 1,
                                             ## thresholds
                                             degreeThresholdForDrawing = 0.5,
                                             minimumAngleToShowSegmentText      = 15,
                                             minimumAngleToShowSegmentTextSmall =  3,
                                             minimumAngleToShowSegmentTextOuter  =  1,
                                             ## cex
                                             plotCexSegmentText  = 0.8,
                                             plotCexSegmentTextSmall = 1,
                                             plotCexTopLevelText = 1,
                                             max.level = 10,
                                             #is.redundant = FALSE
                                             show.counts = TRUE,
                                             qualityColor = "darkgreen",
                                             pValueColor = "green"
){
  if(FALSE){
    pValues = NULL
    qualities = NULL
    qualitiesMinMax = NULL
    colorStart = 0
    colorAlpha = 1
    degreeThresholdForDrawing = 0.5
    plotCexSegmentText  = 0.8
    plotCexSegmentTextSmall = 1
    plotCexTopLevelText = 1
    minimumAngleToShowSegmentText      = 15
    minimumAngleToShowSegmentTextSmall =  3
    minimumAngleToShowSegmentTextOuter  =  1
    max.level = 10
    #is.redundant = FALSE
    show.counts = TRUE
    qualityColor = "darkgreen"
    pValueColor = "green"
  }
  
  #pValues = rep(x = 1, times = length(classifierClasses))
  #pValues[sample(x = length(classifierClasses), size = 50, replace = F)] = runif(min = 0, max = 0.05, n = 50)
  noPValues   <- is.null(pValues)
  noQualities <- is.null(qualities)
  level <- unlist(lapply(X = strsplit(x = classifierClasses, split = "; "), FUN = length))
  
  maximumImportance <- 1 / 0.001
  palette <- NULL
  if(!noPValues)
    palette <- colorRampPalette(c('blue', ifelse(test = !noQualities, yes = qualityColor, no = pValueColor)))
  if(!noQualities)
    palette <- colorRampPalette(c('white', ifelse(test = !noQualities, yes = qualityColor, no = pValueColor)))
  if(!is.null(palette))
    colors2  <- palette(maximumImportance)
  
  #if(is.redundant){
  #  for(idx in seq_along(classifierClasses)){
  #    superClasses <- classifierClasses %in% sapply(X = seq_len(length(strsplit(x = classifierClasses[[idx]], split = "; ")[[1]]) - 1), FUN = function(idx2){paste(strsplit(x = classifierClasses[[idx]], split = "; ")[[1]][seq_len(idx2)], collapse = "; ")})
  #    ## ...
  #  }
  #}
  
  ## all class levels
  classesAndSubClasses <- lapply(X = strsplit(x = classifierClasses, split = "; "), FUN = function(x){
    sapply(X = seq_along(x), FUN = function(y){paste(x[1:y], collapse = "; ")})
  })
  classesByLevel <- list()
  labelsByLevel <- list()
  for(levelHere in seq_len(max(level))){
    classesByLevel[[levelHere]] <- sort(unique(unlist(lapply(X = classesAndSubClasses, FUN = function(y){
      if(length(y) < levelHere) return(NULL)
      else return(y[[levelHere]])
    }))))
    labelsByLevel[[levelHere]] <- unlist(lapply(X = strsplit(x = classesByLevel[[levelHere]], split = "; "), FUN = tail, n=1))
  }
  
  ## class counts
  countsByLevel <- list()
  for(levelHere in rev(seq_len(max(level)))){
    countsByLevel[[levelHere]] <- unlist(lapply(X = classesByLevel[[levelHere]], FUN = function(class){
      newSpectra <- ifelse(test = class %in% classifierClasses, yes = numberOfSpectra[[which(class == classifierClasses)]], no = 0)
      oldSpectra <- ifelse(test = levelHere < max(level), yes = sum(countsByLevel[[levelHere+1]][grepl(x = classesByLevel[[levelHere+1]], pattern = paste("^", class, sep = ""))]), no = 0)
      return(newSpectra + oldSpectra)
    }))
  }
  rootCount <- sum(countsByLevel[[1]])
  
  ## coordinates
  colors <- rainbow(n = 1000, start = colorStart, alpha = colorAlpha)
  startDegreeByLevel <- list()
  spanDegreeByLevel <- list()
  colorByLevel <- list()
  importanceByLevel <- list()
  for(levelHere in seq_len(max(level))){
    startDegreeByLevel[[levelHere]] <- list()
    spanDegreeByLevel[[levelHere]] <- list()
    colorByLevel[[levelHere]] <- list()
    importanceByLevel[[levelHere]] <- list()
    
    classesToProcess <- classesByLevel[[levelHere]]
    precursorClasses <- NULL
    if(levelHere == 1)  precursorClasses <- ""
    else                precursorClasses <- classesByLevel[[levelHere-1]]
    
    for(precursorClassIdx in seq_along(precursorClasses)){
      precursorClass <- precursorClasses[[precursorClassIdx]]
      
      print(paste(levelHere, precursorClassIdx, precursorClass))
      
      classesToProcessHereSelection <- grepl(x = classesToProcess, pattern = precursorClass)
      classesToProcessHere <- classesToProcess[classesToProcessHereSelection]
      startDegree <- ifelse(test = levelHere == 1, yes = 0, no = startDegreeByLevel[[levelHere-1]][[precursorClassIdx]])
      scalingFactor <- 1
      for(classToProcessHere in classesToProcessHere){
        classIdx <- which(classesByLevel[[levelHere]] == classToProcessHere)
        degreeSpan <- 360 * countsByLevel[[levelHere]][[classIdx]] / rootCount * scalingFactor
        startDegreeByLevel[[levelHere]][[classIdx]] <- startDegree
        spanDegreeByLevel [[levelHere]][[classIdx]] <- degreeSpan
        colorByLevel      [[levelHere]][[classIdx]] <- colors[[(floor(startDegree + degreeSpan / 2) / 360 * length(colors)) + ifelse(test = (floor(startDegree + degreeSpan / 2) / 360 * length(colors))==length(colors), yes = 0, no = 1) ]]
        importanceByLevel [[levelHere]][[classIdx]] <- ifelse(test = noPValues & noQualities, 
                                                              yes = NA, 
                                                              no =  ifelse(test = !noPValues, 
                                                                           yes = ifelse(test = classToProcessHere %in% classifierClasses, yes = 1 / pValues[classifierClasses == classToProcessHere], no = NA),
                                                                           no =  ifelse(test = !noQualities, 
                                                                                        yes = ifelse(test = classToProcessHere %in% classifierClasses, yes = qualities[classifierClasses == classToProcessHere] / qualitiesMinMax[[2]], no = NA),
                                                                                        no = NA
                                                                           )
                                                              )
        )
        if((!noPValues) & !is.na(importanceByLevel[[levelHere]][[classIdx]]))
          if(importanceByLevel[[levelHere]][[classIdx]] > maximumImportance)
            importanceByLevel [[levelHere]][[classIdx]] <- maximumImportance
        startDegree <- startDegree + degreeSpan
      }
    }
  }
  spanDegreeByLevel <- lapply(X = spanDegreeByLevel, FUN = unlist)
  
  thereIsNextLevelByLevel <- list()
  for(levelHere in seq_len(max(level))){
    thereIsNextLevelByLevel[[levelHere]] <- list()
    if(levelHere == max(level)){
      thereIsNextLevelByLevel[[levelHere]] <- rep(x = FALSE, times = length(classesByLevel[[levelHere]]))
    } else {
      for(classIdx in seq_along(classesByLevel[[levelHere]]))
        #thereIsNextLevelByLevel[[levelHere]][[classIdx]] <- any(grepl(x = classesByLevel[[levelHere+1]], pattern = classesByLevel[[levelHere]][[classIdx]]) & spanDegreeByLevel[[levelHere]][[classIdx]] < degreeThresholdForDrawing)
        thereIsNextLevelByLevel[[levelHere]][[classIdx]] <- any(grepl(x = classesByLevel[[levelHere+1]][  spanDegreeByLevel[[levelHere+1]] >= degreeThresholdForDrawing  ], pattern = classesByLevel[[levelHere]][[classIdx]]))
    }
  }
  
  ## If you use it in published research, please cite:
  ##  Gu, Z. circlize implements and enhances circular visualization in R. Bioinformatics 2014.
  library("circlize")
  library("plotrix")
  desat <- function(cols, sat=0.5) {
    X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
  }
  
  plotRadius <- min(max(level), max.level) + 1.5
  wordWrapExpansionFactor <- 1.2
  
  if(any(level > max.level)){
    classesByLevel          <- classesByLevel         [seq_len(max.level)]
    labelsByLevel           <- labelsByLevel          [seq_len(max.level)]
    countsByLevel           <- countsByLevel          [seq_len(max.level)]
    startDegreeByLevel      <- startDegreeByLevel     [seq_len(max.level)]
    spanDegreeByLevel       <- spanDegreeByLevel      [seq_len(max.level)]
    colorByLevel            <- colorByLevel           [seq_len(max.level)]
    thereIsNextLevelByLevel <- thereIsNextLevelByLevel[seq_len(max.level)]
    for(idx in seq_along(thereIsNextLevelByLevel[[max.level]]))
      thereIsNextLevelByLevel[[max.level]][[idx]] <- FALSE
  }
  
  par("mar"=c(0,0,0,0))
  plot(1, type="n", xlab="", ylab="", xlim=c(-plotRadius, plotRadius), ylim=c(-plotRadius, plotRadius), axes = FALSE)
  
  ## collect data
  tmp <- unlist(lapply(X = seq_along(classesByLevel), FUN = function(levelHere){
    lapply(X = seq_along(classesByLevel[[levelHere]]), FUN = function(classIdx){
      if(spanDegreeByLevel[[levelHere]][[classIdx]] < degreeThresholdForDrawing)  return(NULL)
      list(
        start.degree = startDegreeByLevel[[levelHere]][[classIdx]], 
        end.degree = startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]], 
        rou1 = levelHere - 1, 
        rou2 = levelHere, 
        center = c(0,0), 
        clock.wise = FALSE, 
        col = colorByLevel[[levelHere]][[classIdx]], 
        importance = ifelse(test = !noPValues, yes = as.integer(importanceByLevel[[levelHere]][[classIdx]]), no = importanceByLevel[[levelHere]][[classIdx]])#"white"
      )
    })
  }), recursive = FALSE)
  tmp <- tmp[!unlist(lapply(X = tmp, FUN = is.null))]
  tmp <- tmp[order(unlist(lapply(X = tmp, FUN = "[", 8)))]
  
  ## draw sectors
  if(noPValues & noQualities)
    tmp <- lapply(X = tmp, FUN = function(list){
      draw.sector(
        start.degree = list[["start.degree"]], 
        end.degree   = list[["end.degree"]], 
        rou1         = list[["rou1"]], 
        rou2         = list[["rou2"]], 
        center       = list[["center"]], 
        clock.wise   = list[["clock.wise"]], 
        col    = list[["col"]],
        lwd    = 1,
        border = "white"
      )
    })
  
  pValueImportanceToColor <- function(importance){
    if(importance == 1) return("grey")
    if(importance == 1000) return(ifelse(test = !noPValues, yes = pValueColor, no = qualityColor))
    #desat(cols = pValueColor, sat = (importance + maximumImportance*0.35) / (maximumImportance + maximumImportance*0.35))
    #desat(cols = pValueColor, sat = log10(importance) / log10(maximumImportance))
    #colors2[[importance]]
    colors2[[ (log10(importance) / log10(maximumImportance)) * maximumImportance ]]
  }
  if(!noPValues)
    tmp <- lapply(X = tmp, FUN = function(list){
      draw.sector(
        start.degree = list[["start.degree"]], 
        end.degree   = list[["end.degree"]], 
        rou1         = list[["rou1"]], 
        rou2         = list[["rou2"]], 
        center       = list[["center"]], 
        clock.wise   = list[["clock.wise"]], 
        col    = ifelse(test = is.na(list[["importance"]]), yes = "grey", no = 
                          ifelse(test = list[["importance"]] == 1, 
                                 #yes = desat(cols = list[["col"]], sat = 0.1), 
                                 yes = "grey", 
                                 no  = pValueImportanceToColor(list[["importance"]])
                          )),
        lwd    = 1,#ifelse(test = list[["importance"]] == 1, yes = 1, no = ceiling(log10(list[["importance"]]))),
        border = "black"#ifelse(test = list[["importance"]] == 1, yes = "white", no = "black")#ifelse(test = list[["importance"]] == 1, yes = "white", no = colorGradient[[list[["importance"]]]])#"white"
      )
    })
  if(!noQualities)
    tmp <- lapply(X = tmp, FUN = function(list){
      draw.sector(
        start.degree = list[["start.degree"]], 
        end.degree   = list[["end.degree"]], 
        rou1         = list[["rou1"]], 
        rou2         = list[["rou2"]], 
        center       = list[["center"]], 
        clock.wise   = list[["clock.wise"]], 
        col    = ifelse(test = is.na(list[["importance"]]), yes = "grey", no = 
                          ifelse(test = list[["importance"]] == 1, 
                                 yes = "white", 
                                 no = colors2[[list[["importance"]] * maximumImportance]]
                                 #no = desat(cols = qualityColor, sat = list[["importance"]] / qualitiesMinMax[[2]])
                          )),
        lwd    = 1,#ifelse(test = list[["importance"]] == 1, yes = 1, no = ceiling(log10(list[["importance"]]))),
        border = "black"#ifelse(test = list[["importance"]] == 1, yes = "white", no = colorGradient[[list[["importance"]]]])#"white"
      )
    })
  
  ## segment text
  tmp <- sapply(X = seq_along(classesByLevel), FUN = function(levelHere){
    sapply(X = seq_along(classesByLevel[[levelHere]]), FUN = function(classIdx){
      if(spanDegreeByLevel[[levelHere]][[classIdx]] < minimumAngleToShowSegmentTextSmall)  return()
      verboseText <- spanDegreeByLevel[[levelHere]][[classIdx]] >= minimumAngleToShowSegmentText
      if(verboseText) {
        ## class plus class count
        textTokens <- strwrap(x = labelsByLevel[[levelHere]][[classIdx]], width = max(nchar(strsplit(x = "Something", split = " ")[[1]])))
        if(show.counts) textTokens <- c(textTokens, paste("(", as.character(countsByLevel[[levelHere]][[classIdx]]), ")", sep = "", collapse = ""))
        
        ## merge short words
        while(all(length(textTokens) > 1, any(sapply(X = seq(from=1, to = length(textTokens) - 1), FUN = function(idx){nchar(paste(textTokens[c(idx, idx+1)], collapse = " "))}) <= (max(nchar(textTokens)) * wordWrapExpansionFactor)))){
          lengths <- sapply(X = seq(from=1, to = length(textTokens) - 1), FUN = function(idx){nchar(paste(textTokens[c(idx, idx+1)], collapse = " "))})
          lengthsOK <- lengths <= (max(nchar(textTokens)) * wordWrapExpansionFactor)
          idx <- which(lengths == max(lengths[lengthsOK]))[[1]]
          textTokens[c(idx, idx+1)] <- paste(textTokens[c(idx, idx+1)], collapse = " ")
          textTokens <- textTokens[-idx]
        }
        
        middle <- (startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2)/360 * 2*pi
        for(idx in seq_along(textTokens)){
          offset <-  ifelse(test = middle > pi, yes = (length(textTokens) - idx + 1 + 0.) / (length(textTokens) + 1), no = (idx + 0.) / (length(textTokens) + 1))
          if(levelHere == 1 & length(classesByLevel[[levelHere]]==1)){
            text(
              x = 0, y = 0.5 - offset, labels = textTokens[[idx]], adj = c(0.5,0.5), 
              #srt=srt,
              cex = plotCexTopLevelText
            )
          } else {
            isSwitched <- middle > pi
            arctext(
              x = textTokens[[idx]], 
              center = c(ifelse(test = isSwitched, yes = 0, no = 0), 0), radius = levelHere - offset,# - 0.04, 
              middle = middle, 
              cex = plotCexSegmentText, stretch = 1,
              clockwise = !isSwitched
            )
          }
        }
      } else {
        if(show.counts){
          ## only class count
          radius <- levelHere - 0.5# + 0.2
          angle <- (startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2)/360 * 2*pi
          x <- radius * sin(angle+pi/2)
          y <- radius * cos(angle+pi/2)
          srt <- startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2
          isSwitched <- srt > 90 & srt < 270
          srt <- ifelse(test = isSwitched, yes = srt + 180, no = srt)
          text(
            x = x, y = -y, labels = as.character(countsByLevel[[levelHere]][[classIdx]]), adj = c(0.5,0.5), 
            srt=srt,
            cex = plotCexTopLevelText
          )
        }
      }
    })
  })
  ## outer text
  levelMaxHere <- max(level)
  
  tmp <- unlist(lapply(X = seq_along(classesByLevel), FUN = function(levelHere){
    #if(levelHere > levelMaxHere)  return(NULL)
    lapply(X = seq_along(classesByLevel[[levelHere]]), FUN = function(classIdx){
      if(spanDegreeByLevel[[levelHere]][[classIdx]] < degreeThresholdForDrawing)  return(NULL)
      if(spanDegreeByLevel[[levelHere]][[classIdx]] >= minimumAngleToShowSegmentText)  return(NULL)
      #if(spanDegreeByLevel[[levelHere]][[classIdx]] < minimumAngleToShowSegmentTextOuter)  return()
      if(thereIsNextLevelByLevel[[levelHere]][[classIdx]]){
        if(sum(spanDegreeByLevel[[levelHere+1]][  spanDegreeByLevel[[levelHere+1]] >= degreeThresholdForDrawing  ][grepl(x = classesByLevel[[levelHere+1]][  spanDegreeByLevel[[levelHere+1]] >= degreeThresholdForDrawing  ], pattern = classesByLevel[[levelHere]][[classIdx]])]) / spanDegreeByLevel[[levelHere]][[classIdx]] > 0.4)
          return(NULL)
      }
      radius <- levelHere + 0.2
      angle <- (startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2)/360 * 2*pi
      x <- radius * sin(angle+pi/2)
      y <- radius * cos(angle+pi/2)
      srt <- startDegreeByLevel[[levelHere]][[classIdx]] + spanDegreeByLevel[[levelHere]][[classIdx]]/2
      isSwitched <- srt > 90 & srt < 270
      adj <- c(ifelse(test = isSwitched, yes = 1, no = 0), 0.5)
      srt <- ifelse(test = isSwitched, yes = srt + 180, no = srt)
      return(list(
        angle = angle,
        x = x, y = -y, labels = labelsByLevel[[levelHere]][[classIdx]], adj = adj, 
        srt=srt#,
        #cex = plotCexTopLevelText
      ))
    })
  }), recursive = FALSE)
  tmp <- tmp[!unlist(lapply(X = tmp, FUN = is.null))]
  angles <- unlist(lapply(X = tmp, FUN = "[[", "angle")) * 2 * pi
  
  order <- order(angles)
  tmp <- tmp[order]
  angles <- angles[order]
  
  anglesDiff <- diff(angles)
  retain <- rep(x = TRUE, times = length(angles))
  for(idx in seq_along(anglesDiff)){
    if(anglesDiff[[idx]] < minimumAngleToShowSegmentTextOuter){
      retain[[idx+1]] <- FALSE
      if(idx < length(anglesDiff))
        anglesDiff[[idx + 1]] <- anglesDiff[[idx]] + anglesDiff[[idx + 1]]
    }
  }
  tmp <- tmp[retain]
  lapply(X = tmp, FUN = function(list){
    text(
      x = list[["x"]], y = list[["y"]], labels = list[["labels"]], adj = list[["adj"]], 
      srt=list[["srt"]],
      cex = plotCexTopLevelText
    )
  })
  
  ## legend
  if(!noQualities | !noPValues){
    library("squash")
    #colfunc <- colorRampPalette(c("white", "red"))
    #legend_image <- as.raster(matrix(colfunc(20), ncol=1))
    #rasterImage(legend_image, 0, 0, 1,1)
    legendAnchor <- levelMaxHere - 1.
    colorMap <- makecmap(x = c(0, 1), n = 100, colFn = palette)
    
    #legend_imageLFC <- as.raster(x = t(x = matrix(data = cmap(x = seq(from = qualitiesMinMax[[1]], to = qualitiesMinMax[[2]], length.out = 100), map = colorMap), nrow=1)))
    legend_imageLFC <- as.raster(x = t(x = matrix(data = cmap(x = seq(
      from = ifelse(test = !noPValues, yes = 1,                                         no = qualitiesMinMax[[2]]), 
      #to   = ifelse(test = !noPValues, yes = -(log10(0.05) / log10(maximumImportance)), no = qualitiesMinMax[[1]]), 
      to   = ifelse(test = !noPValues, yes = 0, no = qualitiesMinMax[[1]]), 
      length.out = 100), map = colorMap), nrow=1)))
    rasterImage(image = legend_imageLFC, xleft = legendAnchor-2, ybottom = legendAnchor-1, xright = legendAnchor-1, ytop = legendAnchor)
    
    if(!noPValues){
      lfcLegendPositionsLabels <- c(0.05, 0.025, 0.01, 0.005, 0.001)
      lfcLegendPositions       <- (((log10(lfcLegendPositionsLabels) - log10(0.001)) / (log10(maximumImportance) + log10(0.05))) )
    }
    if(!noQualities){
      lfcLegendPositionsLabels <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
      lfcLegendPositions       <- rev(c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
    }
    graphics::text(x = legendAnchor - 1 + .2, y = legendAnchor - 1 + rev(lfcLegendPositions), labels = lfcLegendPositionsLabels, pos = 4)
    
    segments(## lower hori; upper hori; left vert; right vert
      x0  = c(0, 0, 0, 1) + legendAnchor - 2,
      x1  = c(1, 1, 0, 1) + legendAnchor - 2,
      y0  = c(0, 1, 0, 0) + legendAnchor - 1,
      y1  = c(0, 1, 1, 1) + legendAnchor - 1
    )
    ## axis marks
    segments(
      x0  = rep(x = legendAnchor - 1 - 0.1, times = length(lfcLegendPositions)),
      x1  = rep(x = legendAnchor - 1 + 0.2, times = length(lfcLegendPositions)),
      y0  = legendAnchor - 1 + lfcLegendPositions,
      y1  = legendAnchor - 1 + lfcLegendPositions
    )
    graphics::text(x = legendAnchor - 2, y = legendAnchor + 0.09, labels = ifelse(test = !noPValues, yes = "p-value", no = qualitiesName), adj = c(0,0))
  }
}

##Taken from MassBankAdministrationScripts/R projects/MBrecordToFile/getInfo_and_InChIKeys_extended_fullrecord.R
##
##Directory is the name of the directory
##csvname is the designated name that the csv-file will get
##babel_dir is the directory path (no spaces!) containing obabel.exe
getInfoFixKey <- function(Directory, csvname, babel_dir, verbose.output=TRUE){
  Files <- list.files(Directory, pattern="*.txt", full.names=TRUE, recursive=TRUE)
  if(length(Files) == 0){
    cat(paste("### Warning ### No tsv generated: No record files in folder", Directory, "\n"))
    return()
  }
  #recursive=TRUE should get all sub-dirs
  #need to add pattern to skip the mols and tsvs
  #
  ## Good to know how many files will be processed
  if(verbose.output) print(paste0("We will process ",length(Files), " Files!"))
  
  wantedmat <- matrix(0,length(Files),(38))
  for(i in seq_along(Files)){
    fileConnection <- file(normalizePath(Files[i]))
    record <- readLines(fileConnection)
    close(fileConnection)
    
    ## Check if fields contain NAs
    CSIDTRUE <- grep('CH$LINK: CHEMSPIDER',record, value = TRUE, fixed = TRUE)
    CSIDFALSE <- "N/A"
    CASTRUE <- grep('CH$LINK: CAS',record, value = TRUE, fixed = TRUE)
    CASFALSE <- "N/A"
    CIDTRUE <- grep('CH$LINK: PUBCHEM',record, value = TRUE, fixed = TRUE)
    CIDFALSE <- "N/A"
    CSIDTRUE <- grep('CH$LINK: CHEMSPIDER',record, value = TRUE, fixed = TRUE)
    CSIDFALSE <- "N/A"
    INCHIKEYTRUE <- grep('CH$LINK: INCHIKEY',record, value = TRUE, fixed = TRUE)
    INCHIKEYFALSE <- "N/A"
    INCHIKEY <- ifelse(length(INCHIKEYTRUE)==1, substring(grep('CH$LINK: INCHIKEY',record, value = TRUE, fixed = TRUE),19), INCHIKEYFALSE)
    SMILES <- substring(grep('CH$SMILES:',record, value = TRUE, fixed = TRUE),12)
    SMILES_NA <- grep("N/A",SMILES, value = TRUE, fixed = TRUE)
    INCHIKEY_NA <- grep("N/A",INCHIKEY, value = TRUE, fixed = TRUE)
    INSTRUMENT_TYPE_TRUE <- grep('AC$INSTRUMENT_TYPE:',record, value = TRUE, fixed = TRUE)
    INSTRUMENT_TYPE_FALSE <- "N/A"
    INSTRUMENT_TRUE <- grep('AC$INSTRUMENT:',record, value = TRUE, fixed = TRUE)
    INSTRUMENT_FALSE <- "N/A"
    RESOLUTION_TRUE <- grep('AC$MASS_SPECTROMETRY: RESOLUTION',record, value = TRUE, fixed = TRUE)
    RESOLUTION_FALSE <- "N/A"
    MS_TYPE_TRUE <- grep('AC$MASS_SPECTROMETRY: MS_TYPE',record, value = TRUE, fixed = TRUE)
    MS_TYPE_FALSE <- "N/A"
    PRECURSOR_TYPE_TRUE <- grep('MS$FOCUSED_ION: PRECURSOR_TYPE',record, value = TRUE, fixed = TRUE)
    PRECURSOR_TYPE_FALSE <- "N/A"
    PRECURSOR_MZ_TRUE <- grep('MS$FOCUSED_ION: PRECURSOR_M/Z',record, value = TRUE, fixed = TRUE)
    PRECURSOR_MZ_FALSE <- "N/A"
    BASE_PEAK_TRUE <- grep('MS$FOCUSED_ION: BASE_PEAK',record, value = TRUE, fixed = TRUE)
    BASE_PEAK_FALSE <- "N/A"
    IONIZATION_TRUE <- grep('AC$MASS_SPECTROMETRY: IONIZATION',record, value = TRUE, fixed = TRUE)
    IONIZATION_FALSE <- "N/A"
    FRAGMENTATION_MODE_TRUE <- grep('AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE',record, value = TRUE, fixed = TRUE)
    FRAGMENTATION_MODE_FALSE <- "N/A"
    COLL_E_TRUE <- grep('AC$MASS_SPECTROMETRY: COLLISION_ENERGY',record, value = TRUE, fixed = TRUE)
    COLL_E_FALSE <- "N/A"
    COLUMN_NAME_TRUE <- grep('AC$CHROMATOGRAPHY: COLUMN_NAME',record, value = TRUE, fixed = TRUE)
    COLUMN_NAME_FALSE <- "N/A"
    FLOW_GRADIENT_TRUE <- grep('AC$CHROMATOGRAPHY: FLOW_GRADIENT',record, value = TRUE, fixed = TRUE)
    FLOW_GRADIENT_FALSE <- "N/A"
    FLOW_RATE_TRUE <- grep('AC$CHROMATOGRAPHY: FLOW_RATE',record, value = TRUE, fixed = TRUE)
    FLOW_RATE_FALSE <- "N/A"
    RETENTION_TIME_TRUE <- grep('AC$CHROMATOGRAPHY: RETENTION_TIME',record, value = TRUE, fixed = TRUE)
    RETENTION_TIME_FALSE <- "N/A"
    SOLVENT_A_TRUE <- grep('AC$CHROMATOGRAPHY: SOLVENT A',record, value = TRUE, fixed = TRUE)
    SOLVENT_A_FALSE <- "N/A"
    SOLVENT_B_TRUE <- grep('AC$CHROMATOGRAPHY: SOLVENT B',record, value = TRUE, fixed = TRUE)
    SOLVENT_B_FALSE <- "N/A"
    
    SPLASH_TRUE <- grep('PK$SPLASH:',record, value = TRUE, fixed = TRUE)
    SPLASH_FALSE <- "N/A"
    
    # Collapse MS$DATA_PROCESSING
    processing <- list()
    processing <- as.list(substring(grep('MS$DATA_PROCESSING:',record, value = TRUE, fixed = TRUE),21))
    
    DATA_PROCESSING_TRUE <- paste(processing, '#', collapse = '')
    DATA_PROCESSING_FALSE <-"N/A"
    
    # Collapse COMMENT
    comments <- list()
    comments <- as.list(substring(grep('COMMENT:',record, value = TRUE, fixed = TRUE),10))
    
    COMMENT_TRUE <- paste(comments, '#', collapse = '')
    COMMENT_FALSE <-"N/A"
    
    
    #fill in missing InChI Key where possible with Open Babel conversion from SMILES
    # Can only attempt this if SMILES exists; takes a while so don't recalculate unless necessary
    if((length(INCHIKEY_NA)==1)&&(length(SMILES_NA)!=1)) {
      new_inchikey <- create.inchikey(SMILES, babel_dir)
      if(length(new_inchikey)>=1) {
        INCHIKEY <- new_inchikey
      }
    }
    ## Good to know what R is doing
    if(verbose.output) print(paste0("In progress with #",i," of ",length(Files)," records with ACCESSION ",substring(grep('ACCESSION:',record, value = TRUE, fixed = TRUE),12),"."))
    
    ## Parse the fields from the records	  
    colnames(wantedmat) <- c("ACCESSION","AUTHORS","LICENSE","DATE","COMMENT","NAME","FORMULA","EXACT_MASS","IUPAC","INCHIKEY","SMILES","CSID","CID","CAS","INSTRUMENT","INSTRUMENT_TYPE","MS_TYPE","IONIZATION","ION_MODE","FRAGMENTATION_MODE","COLL_E","COLL_E_UNIT","RESOLUTION","COLUMN_NAME","FLOW_GRADIENT","FLOW_RATE","FLOW_RATE_UNIT","RETENTION_TIME","RETENTION_TIME_UNIT","SOLVENT_A","SOLVENT_B","BASE_PEAK","PRECURSOR_MZ","PRECURSOR_TYPE","DATA_PROCESSING","SPLASH","EULINK","JPLINK")
    ## The information block
    wantedmat[i,'ACCESSION'] <- substring(grep('ACCESSION:',record, value = TRUE, fixed = TRUE),12)
    wantedmat[i,'AUTHORS'] <- substring(grep('AUTHORS:',record, value = TRUE, fixed = TRUE),10)
    wantedmat[i,'LICENSE'] <- substring(grep('LICENSE:',record, value = TRUE, fixed = TRUE),10)
    wantedmat[i,'DATE'] <- substring(grep('DATE:',record, value = TRUE, fixed = TRUE),7)
    ifelse(length(COMMENT_TRUE)==1, wantedmat[i,'COMMENT'] <- COMMENT_TRUE, wantedmat[i,'COMMENT'] <- COMMENT_FALSE)
    chnames <- list()
    chnames <- as.list(substring(grep('CH$NAME:',record, value = TRUE, fixed = TRUE),10))
    wantedmat[i,'NAME'] <- chnames[[1]]
    #wantedmat[i,'SMILES'] <- substring(grep('CH$SMILES:',record, value = TRUE, fixed = TRUE),12)
    wantedmat[i,'FORMULA'] <- substring(grep('CH$FORMULA:',record, value = TRUE, fixed = TRUE),13)
    wantedmat[i,'EXACT_MASS'] <- substring(grep('CH$EXACT_MASS',record, value = TRUE, fixed = TRUE),16)
    wantedmat[i,'SMILES'] <- SMILES
    wantedmat[i,'IUPAC'] <- substring(grep('CH$IUPAC:',record, value = TRUE, fixed = TRUE),11)
    
    ## The next lines check if field is NA or not (for optional fields)
    #ifelse(is.na(INCHIKEYTRUE) == TRUE, wantedmat[i,'INCHIKEY'] <- INCHIKEYFALSE, wantedmat[i,'INCHIKEY'] <- substring(grep('CH$LINK: INCHIKEY',record, value = TRUE, fixed = TRUE),19))
    #ifelse(is.na(CSIDTRUE) == TRUE, wantedmat[i,'CSID'] <- CSIDFALSE, wantedmat[i,'CSID'] <- substring(grep('CH$LINK: CHEMSPIDER',record, value = TRUE, fixed = TRUE),21))
    #INCHIKEY <- ifelse(length(INCHIKEYTRUE)==1, wantedmat[i,'INCHIKEY'] <- substring(grep('CH$LINK: INCHIKEY',record, value = TRUE, fixed = TRUE),19), wantedmat[i,'INCHIKEY'] <- INCHIKEYFALSE)
    wantedmat[i,'INCHIKEY'] <- INCHIKEY
    ifelse(length(CSIDTRUE)==1, wantedmat[i,'CSID'] <- substring(grep('CH$LINK: CHEMSPIDER',record, value = TRUE, fixed = TRUE),21), wantedmat[i,'CSID'] <- CSIDFALSE)
    ifelse(length(CIDTRUE)==1, wantedmat[i,'CID'] <- gsub("[a-z, A-Z, ,:]","", substring(grep('CH$LINK: PUBCHEM',record, value = TRUE, fixed = TRUE),17)), wantedmat[i,'CID'] <- CIDFALSE)
    ifelse(length(CASTRUE)==1, wantedmat[i,'CAS'] <- as.character(substring(grep('CH$LINK: CAS',record, value = TRUE, fixed = TRUE),14)), wantedmat[i,'CAS'] <- CASFALSE)
    
    ## The instrument block
    ifelse(length(INSTRUMENT_TRUE)==1, wantedmat[i,'INSTRUMENT'] <- substring(grep('AC$INSTRUMENT:',record, value = TRUE, fixed = TRUE),16), wantedmat[i,'INSTRUMENT'] <- INSTRUMENT_FALSE)
    ifelse(length(INSTRUMENT_TYPE_TRUE)==1, wantedmat[i,'INSTRUMENT_TYPE'] <- substring(grep('AC$INSTRUMENT_TYPE:',record, value = TRUE, fixed = TRUE),21), wantedmat[i,'INSTRUMENT_TYPE'] <- INSTRUMENT_TYPE_FALSE)
    
    ## The spectroscopy block
    ifelse(length(MS_TYPE_TRUE)==1, wantedmat[i,'MS_TYPE'] <- substring(grep('AC$MASS_SPECTROMETRY: MS_TYPE',record, value = TRUE, fixed = TRUE),31), wantedmat[i,'MS_TYPE'] <- MS_TYPE_FALSE)
    ifelse(length(IONIZATION_TRUE)==1, wantedmat[i,'IONIZATION'] <- substring(grep('AC$MASS_SPECTROMETRY: IONIZATION',record, value = TRUE, fixed = TRUE),34), wantedmat[i,'IONIZATION'] <- IONIZATION_FALSE)
    wantedmat[i,'ION_MODE'] <- substring(grep('AC$MASS_SPECTROMETRY: ION_MODE',record, value = TRUE, fixed = TRUE),32)
    ifelse(length(FRAGMENTATION_MODE_TRUE)==1, wantedmat[i,'FRAGMENTATION_MODE'] <- substring(grep('AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE',record, value = TRUE, fixed = TRUE),42), wantedmat[i,'FRAGMENTATION_MODE'] <- FRAGMENTATION_MODE_FALSE)
    ifelse(length(COLL_E_TRUE)==1, wantedmat[i,'COLL_E'] <- gsub("[a-z,A-Z, ,(,),%,-]","",substring(grep('AC$MASS_SPECTROMETRY: COLLISION_ENERGY',record, value = TRUE, fixed = TRUE),40)), wantedmat[i,'COLL_E'] <- COLL_E_FALSE)
    ifelse(length(COLL_E_TRUE)==1, wantedmat[i,'COLL_E_UNIT'] <- gsub("[0-9, ,(,),.]","",substring(grep('AC$MASS_SPECTROMETRY: COLLISION_ENERGY',record, value = TRUE, fixed = TRUE),40)), wantedmat[i,'COLL_E_UNIT'] <- COLL_E_FALSE)
    ifelse(length(RESOLUTION_TRUE)==1, wantedmat[i,'RESOLUTION'] <- substring(grep('AC$MASS_SPECTROMETRY: RESOLUTION',record, value = TRUE, fixed = TRUE),34), wantedmat[i,'RESOLUTION'] <- RESOLUTION_FALSE)
    
    ## The chromatography block
    ifelse(length(COLUMN_NAME_TRUE)==1, wantedmat[i,'COLUMN_NAME'] <- substring(grep('AC$CHROMATOGRAPHY: COLUMN_NAME',record, value = TRUE, fixed = TRUE),32), wantedmat[i,'COLUMN_NAME'] <- COLUMN_NAME_FALSE)
    ifelse(length(FLOW_GRADIENT_TRUE)==1, wantedmat[i,'FLOW_GRADIENT'] <- substring(grep('AC$CHROMATOGRAPHY: FLOW_GRADIENT',record, value = TRUE, fixed = TRUE),34), wantedmat[i,'FLOW_GRADIENT'] <- FLOW_GRADIENT_FALSE)
    ifelse(length(FLOW_RATE_TRUE)==1, wantedmat[i,'FLOW_RATE'] <- gsub("[a-z,A-Z, ,/]", "", substring(grep('AC$CHROMATOGRAPHY: FLOW_RATE',record, value = TRUE, fixed = TRUE),30)), wantedmat[i,'FLOW_RATE'] <- FLOW_RATE_FALSE)
    ifelse(length(FLOW_RATE_TRUE)==1, wantedmat[i,'FLOW_RATE_UNIT'] <- gsub("[0-9, ,.]", "", substring(grep('AC$CHROMATOGRAPHY: FLOW_RATE',record, value = TRUE, fixed = TRUE),30)), wantedmat[i,'FLOW_RATE_UNIT'] <- FLOW_RATE_FALSE)
    ifelse(length(RETENTION_TIME_TRUE)==1, wantedmat[i,'RETENTION_TIME'] <- gsub("[a-z,A-Z, ]", "", substring(grep('AC$CHROMATOGRAPHY: RETENTION_TIME',record, value = TRUE, fixed = TRUE),35)), wantedmat[i,'RETENTION_TIME'] <- RETENTION_TIME_FALSE)
    ifelse(length(RETENTION_TIME_TRUE)==1, wantedmat[i,'RETENTION_TIME_UNIT'] <- gsub("[0-9, ,.]","", substring(grep('AC$CHROMATOGRAPHY: RETENTION_TIME',record, value = TRUE, fixed = TRUE),35)), wantedmat[i,'RETENTION_TIME_UNIT'] <- RETENTION_TIME_FALSE)
    ifelse(length(SOLVENT_A_TRUE)==1, wantedmat[i,'SOLVENT_A'] <- substring(grep('AC$CHROMATOGRAPHY: SOLVENT A',record, value = TRUE, fixed = TRUE),30), wantedmat[i,'SOLVENT_A'] <- SOLVENT_A_FALSE)
    ifelse(length(SOLVENT_B_TRUE)==1, wantedmat[i,'SOLVENT_B'] <- substring(grep('AC$CHROMATOGRAPHY: SOLVENT B',record, value = TRUE, fixed = TRUE),30), wantedmat[i,'SOLVENT_B'] <- SOLVENT_A_FALSE)
    
    ## The MS block
    ifelse(length(BASE_PEAK_TRUE)==1, wantedmat[i,'BASE_PEAK'] <- substring(grep('MS$FOCUSED_ION: BASE_PEAK ',record, value = TRUE, fixed = TRUE),26), wantedmat[i,'BASE_PEAK'] <- BASE_PEAK_FALSE)
    ifelse(length(PRECURSOR_MZ_TRUE)==1, wantedmat[i,'PRECURSOR_MZ'] <- substring(grep('MS$FOCUSED_ION: PRECURSOR_M/Z ',record, value = TRUE, fixed = TRUE),30), wantedmat[i,'PRECURSOR_MZ'] <- PRECURSOR_MZ_FALSE)
    ifelse(length(PRECURSOR_TYPE_TRUE)==1, wantedmat[i,'PRECURSOR_TYPE'] <- substring(grep('MS$FOCUSED_ION: PRECURSOR_TYPE ',record, value = TRUE, fixed = TRUE),31), wantedmat[i,'PRECURSOR_TYPE'] <- PRECURSOR_TYPE_FALSE)
    ifelse(length(DATA_PROCESSING_TRUE)==1, wantedmat[i,'DATA_PROCESSING'] <- DATA_PROCESSING_TRUE, wantedmat[i,'DATA_PROCESSING'] <- DATA_PROCESSING_FALSE)
    
    
    # The Peak block
    ifelse(length(SPLASH_TRUE)==1, wantedmat[i,'SPLASH'] <- substring(grep('PK$SPLASH:',record, value = TRUE, fixed = TRUE),12), wantedmat[i,'SPLASH'] <- SPLASH_FALSE)
    
    
    # wantedmat[i,'SPLASH'] <- substring(grep('PK$SPLASH:',record, value = TRUE, fixed = TRUE),12)
    
    ## The deep links
    wantedmat[i,'EULINK'] <- paste("http://massbank.eu/MassBank/jsp/FwdRecord.jsp?id=", substring(grep('ACCESSION:',record, value = TRUE, fixed = TRUE),12), sep="")
    wantedmat[i,'JPLINK'] <- paste("http://www.massbank.jp/jsp/FwdRecord.jsp?id=", substring(grep('ACCESSION:',record, value = TRUE, fixed = TRUE),12), sep="")
    if(verbose.output) print(wantedmat[i,c('ACCESSION','SPLASH')])
  }
  
  ## Write the csv file
  write.csv(wantedmat,csvname)
  return("Successfully wrote the csv")
}













# Run
preprocessContributorToMassBankWorkflow(folder, accessionPrefix, xlsxFile, ms_type, takeRecordedNames)
runContributorToMassBankWorkflow(folder, applyIntensityThreshold, reprocess = FALSE)

# Utility
#cleanContribtorDirectoryForReprocessing(folder)
#correctRecords(folder)

#aggregateSpectra_all(parentFolder, aggregationFolder)
#aggregateSpectra_ready(processedFolders, aggregationFolder, tag = "Validated")
