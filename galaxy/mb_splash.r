#!/usr/bin/env Rscript

# ---------- Preparations ----------
# Load libraries
library(RMassBank)

# Setup R error handling to go to stderr
options(show.error.messages=F, error=function() { cat(geterrmessage(), file=stderr()); q("no",1,F) } )

# Set encoding
#options(encoding="UTF-8")

# Set options
options(stringAsfactors=FALSE, useFancyQuotes=FALSE)

# Take in trailing command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
	print("Error! No or not enough arguments given.")
	print("Usage: $0 mb_record_input.txt mb_record_output.txt")
	quit(save="no", status=1, runLast=FALSE)
}

# Data files
mb_input <- args[1]
mb_output <- args[2]



# ---------- MAIN ----------
# Read input file
con <- file(args[1], "r", blocking=FALSE)
record <- readLines(con=con)
close(con)

# Test whether splash exists
splash_index <- which(grepl(x=record, pattern="^PK\\$SPLASH"))

# Splash does not exist in input file
if (length(splash_index) == 0) {
	# Insert PK$SPLASH just before PK$ANNOTATION
	annotation_index <- which(grepl(x=record, pattern="^PK\\$ANNOTATION"))
	if (length(annotation_index) == 0) {
		stop("Error! Input file has no annotated peak list (PK$ANNOTATION).")
		quit(save=no, status=3)
	}

	# Insert PK$SPLASH
	record <- c(record[c(1:annotation_index-1)], "PK$SPLASH: 0", record[c(annotation_index:length(record))])
	splash_index <- annotation_index
}

# Get peak list
peak_beg <- which(grepl(x=record, pattern="^PK\\$PEAK"))
peak_end <- which(grepl(x=record, pattern="^//"))

if ((length(peak_beg) == 0) | (length(peak_end) == 0)) {
	stop("Error! Input file has no peak list (PK$PEAK).")
	quit(save=no, status=4)
}

peak_list <- record[c((peak_beg+1):(peak_end-1))]

# Make data frame
peak_mz  <- as.numeric(lapply(X=strsplit(x=trimws(peak_list), split=' '), FUN=function(x) { x <- x[1] })) 
peak_int <- as.numeric(lapply(X=strsplit(x=trimws(peak_list), split=' '), FUN=function(x) { x <- x[2] }))
peaks <- data.frame("mz" = peak_mz, "int"=peak_int)

# Get splash
splash <- RMassBank:::getSplash(peaks=peaks)

# Write splash into record
record[splash_index] <- paste0("PK$SPLASH: ", splash)

# Write output file
con <- file(args[2], "w", blocking=FALSE)
writeLines(text=record, con=con)
close(con)

