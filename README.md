
# Pipeline for the compilations of MassBank records from raw data

## Information to be supplied

### A meta data spreadsheet containing information for the measured substances

#### Sheet *Compounds*
Sheet for the meta data of the sample, compound, spectrum, and authors. The table contains one row for each spectrum.

- **File** MSP file containing the spectrum; located under path `raw data/exported as raw msp` in the project folder (mandatory, e.g. Rutin.msp)
- **Name** Preferred compound name (mandatory, e.g. Rutin)
- **Synonyms** Common compound name synonyms (optional, e.g. Rutinoside; Eldrin)
- **Structure / Unique ID** The compound structure - one of the following three options (mandatory, e.g. 5280805)
  - **InChI** 
  - **SMILES** 
  - **PubChem CID** 
- **Formula** Molecular formula of the compound (optional, e.g. C27H30O16)
- **Exact mass** The exact mass of the compound (optional, e.g. 610.15339)
- **CHROMATOGRAPHY** The name of the chromatography instrumentation (mandatory, e.g. Symmetry C18 Column, Waters)
- **RT (min)** The retention time of the parent ion (mandatory, e.g. 9.5)
- **INSTRUMENT** The name of the mass spectrometry instrumentation (mandatory, e.g.  LTQ Orbitrap XL, Thermo Scientfic; HP-1100 HPLC, Agilent)
- **INSTRUMENT_TYPE** The type of the mass spectrometry instrumentation (mandatory, e.g.  LC-ESI-ITFT)
- **IONIZATION** The kind of ionization (mandatory, e.g. ESI) (mandatory, e.g. ESI)
- **Collision energy** Collision energy for fragmentation (mandatory, e.g. 15eV)
- **Adduct** The parent ion species (mandatory, e.g. [M-H]-)
- **MS/MS Acquisition mode** The acquisition mode (mandatory, e.g. DDA)
- **Ionization mode** The ionization mode (mandatory, e.g. Negative)
- **Confidence** The confidence of compound identification (mandatory, Pure standard or Predicted)
- **Compound class** The compound class (optional, e.g. Flavonoid)
- **Found in** The sample origin (optional, e.g. Tomato leaf (Solanum habrochaites LA1777))
- **Database links** Links to compound databases (optional, e.g. KEGG C05625; ChEBI 28527)
- **Authors** All authors who substantially contributed to the spectrum (mandatory, e.g. Ales Svatos, Ravi Kumar Maddula, MPI for Chemical Ecology, Jena, Germany)
- **Publication** Relevant pblication related to the spectrum (optional, e.g. F. Rasche, A. Svatos, R.K. Maddula, C. Boettcher and S. Boecker. Computing fragmentation trees from tandem mass spectrometry data. Anal. Chem., 2011, 83, 1243-1251)

#### Sheet *Chromatography*

Sheet for properties describing the chromatography instrumentation. Composed of two columns, namely `Property` with the property names and `Value` with the value corresponding to the property. It is recommended to use MassBank tags as property names. Please find the recommended MassBank record tags [here](https://github.com/MassBank/MassBank-web/blob/master/Documentation/MassBankRecordFormat.md#246-acchromatography-subtag-description).

#### Sheet *Mass_Spectrometry*

Sheet for properties describing the chromatography instrumentation. Composed of two columns, namely `Property` with the property names and `Value` with the value corresponding to the property. It is recommended to use MassBank tags as property names. Please find the recommended MassBank record tags [here](https://github.com/MassBank/MassBank-web/blob/master/Documentation/MassBankRecordFormat.md#243-acmass_spectrometry-ms_type).

## Install

This pipeline is designed for a Linux machine.
In command line clone the slightly adapted code of the RMassBank package.
`git clone -b branch_treutler https://github.com/MassBank/RMassBank.git`
Install the package in RStudio.
`install.packages("[path to RMassBank source]", repos = NULL, type="source")`
Install openbabel.
`sudo apt-get install openbabel`

## Creating a project

A project is a self-contained data set for the generation of a set of MassBank records from raw data. Projects are currently stored at `/IPB/Projects/2017_005_MS-databases/mFam contributions/`. The project folder name is arbitrary and usually comprises the name and affiliation of a data contributor. The project folder structure includes the following folders.

- **converted to MassBank** target for the compiled MassBank records
- **converted to msp** target for the compiled MassBank records converted to msp format
- **meta data** meta data spreadsheet(s in different versions)
- **raw data** raw data section
- **raw data/converted to abf** raw data in abf format
- **raw data/exported as raw msp** raw msp files (exported from MS-DIAL)
- **raw data/raw** raw data in vendor format
- **script** arbitrary scripts or misc. supplemental information

The recommended workflow for the raw msp files is as follows.

- **Conversion** of raw data files in vendor format to abf format (see Abf converter web page [here](https://www.reifycs.com/AbfConverter/)
- **Raw data processing** using MS-DIAL (see MS-DIAL web page [here](http://prime.psc.riken.jp/Metabolomics_Software/MS-DIAL/)
- **Export** of measured MS/MS spectra from MS-DIAL to raw msp format (see settings `msp format` and `deconvoluted` in the alignment export dialog)

The meta data spreadsheet is filled independently by the experimentalist starting from the template in project `template` in folder `meta data`. 

## Running a project

### Prerequisites

- Follow the installation guidelines above
- meta data spreadsheet under path `meta data` in the project folder
- raw msp files under path `raw data/exported as raw msp` in the project folder

### Run

The conversion of raw data and meta data to MassBank records is done by RMassBank with some pre- and post-processing. To run RMassBank adapt and run the following code in the R console.
```
source("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam_code/RMassBank_mFam/CommunityToMassBank.R")

parentFolder <- "/mnt/data/IPB/Projects/2017_005_MS-databases/mFam contributions/"

ms_type <- "MS2"
accessionPrefix <- "[TODO]"
xlsxFile <- "[TODO]"
takeRecordedNames <- TRUE
applyIntensityThreshold <- FALSE

preprocessContributorToMassBankWorkflow(folder, accessionPrefix, xlsxFile, ms_type, takeRecordedNames)
runContributorToMassBankWorkflow(folder, applyIntensityThreshold, reprocess = FALSE)
```

The console will print the progress of the RMassBank workflow. The MassBank records are generated in folder path `converted to MassBank` in the project folder. For an example of an MassBank record see [here](https://massbank.eu/MassBank/jsp/RecordDisplay.jsp?id=EQ310401&dsn=Eawag).

There are some utility function for the management of projects as follows.
```
source("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam_code/RMassBank_mFam/CommunityToMassBank.R")

parentFolder <- "/mnt/data/IPB/Projects/2017_005_MS-databases/mFam contributions/"
aggregationFolder <- paste(parentFolder, "mFam Aggregation", sep = "")

## if a RMassBank run failes, this command cleans the project folder to the original state
#cleanContribtorDirectoryForReprocessing(folder)

## utility function for the correction of common record format problems
#correctRecords(folder)

## functions for the aggregation of spectra from multiple projects - either for all projects or a specified set of projects
#aggregateSpectra_all(parentFolder, aggregationFolder)
#aggregateSpectra_ready(processedFolders, aggregationFolder, tag = "Validated")
```

### Test in Docker 

Start the container with `docker run -v $PWD/data:/data -w /data -it --rm sneumann/msp2massbank
bash`

and run the example data:

```

mkdir /tmp/msp2massbank/
cp -avx /data/* /tmp/msp2massbank/

/usr/local/bin/msp2massbank.r "XY" test.xlsx

```
