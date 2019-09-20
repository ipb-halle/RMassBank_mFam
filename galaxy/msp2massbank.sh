#!/bin/sh
WORKDIR="/tmp/msp2massbank"
MSP="${1}"
MSP_NAME="${2}"
EXCEL="${3}"
EXCEL_NAME="${4}"

# Create working directory
mkdir -p "${WORKDIR}"

# Create subdirectories
mkdir -p "${WORKDIR}/converted to msp"
mkdir -p "${WORKDIR}/converted to MassBank"
mkdir -p "${WORKDIR}/meta data"
mkdir -p "${WORKDIR}/raw data"
mkdir -p "${WORKDIR}/raw data/exported as raw msp"
mkdir -p "${WORKDIR}/RMassBank"

cp "${MSP}" "${WORKDIR}/raw data/exported as raw msp/${MSP_NAME}"
cp "${EXCEL}" "${WORKDIR}/meta data/${EXCEL_NAME}"

