#!/bin/bash

INPUT="${1}"
OUTPUT="${2}"

# Copy file
cat ${INPUT} > ${OUTPUT}

# Remove lines with duplicate peaks
perl -i -ne 'print if ! $a{$_}++' ${OUTPUT}

# Fix PK$NUM_PEAK
NUM_PEAK=$(cat ${OUTPUT} | perl -pe 's/\n/\|/' | perl -pe 's/.*PEAK//' | perl -pe 's/.*^(.+?)\|//' | perl -pe 's/\|\/\/.*/\|/' | perl -pe 's/\|/\n/g' | wc -l)
NUM_PEAK=($NUM_PEAK)
perl -i -pe "s/PK\\\$NUM_PEAK.*/PK\\\$NUM_PEAK: ${NUM_PEAK}/" ${OUTPUT}

