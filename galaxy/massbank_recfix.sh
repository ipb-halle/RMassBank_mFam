#!/bin/sh
WORKDIR="/tmp/msp2massbank/converted to MassBank/"

cd "${WORKDIR}"

# Removes lines with duplicate peaks
for i in *; do
       perl -i -ne 'print if ! $a{$_}++' $i
done

# Fix PK$NUM_PEAK
for i in *; do
	NUM_PEAK=$(cat $i | perl -pe 's/\n/\|/' | perl -pe 's/.*PEAK//' | perl -pe 's/.*^(.+?)\|//' | perl -pe 's/\|\/\/.*/\|/' | perl -pe 's/\|/\n/g' | wc -l)
	NUM_PEAK=($NUM_PEAK)
	perl -i -pe "s/PK\\\$NUM_PEAK.*/PK\\\$NUM_PEAK: ${NUM_PEAK}/" $i
done

