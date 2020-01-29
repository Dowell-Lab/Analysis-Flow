#!/bin/bash

set -euxo pipefail

#	Assume only	utf-8
export LC_ALL=C

module load python/3.6.3
module load subread

## Read from exported variable
regionFile="$InFile"

function logr {
    echo "[""$(date -d@$SECONDS -u +%H:%M:%S)""]: $*"
}

logr "Parsed Params: ""$regionFile"

logr "Starting Analysis"
numRegions=100

NUM_CORES=8
safFile=region_split.saf
countsSenseOut=metagene_counts_sense.txt
countsSenseFix=metagene_counts_sense_fix.txt
countsAntiSenseOut=metagene_counts_antisense.txt
countsAntiSenseFix=metagene_counts_antisense_fix.txt

logr "Generating Segmented SAF File"
bedgraph_split_for_metagene.py \
		-f "$regionFile" \
		-n "$numRegions" \
		-o "$safFile"

## DO NOT QUOTE $BamFiles. Featurecounts fails on processing the
## resulting quotes.
logr "Building Sense Counts Table from SAF"
featureCounts \
		-T "$NUM_CORES" \
		-s 1 \
		--fracOverlap 0.51 \
		-F 'SAF' \
		-a "$safFile" \
		-o "$countsSenseOut" \
		${BamFiles[@]}

logr "Building Antisense Counts Table from SAF"
featureCounts \
		-T "$NUM_CORES" \
		-s 2 \
		--fracOverlap 0.51 \
		-F 'SAF' \
		-a "$safFile" \
		-o "$countsAntiSenseOut" \
		${BamFiles[@]}

logr "Fixing Counts File"
tail -n+2 "$countsSenseOut" > "$countsSenseFix"
tail -n+2 "$countsAntiSenseOut" > "$countsAntiSenseFix"

logr "Done..."
