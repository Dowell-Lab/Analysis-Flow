#!/bin/bash
# gen_featurecounts.sbatch --- Generate Featurecounts
#
# Filename: gen_featurecounts.sbatch
# Description: Make Featurcounts with Coordinate Modification
# Author: Student Zachary Maas <zama8258@colorado.edu>
# Maintainer: Student Zachary Maas <zama8258@colorado.edu>
#

# Commentary:
#
# This file contains code for generating counts tables using
# featureCounts, with builtin support for coordinate modification for
# modifying regions in the bedfile.
#

# Code:

# -e - exits on the first error
# -u - fails when unset variables are called
# -o pipefail - fails when error in a pipe
set -euo pipefail

# Recently tweaked for better POSIX compliance.
function logr() {
    echo "[""$(date -d@$SECONDS -u +%H:%M:%S)""]: $*"
}

# Dynamically load bedtools if it isn't available
if ! type -t bedtools; then
		module load bedtools
fi
if ! type -t featureCounts; then
		module load subread
fi

# Set the number of available cores
NUM_CORES=8

# TmpDir=/scratch/Users/zama8258/processed_nascent/scratch/features
BaseDir="$PWD"
outFull="$BaseDir"/counts_"$1".txt
safFull="$BaseDir"/full.saf

## Generate the SAF File
awk -v OFS='\t' '{print $4, $1, $2, $3, $6}' "$InFile" > "$safFull"

## Change directory because FeatureCounts is picky about running in
## the same directory as the bams.
logr "Intersecting the Reference Sequence"

# Finally, run featurecounts to actually do the counting.
logr "Performing Featurecounts (Full Gene)"
featureCounts \
		-T "$NUM_CORES" \
		-s 1 \
		-F 'SAF' \
		-a "$safFull" \
		-o "$outFull" \
		${BamFiles[@]}

tail -n +2 "$outFull" > "$outFull"_without_header

logr "Done"

#
# gen_featurecounts.sbatch ends here
