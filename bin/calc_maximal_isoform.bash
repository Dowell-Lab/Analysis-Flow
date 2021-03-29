#!/bin/bash
## This script contains code to filter out the maximal isoform of all
## genes in a bedfile using bedops bedmap.

set -euxo pipefail

## Start timing the script
SECONDS=0

## This is a logrging function that will display the elapsed time of
## the script and how long the script has been executing for on the
## server. Useful for tracking execution time in various scripts.
function logr {
    echo "[""$(date -d@$SECONDS -u +%H:%M:%S)""]: $*"
}

# Make sure we have the necessary modules on the cluster
if ! type -t bedtools
then module load bedtools
fi
if ! type -t featureCounts
then module load subread
fi

##############################
## Variables we always need ##
##############################

## This is a scratch directory that we use to hold temporary files.
## This is kept for legacy reasons, from before the script was used in
## a pipeline environment.
ScratchDir="$PWD"

## This is the output directory
OutDir="$ScratchDir"

#####################################
## Generate	Temporary Files in RAM ##
#####################################

## Here we generate the temporary directory used to store the temp
## files we use during file processing. This is a fast process because
## we're only reading and writing to RAM, which gives us IO on the
## order of 60Gb/s
logr "Running ""$InterestFile"" in Production Mode"
TmpDir=.

PosRefFile="$TmpDir"/pos_ref.saf
NegRefFile="$TmpDir"/neg_ref.saf

RootName=$(basename "$InterestFile" .bam)

# Cache the sum files	because they take the longest to generate
PosStrandSums="$ScratchDir"/"$RootName".pos.count
NegStrandSums="$ScratchDir"/"$RootName".neg.count

PosCountsMapped=0
NegCountsMapped=0

PosScalingFactor=0
NegScalingFactor=0

PosStrandMillionsNormalizedSums="$TmpDir"/pos_mil_norm_sums
NegStrandMillionsNormalizedSums="$TmpDir"/neg_mil_norm_sums

PosStrandFpkmNormalizedSums="$TmpDir"/pos_fpkm_norm_sums
NegStrandFpkmNormalizedSums="$TmpDir"/neg_fpkm_norm_sums

PosStrandFinalSorted="$ScratchDir"/"$RootName".pos.maxiso
NegStrandFinalSorted="$ScratchDir"/"$RootName".neg.maxiso

StrandsMerged="$ScratchDir"/"$RootName".merge.maxiso
FPKMCommonID="$ScratchDir"/"$RootName".merge.maxiso.common
FinalFPKM="$OutDir"/"$RootName".isoform_max
FinalOut="$FinalFPKM".bed

###############################
## Start Performing Analysis ##
###############################

## Here, we perform strand-specific filtering to account for `bedtools
## multicov` being unable to write out separate results for each
## strand. This strips down the RefSeq annotation to a bed6 format,
## and drops out any non NM_ coded genes, since for the purpose of
## this script we only care about protein coding genes.
logr "Filtering RefSeq Annotations by Strand for ""$InterestFile"	&
<"$RefSeq" awk -v OFS='\t' '{if ($6 == "+") print $4, $1, $2, $3, "+"}'\
 > "$PosRefFile" &
<"$RefSeq" awk -v OFS='\t' '{if ($6 == "-") print $4, $1, $2, $3, "-"}'\
 > "$NegRefFile" &
wait

## Here, we make sure that our bam file is indexed and if not, perform
## the indexing
if [ ! -f "$InterestFile".bai ]; then
		logr "No Index File found. Generating One."
		module load samtools
		samtools index -@ "$NUM_CORES" "$InterestFile"
fi

if [ ! -f "$PosStrandSums" ] || [ ! -f "$NegStrandSums" ]; then
		## Here, we construct a bedfile from our input BAM, in order to
		## normalize it by FPKM. We do this using bedtools multicov on the
		## sample in a strand specific manner. Note that bedtools multicov
		## will not account for duplicate reads by default, so if we want RPKM
		## instead of FPKM we need to add	the	`-D` flag
		logr "Calculating Strand-Specific Coverage for BAM"
		featureCounts \
				-T "$NUM_CORES" \
				-s "$1" \
				"$2" \
				-O \
				-F 'SAF' \
				-a "$PosRefFile" \
				-o "$PosStrandSums"_tmp \
				$InterestFile
		tail -n +3 "$PosStrandSums"_tmp | \
				awk -v OFS='\t' '{print $2, $3, $4, $1, 0, $5, $7}' | \
				grep -v ';' > "$PosStrandSums"
		featureCounts \
				-T "$NUM_CORES" \
				-s "$1" \
				"$2" \
				-O \
				-F 'SAF' \
				-a "$NegRefFile" \
				-o "$NegStrandSums"_tmp \
				$InterestFile
		tail -n +3 "$NegStrandSums"_tmp | \
				awk -v OFS='\t' '{print $2, $3, $4, $1, 0, $5, $7}' | \
				grep -v ';' > "$NegStrandSums"
fi


## Next, we need to find the sum of all reads mapped in our regions of
## interest, for FPKM. We do this by adding up every mapped count in
## the 5th column of our bedfile from bedtools multicov
read -r PosCountsMapped <<< "$(awk	-F '\t' -v OFS='\t'	'{sum += $7} END {print sum}' "$PosStrandSums")"
read -r NegCountsMapped <<< "$(awk	-F '\t' -v OFS='\t'	'{sum += $7} END {print sum}' "$NegStrandSums")"

## Then, we generate a scaling factor for RPKM by dividing our mapped
## counts by 1 million
read -r PosScalingFactor <<< "$((PosCountsMapped / 1000000))"
read -r NegScalingFactor <<< "$((NegCountsMapped / 1000000))"

## We proceed to divide every count in our bedfile by the
## strand-specific scaling factor
awk	-F '\t' -v OFS='\t' -v norm="$PosScalingFactor" \
		'{print $1, $2, $3, $4, $5, $6, $7 / norm}' "$PosStrandSums" > "$PosStrandMillionsNormalizedSums"
awk	-F '\t' -v OFS='\t' -v norm="$NegScalingFactor" \
		'{print $1, $2, $3, $4, $5, $6, $7 / norm}' "$NegStrandSums" > "$NegStrandMillionsNormalizedSums"

## We finally take the sums (now normalized for millions mapped), then
## divide those values by gene length.
awk -F '\t' -v OFS='\t'	'{print $1, $2, $3, $4, $5, $6, $7 / ($3 - $2)}' \
		"$PosStrandMillionsNormalizedSums" >	"$PosStrandFpkmNormalizedSums"
awk -F '\t' -v OFS='\t'	'{print $1, $2, $3, $4, $5, $6, $7 / ($3 - $2)}' \
		"$NegStrandMillionsNormalizedSums" >	"$NegStrandFpkmNormalizedSums"

## With the heavy lifting of normalization done, we can now easily
## take the maximum FPKM value for each isoform.
logr	"Re-Sorting Files for Merging"
sort -rnk7 "$PosStrandFpkmNormalizedSums" | sort -u -k4 | sort -k 1,1 -k2,2n > "$PosStrandFinalSorted"
sort -rnk7 "$NegStrandFpkmNormalizedSums" | sort -u -k4 | sort -k 1,1 -k2,2n > "$NegStrandFinalSorted"
wait

## Merge the two strands together. No sorting	is necessary yet, so save the bandwidth.
cat	"$PosStrandFinalSorted"	"$NegStrandFinalSorted" | sed -e 's/\.[0-9]\t/\t/' > "$StrandsMerged"

## Next, we use a supplementary python script to add common gene ID's
## so that we can filter out isoforms from our final file.
convert_isoform.py -l "$ConversionFile" \
									 -f "$StrandsMerged" \
									 -o "$FPKMCommonID"

## At long last, we sort in descending order and filter
sort -rnk7 "$FPKMCommonID" | sort -usk8 |
		sort -k1,1 -k2,2 > "$FinalFPKM"

awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6}' "$FinalFPKM" > "$FinalOut"

logr "Final Genes: ""$(wc -l "$FinalOut" | awk '{print $1}')"

## If we made it this far, we are done
logr "Done executing"
