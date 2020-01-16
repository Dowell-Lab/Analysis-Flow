#!/bin/bash
# Pausing Ratio Calculator
# Author: Zachary Maas <zama8258@colorado.edu>
# TODO
# -	Combine all modes into 1 script using an argument
#	- Use	tempfiles instead	of a fixed directory
# - Add better error handling
#	-	Make it all faster by	pipelining stuff better...
# - Strip out	redundant	TSS's
# -	Update all logging messages to include the BEDFILE
#	-	Update usage message...

# Strict error checking
set -e
# set -o nounset
set -o errexit
#	Assume only	utf-8
export LC_ALL=C
# Echo start time

# Argument Parsing
# -pus/-pds = pause upstream/downstream
# -gds = gene upstream
# FIXME - add back in	-gus with new modes

usage()
{
		echo "calc_pausing_indices_polya.sh - a script for calculating fixed-width pausing indices"
		echo "Example:"
		echo "    ./calc_pausing_indices_fixwin.sh --pus=-100 --pds=300 --gds=2000 --bedfile=example.bedgraph"
		echo "Usage:"
		echo "    -h/--help -- Display this help message."
		echo "    --ref     -- Maximal Isoform FPKM Bedfile to Use"
		echo "    --pus     -- Pausing bases upstream"
		echo "    --pds     -- Pausing bases downstream"
		echo "    --gds     -- Gene bases downstream from polyA"
		echo "    --bedfile -- The filename of the bedfile to use"
		exit 0
}

while [ "$1" != "" ]; do
		PARAM=$(echo "$1" | awk -F= '{print $1}')
		VALUE=$(echo "$1" | awk -F= '{print $2}')
		case $PARAM in
				-h | --help)
						usage
						exit
						;;
				--ref)
						ref=$VALUE
						;;
				--pus)
						pus=$VALUE
						;;
				--pds)
						pds=$VALUE
						;;
				--gds)
						gds=$VALUE
						;;
				--bedfile)
						bedfile=$VALUE
						;;
				--outdir)
						outdir=$VALUE
						;;
				*)
						echo "ERROR: unknown parameter \"$PARAM\""
						usage
						exit 1
						;;
		esac
		shift
done

# Set Gene Downstream. This will need to change once we add support
# for multiple modes of operation in a single script.
gus=$(echo "$pds+1" | bc)

echo "[PARAM] Up: $pus, Down: $pds, Gene Up: $gus, Gene Down: $gds, BEDFILE: $bedfile, Outdir: $outdir"
echo "[WARNING] Gene Down Parameter will not be used."

# Make sure we have the necessary modules on the cluster
if ! type -t bedtools 
then module load bedtools
fi

################################################################################
################################################################################

# Variables we always need
InterestFile="$bedfile" # The location of the bedfiles
Infile=$ref
baseBed=$(basename "$bedfile")
OutFile="$outdir"/"$baseBed"_pause_ratios_"$gds".data

#	During debugging, we write out all output to disk so that we can
#	examine it and see what's going on with our script changes. This is
#	not necessary during production.
testing=false

if $testing; then
		# Variables	-	DEBUG
		echo "[LOG] Running ""$bedfile"" in Debug Mode"
		TmpDir=charli_pi

		OutGenePosFile=$DirPrefix/$TmpDir/"$bedfile"_"$gds"_pos_tss.bed
		OutBodyPosFile=$DirPrefix/$TmpDir/"$bedfile"_"$gds"_pos_body.bed
		OutGeneNegFile=$DirPrefix/$TmpDir/"$bedfile"_"$gds"_neg_tss.bed
		OutBodyNegFile=$DirPrefix/$TmpDir/"$bedfile"_"$gds"_neg_body.bed

		InterestFilePos=$DirPrefix/$TmpDir/"$bedfile"_"$gds"_interest_pos.bed
		InterestFileNeg=$DirPrefix/$TmpDir/"$bedfile"_"$gds"_interest_neg.bed

		GeneOutPos=$DirPrefix/$TmpDir/"$bedfile"_"$gds"_out_gene_pos.bed
		GeneOutNeg=$DirPrefix/$TmpDir/"$bedfile"_"$gds"_out_gene_neg.bed
		BodyOutPos=$DirPrefix/$TmpDir/"$bedfile"_"$gds"_out_body_pos.bed
		BodyOutNeg=$DirPrefix/$TmpDir/"$bedfile"_"$gds"_out_body_neg.bed

		FinalPos=$DirPrefix/$TmpDir/"$bedfile"_"$gds"_out_final_pos.bed
		FinalNeg=$DirPrefix/$TmpDir/"$bedfile"_"$gds"_out_final_neg.bed

else

		echo "[LOG] Running ""$bedfile"" in Production Mode."
		TmpDir=$(mktemp -d)

		OutGenePosFile="$TmpDir""/""$(uuidgen)"
		OutBodyPosFile="$TmpDir""/""$(uuidgen)"
		OutGeneNegFile="$TmpDir""/""$(uuidgen)"
		OutBodyNegFile="$TmpDir""/""$(uuidgen)"

		InterestFilePos="$TmpDir""/""$(uuidgen)"
		InterestFileNeg="$TmpDir""/""$(uuidgen)"

		GeneOutPos="$TmpDir""/""$(uuidgen)"
		GeneOutNeg="$TmpDir""/""$(uuidgen)"

		BodyOutPos="$TmpDir""/""$(uuidgen)"
		BodyOutNeg="$TmpDir""/""$(uuidgen)"

		FinalPos="$TmpDir""/""$(uuidgen)"
		FinalNeg="$TmpDir""/""$(uuidgen)"

		# Clean up temp files on exit
		function cleanup {
				rm -rf "$TmpDir"
				echo "[LOG] Deleted temporary directory $TmpDir"
		}
		# Register the cleanup function to be called on the EXIT signal
		trap cleanup EXIT

fi

################################################################################
################################################################################

#	First we pre-filter refseq according to the input parameters, so
#	that we can use bedtools to do the necessary genome arithmetic
#	later. At the same time, we split our datafile into a stranded
#	format to accomplish the same thing.

# TODO This step is the prime location for breaking apart this thing
# to have mode support. I think the easiest way to do this might just
# be to break apart the awk commands for each respective mode and use
# variable expansion inside other variable expansion to replace the
# commands on the fly. A sort of bizarre metaprogramming, if you will.
echo "[LOG] Prefiltering ""$bedfile"
awk -v OFS='\t' -v pus="$pus" -v pds="$pds" \
		'{if ($6 == "+") print $1, $2-pus, $2+pds, $4, $5, $6}' "$Infile"\
		| sort -k1,1 -k2,2n | \
		awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4}' \
				> "$OutBodyPosFile" &
awk -v OFS='\t' -v gus="$gus" -v gds="$gds" \
		'{if ($6 == "+") print $1, $2+gus, $3, $4, $5, $6}' "$Infile" \
		| sort -k1,1 -k2,2n | \
		awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4}' \
				> "$OutGenePosFile" &
awk -v OFS='\t' -v pus="$pus" -v pds="$pds" \
		'{if ($6 == "-") print $1, $3-pus, $3+pds, $4, $5, $6}' "$Infile" \
		| sort -k1,1 -k2,2n | \
		awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4}' \
				> "$OutBodyNegFile" &
awk -v OFS='\t' -v gus="$gus" -v gds="$gds" \
		'{if ($6 == "-") print $1, $3-gus, $3, $4, $5, $6}' "$Infile" \
		| sort -k1,1 -k2,2n | \
		awk -v OFS='\t' '{if ($2 < $3) print $1, $2, $3, $4}' \
				> "$OutGeneNegFile" &
echo "[LOG] Splitting data file ""$bedfile" &
awk -v OFS='\t' '{if ($4 > 0) print $1, $2, $3, $4}' "$InterestFile" \
		| sort -k1,1 -k2,2n > "$InterestFilePos" &
awk -v OFS='\t' '{if ($4 < 0) print $1, $2, $3, $4}' "$InterestFile" \
		| sort -k1,1 -k2,2n > "$InterestFileNeg" &
wait

# With the first portion of the analysis done, we proceed to calculate
# the sum of reads in the regions we gathered using awk in the above
# procedure. This is the slowest part of the script, based on current
# testing. This should be fine (FIXME, maybe), because we already
# consider the strandedness of the data in the preliminary step of
# separating strands and calculating a modified reference sequence.

echo "[LOG] Calculating Region Sums ""$bedfile"
bedtools map -a "$OutGenePosFile" -b "$InterestFilePos" -c 4 -o sum \
		| awk -F '\t' '($5 != "." && $5 != 0) ' > "$GeneOutPos" &
bedtools map -a "$OutBodyPosFile" -b "$InterestFilePos" -c 4 -o sum \
		| awk -F '\t' '($5 != "." && $5 != 0) ' > "$BodyOutPos" &
bedtools map -a "$OutGeneNegFile" -b "$InterestFileNeg" -c 4 -o sum \
		| awk -F '\t' '($5 != "." && $5 != 0) ' > "$GeneOutNeg" &
bedtools map -a "$OutBodyNegFile" -b "$InterestFileNeg" -c 4 -o sum \
		| awk -F '\t' '($5 != "." && $5 != 0) ' > "$BodyOutNeg" &
wait

# Here we calculate the total length of the region for each gene,
# creating a sum to be used later for normalizing by gene-length in
# pausing methods that vary the length (by doing things like going all
# the way to the polyA site).

# With counts for all genes calculated, we can proceed to calculate
# coverage for every gene that we haven't thrown out (we drop genes
# that lack any reads in the paused region or the gene-body region,
# since those missing gene-body reads leads to division by zero). Here
# we calculate gene read coverage normalizing by length (TODO). Using
# process substitution, we immediately pipe those output values into a
# final file for our last step.

# Finish Coverage Statistics by dividing by length.
echo "[LOG] Calculating Coverage Statistics ""$bedfile"
paste "$BodyOutPos" \
			<(awk -F '\t' 'NR==FNR{a[NR]=$5;next}{print $5+a[FNR]}' \
						"$GeneOutPos" "$BodyOutPos") \
			<(awk -F '\t' 'NR==FNR{a[NR]=$3-$2;next}{print ($3-$2)+a[FNR]}' \
						"$GeneOutPos" "$BodyOutPos") | \
		awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $4, $5, $6/$7}' > "$FinalPos" &
paste "$BodyOutNeg" \
			<(awk -F '\t' 'NR==FNR{a[NR]=$5;next}{print $5+a[FNR]}' \
						"$GeneOutNeg" "$BodyOutNeg") \
			<(awk -F '\t' 'NR==FNR{a[NR]=$3-$2;next}{print ($3-$2)+a[FNR]}' \
						"$GeneOutNeg" "$BodyOutNeg") | \
			awk -F '\t' -v OFS='\t' '{print $1, $2, $3, $4, $5, -($6/$7)}' > "$FinalNeg" &
wait

# This is the last step, and the most complicated awk procedure. We
# use a associative array with the (FIXME?) gene name as the key. Then
# we can calculate the final pausing index while also retaining our
# normalized coverage statistics.

# FIXME figure out a way to	account for	refseq strandedness earlier on...
echo "[LOG] Calculating Pausing Index ""$bedfile"
awk -F '\t' 'FNR==NR{a[$4]=$5; next} ($4 in a) {print $4,"+",$5/a[$4], $6}' \
		"$GeneOutPos" "$FinalPos" | sort -t$'\t' -k1,1 > "$OutFile"
awk -F '\t' 'FNR==NR{a[$4]=$5; next} ($4 in a) {print $4,"-",$5/a[$4], $6}' \
		"$GeneOutNeg" "$FinalNeg" | sort -t$'\t' -k1,1 >> "$OutFile"
echo "[LOG] Done $bedfile $gds"
