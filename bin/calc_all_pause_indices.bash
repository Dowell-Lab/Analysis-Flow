#!/bin/bash
#SBATCH -p short
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=64gb
#SBATCH --mail-user=zama8258@colorado.edu
#SBATCH --output=/Users/zama8258/e_and_o/%x_%j.out
#SBATCH --error=/Users/zama8258/e_and_o/%x_%j.err

# Bed files
beddir=/Users/zama8258/iris_bed
t0="$beddir"/0-1_S11_L006_R1_001.bedGraph
t30_1="$beddir"/30-1_S15_L007_R1_001.bedGraph
t30_2="$beddir"/30_2_S3_R1_001.bedGraph
ca0_1="$beddir"/O-CA-1_S12_L006_R1_001.bedGraph
ca0_2="$beddir"/0_CA_2_S2_R1_001.bedGraph
ca30_1="$beddir"/30-CA-1_S16_L007_R1_001.bedGraph
ca30_2="$beddir"/30_CA_2_S4_R1_001.bedGraph

# FPKM files
fpkmdir=/Users/zama8258/iris_out/scratch
f_t0="$fpkmdir"/0-1_S11_L006_R1_001.sorted.isoform_max.bed
f_t30_1="$fpkmdir"/30-1_S15_L007_R1_001.sorted.isoform_max.bed
f_t30_2="$fpkmdir"/30_2_S3_R1_001.sorted.isoform_max.bed
f_ca0_1="$fpkmdir"/O-CA-1_S12_L006_R1_001.sorted.isoform_max.bed
f_ca0_2="$fpkmdir"/0_CA_2_S2_R1_001.sorted.isoform_max.bed
f_ca30_1="$fpkmdir"/30-CA-1_S16_L007_R1_001.sorted.isoform_max.bed
f_ca30_2="$fpkmdir"/30_CA_2_S4_R1_001.sorted.isoform_max.bed

outdir=/Users/zama8258/iris_out/pause
script=/scratch/Users/zama8258/pause_analysis_src/calculate_pause_index_to_polya.sh
upstream=-30
downstream=300
tag=IRIS

bash "$script" \
		 --ref="$f_t0" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$t0" &
bash "$script" \
		 --ref="$f_t30_1" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$t30_1" &
bash "$script" \
		 --ref="$f_t30_2" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$t30_2" &
bash "$script" \
		 --ref="$f_ca0_1" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$ca0_1" &
bash "$script" \
		 --ref="$f_ca0_2" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$ca0_2" &
bash "$script" \
		 --ref="$f_ca30_1" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$ca30_1" &
bash "$script" \
		 --ref="$f_ca30_2" \
		 --pus="$upstream" \
		 --pds="$downstream" \
		 --gds="$tag" \
		 --outdir="$outdir" \
		 --bedfile="$ca30_2" &
wait
