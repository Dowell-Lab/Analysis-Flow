#!/bin/bash

Rscript metagene_efficient.r \
	-s /home/zach/dowell_lab/pausing_meta_analysis/data/metagene_counts_sense_fix.txt \
	-a /home/zach/dowell_lab/pausing_meta_analysis/data/metagene_counts_antisense_fix.txt \
	-gi "PO_1_S1_R1_001.sorted.bam", "PO_2_S2_R1_001.sorted.bam" \
	-gj "C413_1_S3_R1_001.sorted.bam", "C413_2_S4_R1_001.sorted.bam" \
	-ni "Control" \
	-nj "Treatment" \
	-n 100 \
	-o /home/zach/metagene_test_out.pdf
