#!/bin/bash

Rscript run_deseq.r \
				-c /home/zach/dowell_lab/pausing_meta_analysis/out/counts/counts_full.txt_without_header \
				-gi "C413_1_S3_R1_001.sorted.bam", "C413_2_S4_R1_001.sorted.bam" \
				-gj "PO_1_S1_R1_001.sorted.bam", "PO_2_S2_R1_001.sorted.bam" \
				-ni "Treatment" \
				-nj "Control" \
				-t /home/zach/dowell_lab/pausing_meta_analysis/src/refseq_to_common_id.txt \
				-o /home/zach/dowell_lab/analysis_flow/test/output
