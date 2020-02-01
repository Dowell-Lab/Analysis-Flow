#!/usr/bin/env nextflow

// TODO Change this to a user set parameter
params.conditionsTable = "/scratch/Users/zama8258/analysis_flow/test/conditions_table.txt"
params.designTable = "/scratch/Users/zama8258/analysis_flow/test/design.txt"
params.dataDir = "/scratch/Users/zama8258/nina_rnaseq_processed"
params.bamDir = params.dataDir + "/mapped/bams/"
params.bedgraphs = params.dataDir + "/mapped/bedgraphs/*.bedGraph"
// These should be parameterized in user config
params.refseq = "/scratch/Shares/dowell/genomes/hg38/hg38_refseq.bed"
params.conversionFile = "/scratch/Users/zama8258/pause_analysis_src/refseq_to_common_id.txt"
params.metageneNumRegions=100
params.pauseUpstream=-30
params.pauseDownstream=300
params.pauseTag="pipeline"
params.outdir="FIXME"

// Parse the main conditions table for later use
condTable = Channel
.fromPath(params.conditionsTable)
.splitText() { tuple(it.split()[0], it.split()[1])}
.groupTuple(by: 1)
.map() { [(it[1]): (it[0])] }
.reduce { a, b -> a+b }
.view() {"COND: $it"}

designTable = Channel
.fromPath(params.designTable)
.splitText() { tuple(it.split()[0], it.split()[1])}
.into() { designTableForDESeq; designTableForMetagene }

// We generate bam names from bedgraph names to ensure that they match.
samples = Channel
.fromPath(params.bedgraphs)
.map { file -> tuple(file.baseName, file, "$params.bamDir" + "$file.baseName" + ".sorted.bam")}

// Step 1 -- Filter for maximal isoforms
process filterIsoform {
	cpus 4
	memory '16 GB'
	time '1h'
  tag "$prefix"
  publishDir "${params.outdir}/isoform/", mode: 'copy', pattern: "*.sorted.isoform_max.bed", overwrite: true
	input:
		set val(prefix), file(bedGraph), val(bam) from samples

	output:
		set val(prefix), file(bedGraph), file("*.sorted.isoform_max.bed") into filteredIsoforms

	module 'python/3.6.3'
	module 'bedtools'
	module 'subread'
	script: 
	"""
	export NUM_CORES=4
	export InterestFile=${bam}
	export RefSeq=${params.refseq}
	export ConversionFile=${params.conversionFile}
	calc_maximal_isoform.bash
	"""
}

filteredIsoforms
.into() { filteredIsoformsForSingleRef; filteredIsoformsForPausing }

// Use a seeded random number to take the same sample every time.
singleRef = filteredIsoformsForSingleRef
.randomSample(1, 3432)
.first()

allBam = Channel
.fromPath("$params.bamDir" + "*.bam")
.toSortedList()

// Step 2.1 -- Generate counts for DESeq2
process countsForDESeq {
	cpus 8
	memory '16 GB'
	time '30m'
  tag "$prefix"
  publishDir "${params.outdir}/counts/", mode: 'copy', pattern: "counts*.txt*", overwrite: true
	input:
		set val(prefix), file(bedGraph), val(isoform_max) from singleRef
		file(bam) from allBam

	output:
		file("counts*.txt") into countsTableForDESeq
		file("counts*.txt") into countsTableForPCA

	module 'python/3.6.3'
	module 'bedtools'
	module 'subread'
	script:
	"""
	# Export array as a function to get it into the script.
	function exportArray {
	  BamFiles=(${bam})
	}
	export -f exportArray
	export InFile=${isoform_max}
	bash -c \"exportArray; . counts_for_deseq.bash\"
	"""
}

// Step 2.2 -- Run Differential Expression Analysis
process runDESeq {
	cpus 1
	memory '4 GB'
	time '5m'
  tag "$prefix"
  publishDir "${params.outdir}/deseq/", mode: 'copy', pattern: "*", overwrite: true
	input:
		val(condition_dict) from condTable
		each condition from designTableForDESeq
		each file(counts) from countsTableForDESeq

	output:
	file("*") into deSeqOutput

	script:
	"""
	tail -n +2 ${counts} > ${counts}_fix
	run_deseq.r \
				-c ${counts}_fix \
				-gi \$(echo '${condition_dict.(condition[0])}' | sed -E -e 's/\\[|\\]|,//g') \
				-gj \$(echo '${condition_dict.(condition[1])}' | sed -E -e 's/\\[|\\]|,//g') \
				-ni ${condition[0]} \
				-nj ${condition[1]} \
				-t ${params.conversionFile} \
	"""
}

// Step 2.3 -- Run Principal Component Analysis
process runPCA {
	cpus 1
	memory '4 GB'
	time '2m'
  tag "$prefix"
  publishDir "${params.outdir}/pca/", mode: 'copy', pattern: "*", overwrite: true
	input:
		file(counts) from countsTableForPCA

	output:
	  file("*") into pcaOutput

	script:
	"""
	tail -n +2 ${counts} > ${counts}_fix
	run_pca.r -c ${counts}_fix
	"""
}

// Step 3.1 -- Generate counts for metagene analysis
process countsForMetagene {
	cpus 8
	memory '16 GB'
	time '30m'
  tag "$prefix"
  publishDir "${params.outdir}/counts/", mode: 'copy', pattern: "metagene_counts*.txt*", overwrite: true
	input:
		set val(prefix), file(bedGraph), val(isoform_max) from singleRef
		file(bam) from allBam

	output:
		file("metagene_counts_sense_fix.txt") into countsTableMetageneSense
		file("metagene_counts_antisense_fix.txt") into countsTableMetageneAntiSense

	module 'python/3.6.3'
	module 'bedtools'
	module 'subread'
	script:
	"""
	# Export array as a function to get it into the script.
	function exportArray {
	  BamFiles=(${bam})
	}
	export -f exportArray
	export InFile=${isoform_max}
	export numRegions=${params.metageneNumRegions}
	bash -c \"exportArray; . counts_for_metagene.bash\"
	"""
}

// Step 3.2 -- Run metagene analysis
process runMetagene {
	cpus 1
	memory '4 GB'
	time '5m'
  tag "$prefix"
  publishDir "${params.outdir}/metagene/", mode: 'copy', pattern: "*", overwrite: true
	input:
		file(countsSense) from countsTableMetageneSense
		file(countsAntiSense) from countsTableMetageneAntiSense
		val(condition_dict) from condTable
		each condition from designTableForMetagene

	output:
		file("*") into metageneOutput

	script:
	"""
	metagene_efficient.r \
				-s ${countsSense} \
				-a ${countsAntiSense} \
				-gi \$(echo '${condition_dict.(condition[0])}' | sed -E -e 's/\\[|\\]|,//g') \
				-gj \$(echo '${condition_dict.(condition[1])}' | sed -E -e 's/\\[|\\]|,//g') \
				-ni ${condition[0]} \
				-nj ${condition[1]} \
				-n ${params.metageneNumRegions} \
				-o ${condition[0]}_vs_${condition[1]}_metagene.pdf
	"""
}

// Step 4.1 -- Calculate pause index values
process calcPauseIndices {
	cpus 4
	memory '8 GB'
	time '10m'
  validExitStatus 0
  tag "$name"
  publishDir "${params.outdir}/pausing/", mode: 'copy', pattern: "*.data", overwrite: true
	input:
		set val(prefix), file(bedGraph), file(isoformMax) from filteredIsoformsForPausing

	output:
		set val(prefix), file("*.data") into pauseIndices

	script:
	"""
	calculate_pause_index_to_polya.sh \
	  --ref=${isoformMax} \
		--pus=${params.pauseUpstream} \
		--pds=${params.pauseUpstream} \
		--gds=${params.pauseTag} \
		--outdir=. \
		--bedfile=${bedGraph}
	"""
}

// Step 4.2 -- Generate pause index figures
