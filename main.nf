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
params.pauseUpstream=-30
params.pauseDownstream=300
params.pauseTag="FIXME"
params.outdir="FIXME"

// Parse the main conditions table for later use
condTable = Channel
.fromPath(params.conditionsTable)
.splitText() { tuple(it.split()[0], it.split()[1])}
.groupTuple(by: 1)
.flatMap() { [(it[1]): (it[0])] }
.toList()
.subscribe() { println it }

designTable = Channel
.fromPath(params.designTable)
.splitText() { tuple(it.split()[0], it.split()[1])}
.subscribe() { println it }

// We generate bam names from bedgraph names to ensure that they match.
samples = Channel
.fromPath(params.bedgraphs)
.map { file -> tuple(file.baseName, file, "$params.bamDir" + "$file.baseName" + ".sorted.bam")}

process filterIsoform {
	cpus 8
	memory '16 GB'
	time '1h'
  tag "$prefix"
  publishDir "${params.outdir}/isoform/", mode: 'copy', pattern: "*.sorted.isoform_max.bed"
	input:
		set val(prefix), file(bedGraph), val(bam) from samples

	output:
		set val(prefix), file(bedGraph), file("*.sorted.isoform_max.bed") into filteredIsoforms

	script: 
	"""
	module load python/3.6.3
	export InterestFile=${bam}
	export RefSeq=${params.refseq}
	export ConversionFile=${params.conversionFile}
	calc_maximal_isoform.bash
	"""
}

singleRef = filteredIsoforms.first()

allBam = Channel
.fromPath("$params.bamDir" + "*.bam")
.toSortedList()

process countsForDeSeq {
	cpus 8
	memory '16 GB'
	time '1h'
  tag "$prefix"
  publishDir "${params.outdir}/counts/", mode: 'copy', pattern: "counts*.txt*"
	input:
		set val(prefix), file(bedGraph), val(isoform_max) from singleRef
		file(bam )from allBam

	output:
		file("counts*.txt") into countsTable

	script:
	"""
	module load python/3.6.3
	# Export array as a function to get it into the script.
	function exportArray {
	  BamFiles=(${bam})
	}
	export -f exportArray
	export InFile=${isoform_max}
	bash -c \"exportArray; . counts_for_deseq.bash\"
	"""
}

// process calcPauseIndices {
//   validExitStatus 0
//   tag "$name"
//   publishDir "${params.outdir}/pausing/", mode: 'copy', pattern: "*.data"
// 	input:
// 		set val(prefix), file(bedGraph), file(isoformMax) from filteredIsoforms
//
// 	output:
// 		set val(prefix), file("*.data") into pauseIndices
//
// 	script:
// 	"""
// 	calculate_pause_index_to_polya.sh \
// 	  --ref=${isoformMax} \
// 		--pus=${params.pauseUpstream} \
// 		--pds=${params.pauseUpstream} \
// 		--gds=${params.pauseTag} \
// 		--outdir=$PWD \
// 		--bedfile=${bedGraph}
// 	"""
// }
