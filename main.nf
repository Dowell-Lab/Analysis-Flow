#!/usr/bin/env nextflow

// TODO Change this to a user set parameter
params.dataDir = "/scratch/Users/zama8258/nina_rnaseq_processed"
params.bamDir = params.dataDir + "/mapped/bams/"
params.bedgraphsInit = params.dataDir + "/mapped/bedgraphs/*.bedGraph"
// These should be parameterized in user config
params.refseq = "/scratch/Shares/dowell/genomes/hg38/hg38_refseq.bed"
params.conversionFile = "FIXME"
params.pauseUpstream=-30
params.pauseDownstream=300
params.pauseTag="FIXME"
params.outdir="FIXME"

samples = Channel
.fromPath(params.bedgraphsInit)
.map { file -> tuple(file.baseName, file, params.bamDir + file.baseName + ".sorted.bam")}

process filterIsoform {
  validExitStatus 0
  tag "$name"
  publishDir "${params.outdir}/isoform/", mode: 'copy', pattern: "*.sorted.isoform_max.bed"
	input:
		set val(prefix), file(bedGraph), file(bam) from samples

	output:
		set val(prefix), file(bedGraph), file("*.sorted.isoform_max.bed") into filteredIsoforms

	script: 
	"""
	export InterestFile=${bam}
	export RefSeq=${params.refseq}
	export ConversionFile=${params.conversionFile}
	calc_maximal_isoform.bash
	"""
}

process calcPauseIndices {
  validExitStatus 0
  tag "$name"
  publishDir "${params.outdir}/pausing/", mode: 'copy', pattern: "*.data"
	input:
		set val(prefix), file(bedGraph), file(isoformMax) from filteredIsoforms

	output:
		set val(prefix), file("*.data") into pauseIndices

	script: 
	"""
	calculate_pause_index_to_polya.sh \
	  --ref=${isoformMax} \
		--pus=${params.pauseUpstream} \
		--pds=${params.pauseUpstream} \
		--gds=${params.pauseTag} \
		--outdir=$PWD \
		--bedfile=${bedGraph}
	"""
}
