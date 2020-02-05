#!/usr/bin/env nextflow
/*
========================================================================================
                         AnalysisFlow -- Automated Analysis Pipeline
========================================================================================
 Automated Analysis Pipeline. Started 2020-01-14
 #### Homepage / Documentation
 https://github.com/Dowell-Lab/Analysis-Flow
 #### Authors
 Zachary Maas <zama8258@colorado.edu>
========================================================================================
========================================================================================
*/

def helpMessage() {
	log.info"""
	=========================================
	AnalysisFlow v${params.version}
  =========================================
		This pipeline requires that you manually create a configuration file to
		analyze your data, due to the large number of parameters involved.

	Usage:
  The typical command for running the pipeline is as follows:
		nextflow run main.nf -profile example
		Required arguments:
		-profile                      Configuration profile to use. <base, slurm, example>
		--workdir                     Nextflow working directory where all intermediate files are saved.
		""".stripIndent()
	}

params.help = false
if (params.help){
    helpMessage()
    exit 0
}

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
.map { file -> tuple(file.baseName, file, "$params.cramDir" + "$file.baseName" + ".sorted.cram")}

// Step 0 -- Convert cram to bam
process cramToBam {
	cpus 16
	memory '40 GB'
	time '2h'
	tag "$prefix"

	input:
		set val(prefix), file(bedGraph), val(cram) from samples

	output:
		set val(prefix), file(bedGraph), file("${prefix}.sorted.bam") into bamSamples
		file("${prefix}.sorted.bam") into bamInit

	module 'samtools'
	script:
	"""
	samtools view -@ 16 -b -1 -T ${params.reffasta} ${cram} > ${prefix}.sorted.bam
	"""
}

allBam = bamInit
.toSortedList()

// Step 1 -- Filter for maximal isoforms
process filterIsoform {
	cpus 4
	memory '16 GB'
	time '1h'
  tag "$prefix"
  publishDir "${params.outdir}/isoform/", mode: 'copy', pattern: "*.sorted.isoform_max.bed", overwrite: true
	input:
		set val(prefix), file(bedGraph), val(bam) from bamSamples

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
				-t ${params.conversionFile}
	rm -f Rplots.pdf
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
	rm -f Rplots.pdf
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
	rm -f Rplots.pdf
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
// TODO -- The existing pausing code needs some work...

// Step 5.1 -- Collect figures and put them into the same folder for report generation
/*
process generateReport {
	cpus 1
	memory '4 GB'
	time '5m'
  tag "$name"
  publishDir "${params.outdir}/report/", mode: 'copy', pattern: "*.zip", overwrite: true
  publishDir "${params.outdir}/report/", mode: 'copy', pattern: "*.pdf", overwrite: true
	input:
		file(metagene) from metageneOutput
	  file(pca) from pcaOutput
		file(deseq) from deSeqOutput

	output:
		file("analysis_figures.zip") into generatedArchive
		//file("analysis_report.pdf") into generatedReport

	script:
	"""
	echo \"Unimplemented\"
	zip analysis_figures.zip ${metagene} ${pca} ${deseq}
	"""
}
*/
