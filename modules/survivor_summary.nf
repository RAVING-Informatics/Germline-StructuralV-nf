// generate summary counts for merged VCF
process survivor_summary {
	debug false
	publishDir "${params.outDir}/${params.batchName}/survivor", mode: 'copy'
	container "${params.survivor__container}"

	input:
	tuple val(batchName), path(mergedVCF)

	output:
	tuple val(batchName), path("*")
	
	script:
	"""
	SURVIVOR vcftobed ${batchName}_merged.vcf \
		0 -1 \
		${batchName}_merged.bed
	
	SURVIVOR stats ${batchName}_merged.vcf \
		-1 -1 -1 \
		${batchName}_merged.stats.txt
	"""

}

process survivor_venn {
	debug false
	publishDir "${params.outDir}/${batchName}/survivor", mode: 'copy'
	container "${params.mulled__container}"

	input:
	tuple val(batchName), path(mergedVCF)

	output:
	tuple val(batchName), path("*")

	script:
	"""
	"""
}
