// Merge manta, smoove, tiddit vcfs 
process survivor_merge {
	debug true
	publishDir "${params.outDir}/${params.batchName}/survivor", mode: 'copy'
	container "${params.survivor__container}"
		
	input:
	//tuple val(sampleID), path(mergelist)
	path(mergeFile)

	output:
	tuple val(params.batchName), path("${batchName}_merged.vcf"), emit: mergedVCF

	script:
	// TODO turn $mergeFile variable into input file 
	"""
	echo ${mergeFile} | xargs -n1 > ${params.batchName}_survivor.txt

	SURVIVOR merge ${params.batchName}_survivor.txt \
		1000 1 0 0 0 30 \
		${params.batchName}_merged.vcf
	"""
}
