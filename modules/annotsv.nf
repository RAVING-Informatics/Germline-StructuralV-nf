// run annotSV for variant annotation (human, mouse only)
process annotsv {
	debug false
	publishDir "${params.outDir}/${params.batchName}/annotsv", mode: 'copy'
	container "${params.annotsv__container}"

	input:
<<<<<<< Updated upstream
	tuple val(sampleID), path(mergedVCF)
	path annotsvDir
	val annotsvMode
=======
	tuple val(batchName), path(mergedVCF)
	path("${params.annotsv}")
>>>>>>> Stashed changes

	output:
	tuple val(batchName), path("*_AnnotSV")

	script: 
	// Apply annotation mode flag to command
	def mode = params.annotsvMode
	
	// Change output file name based on annotation mode
	def outputFile = null
	    if (mode == 'full') {
               outputFile = "${sampleID}_full_AnnotSV.tsv"
            } else if (mode == 'split') {
               outputFile = "${sampleID}_split_AnnotSV.tsv"
            } else if (mode == 'both') {
               outputFile = "${sampleID}_both_AnnotSV.tsv"
            } else {
               throw new RuntimeException("Invalid option for --annotSV: ${mode}")}
	
	//Pass any additional flags to the AnnotSV 
	def extraArgs = params.extraAnnotsvFlags ?: ''
	"""
	AnnotSV \
<<<<<<< Updated upstream
		-SVinputFile ${sampleID}_merged.vcf \
		-annotationsDir ${params.annotsvDir} \
=======
		-SVinputFile ${batchName}_merged.vcf \
		-annotationsDir ${params.annotsv} \
>>>>>>> Stashed changes
		-bedtools bedtools -bcftools bcftools \
<<<<<<< HEAD
		-annotationMode both \
		-genomeBuild GRCh38 \
		-includeCI 1 \
		-overwrite 1 \
		-outputFile ${batchName}_AnnotSV.tsv \
		-SVminSize ${params.minSVsize}
=======
		-annotationMode ${mode} \
		-genomeBuild GRCh38 \
		-includeCI 1 \
		-overwrite 1 \
		-outputFile ${outputFile} ${extraArgs}
>>>>>>> upstream/main
	"""
}
