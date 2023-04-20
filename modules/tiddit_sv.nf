// run tiddit structural variant detection
process tiddit_sv {
	debug false
	publishDir "${params.outDir}/${sampleID}/tiddit", mode: 'copy'
	container "${params.tiddit__container}"

	input:
	tuple val(sampleID), file(bam), file(bai)
	path(ref)
	path(ref_fai)

	output:
	tuple val(sampleID), path("Tiddit_${sampleID}_PASSsv.vcf")	, emit: tiddit_vcf
	tuple val(sampleID), path("${sampleID}_sv.ploidies.tab")	, emit: tiddit_ploidy
	tuple val(sampleID), path("${sampleID}_sv_tiddit")		, emit: tiddit_workdir
	
	script:
	// TODO: add optional parameters $args. 
	"""
	tiddit \
		--sv \
		-q 20 \
		--bam ${bam} \
		--ref ${params.ref} \
		-o ${sampleID}_sv \
		--threads ${task.cpus}

	# rename vcf to show its from tiddit 
	mv ${sampleID}_sv.vcf \
		Tiddit_${sampleID}_sv.vcf

	# filter to pass only variants 
	grep -E "#|PASS" Tiddit_${sampleID}_sv.vcf \
		> Tiddit_${sampleID}_PASSsv.vcf
	"""
}

// rehead tiddit SV vcf for merging 
process rehead_tiddit {
	debug false 
	publishDir "${params.outDir}/${sampleID}/tiddit", mode: 'copy'
	container "${params.bcftools__container}"

	input:
	tuple val(sampleID), path(tiddit_vcf)
		
	output:
	tuple val(sampleID), path("Tiddit_${sampleID}_final.vcf")	, emit: tiddit_VCF
		
	script:
	"""
	# bgzip and index tiddit vcf 
	bgzip Tiddit_${sampleID}_PASSsv.vcf
	tabix Tiddit_${sampleID}_PASSsv.vcf.gz

	# create new header for merged vcf
	printf "${sampleID}_tiddit\n" > ${sampleID}_rehead_tiddit.txt

	# replace sampleID with caller_sample for merging
	bcftools reheader \
		Tiddit_${sampleID}_PASSsv.vcf.gz \
		-s ${sampleID}_rehead_tiddit.txt \
		-o Tiddit_${sampleID}_final.vcf.gz
	
	# gunzip vcf
	gunzip Tiddit_${sampleID}_final.vcf.gz
	
	#clean up
	rm -r ${sampleID}_rehead_tiddit.txt
	"""
}

// Merge all TIDDIT SVs
process merge_tiddit {
	debug false
	publishDir "${params.outDir}/tidditsv", mode: 'copy'
	container "${params.svdb__container}"
	
	input:
	path(tiddit_VCF)

	output:
	tuple val(params.batchName), path("tiddit_pass_merged.vcf")		, emit: tiddit_merged_vcf

	script:
	def input_vcfs = tiddit_VCF.collect{"${it}"}.join(' ')
	"""
	echo ${input_vcfs}
	svdb --merge --vcf ${input_vcfs} --bnd_distance 300 --overlap 0.8 > tiddit_pass_merged.vcf
	"""
}
