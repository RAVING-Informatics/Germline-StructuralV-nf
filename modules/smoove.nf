// run smoove structural variant detection
process smoove {
	debug false
	publishDir "${params.outDir}/${batchName}", mode: 'copy'
	container "${params.smoove__container}"

	input:
	//	tuple val(sampleID), file(bam), file(bai)
	val(bams)
	val(batchName)
	path(ref)
	path(ref_fai)

	output:
	tuple val(batchName), path("smoove/${batchName}-smoove.genotyped.vcf.gz")		, emit: smoove_geno
	tuple val(batchName), path("smoove/${batchName}-smoove.genotyped.vcf.gz.csi")	, emit: smoove_geno_csi
	tuple val(batchName), path("smoove/${batchName}.split.bam")			, emit: smoove_split, optional: true
	tuple val(batchName), path("smoove/${batchName}.split.bam.csi")			, emit: smoove_split_csi, optional: true
	tuple val(batchName), path("smoove/${batchName}.disc.bam")			, emit: smoove_disc, optional: true
	tuple val(batchName), path("smoove/${batchName}.disc.bam.csi")			, emit: smoove_disc_csi, optional: true
	tuple val(batchName), path("smoove/${batchName}.histo")				, emit: smoove_histo, optional: true
	
	script:
	// Create sym links to bams
	// links = file(bams).mklink("smoove/bams/")
	// TODO: add optional parameters $args. 
	
	def input_files = bams.collect{"${it}"}.join(' ')	
	"""
	smoove call -d --name ${batchName} \
		--fasta ${params.ref} \
		--outdir smoove \
		--processes ${task.cpus} \
		--genotype ${input_files} 
	"""
} 

// rehead smoove genotyped vcf for merging 
process rehead_smoove {
	debug false 
	publishDir "${params.outDir}/${batchName}/smoove", mode: 'copy'
	container "${params.bcftools__container}" 

	input:
	val(ids)
	tuple val(batchName), path(smoove_geno)
		
	output:
	path("Smoove_${batchName}.vcf")	, emit: smoove_VCF	
		
	script:
	"""
	# create new header for merged vcf
	idlist = [ids.collect{"${it}"}.join('_smoove\n')]
	printf "${idlist}" > ${batchName}_rehead_smoove.txt

	# replace batchName with caller_sample for merging
	bcftools reheader \
		${batchName}-smoove.genotyped.vcf.gz \
		-s ${batchName}_rehead_smoove.txt \
		-o Smoove_${batchName}.vcf.gz
	
	# gunzip vcf
	gunzip Smoove_${batchName}.vcf.gz
	
	#clean up
	#rm -r ${batchName}_rehead_smoove.txt
	"""
}
