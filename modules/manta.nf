// run manta structural variant detection and convert inversions
process manta {
	debug false
	publishDir "${params.outDir}/${batchName}", mode: 'copy'
	container "${params.mulled__container}"

	input:
	val(bams)
	val(batchName)
	path(ref)
	path(ref_fai)

	output:
	tuple val(batchName), path("manta/Manta_${batchName}.candidateSmallIndels.vcf.gz")		, emit: manta_small_indels
	tuple val(batchName), path("manta/Manta_${batchName}.candidateSmallIndels.vcf.gz.tbi")	, emit: manta_small_indels_tbi
	tuple val(batchName), path("manta/Manta_${batchName}.candidateSV.vcf.gz")					, emit: manta_candidate
	tuple val(batchName), path("manta/Manta_${batchName}.candidateSV.vcf.gz.tbi")				, emit: manta_candidate_tbi
	tuple val(batchName), path("manta/Manta_${batchName}.diploidSV.vcf.gz")					, emit: manta_diploid
	tuple val(batchName), path("manta/Manta_${batchName}.diploidSV.vcf.gz.tbi")				, emit: manta_diploid_tbi
	tuple val(batchName), path("manta/Manta_${batchName}.diploidSV_converted.vcf.gz")			, emit: manta_diploid_convert
	tuple val(batchName), path("manta/Manta_${batchName}.diploidSV_converted.vcf.gz.tbi")		, emit: manta_diploid_convert_tbi

	script:
	// TODO: add optional parameters $args. 
	def intervals = params.intervals ? "--callRegions $params.intervals" : ''
	def input_files = bams.collect{"--bam ${it}"}.join(' ')
	"""
	echo $input_files
	# configure manta SV analysis workflow
	configManta.py \
		${input_files} \
		--referenceFasta ${params.ref} \
		--runDir manta \
		$intervals

	# run SV detection 
	manta/runWorkflow.py -m local -j ${task.cpus}

	# clean up outputs
	mv manta/results/variants/candidateSmallIndels.vcf.gz \
		manta/Manta_${batchName}.candidateSmallIndels.vcf.gz
	mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
		manta/Manta_${batchName}.candidateSmallIndels.vcf.gz.tbi
	mv manta/results/variants/candidateSV.vcf.gz \
		manta/Manta_${batchName}.candidateSV.vcf.gz
	mv manta/results/variants/candidateSV.vcf.gz.tbi \
		manta/Manta_${batchName}.candidateSV.vcf.gz.tbi
	mv manta/results/variants/diploidSV.vcf.gz \
		manta/Manta_${batchName}.diploidSV.vcf.gz
	mv manta/results/variants/diploidSV.vcf.gz.tbi \
		manta/Manta_${batchName}.diploidSV.vcf.gz.tbi
	
	# convert multiline inversion BNDs from manta vcf to single line
	convertInversion.py \$(which samtools) ${params.ref} \
		manta/Manta_${batchName}.diploidSV.vcf.gz \
		> manta/Manta_${batchName}.diploidSV_converted.vcf

	# zip and index converted vcf
	bgzip manta/Manta_${batchName}.diploidSV_converted.vcf
	tabix manta/Manta_${batchName}.diploidSV_converted.vcf.gz
	"""
} 

// rehead manta SV vcf for merging 
process rehead_manta {
	debug true
	publishDir "${params.outDir}/${batchName}/manta", mode: 'copy'
	container "${params.bcftools__container}"

	input:
	val(ids)
	tuple val(batchName), path(manta_diploid_convert)
	tuple val(batchName), path(manta_diploid_convert_tbi)

	output:
	path("Manta_*.vcf")	, emit: manta_VCF
		
	script:
	"""
	# create new header for merged vcf
	def idlist = [ids.collect{"${it}"}.join('_manta\n')]
	printf "${idlist}" > ${batchName}_rehead_manta.txt

	# replace batchName with caller_sample for merging
	bcftools reheader \
		Manta_${batchName}.diploidSV_converted.vcf.gz \
		-s ${batchName}_rehead_manta.txt \
		-o Manta_${batchName}.vcf.gz

	# gunzip vcf
	gunzip Manta_${batchName}.vcf.gz
	"""
}
