#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import subworkflows to be run in the workflow
include { checkInputs       }   from './modules/check_cohort'
include { smoove            }   from './modules/smoove' 
include { rehead_smoove     }   from './modules/smoove'
include { manta             }   from './modules/manta'
include { rehead_manta      }   from './modules/manta'
include { tiddit_sv         }   from './modules/tiddit_sv'
include { rehead_tiddit     }   from './modules/tiddit_sv'
include { tiddit_cov        }   from './modules/tiddit_cov'
include { survivor_merge    }   from './modules/survivor_merge'
include { survivor_summary  }   from './modules/survivor_summary'
include { annotsv           }   from './modules/annotsv'

// Print the header to screen when running the pipeline
log.info """\

===================================================================
G E R M L I N E  S T R U C T U R A L  V - N F                   
===================================================================

Created by the Sydney Informatics Hub, University of Sydney

Documentation	@ https://github.com/Sydney-Informatics-Hub/Germline-StructuralV-nf

Cite		@ TODO:INSERT DOI

Log issues	@ https://github.com/Sydney-Informatics-Hub/Germline-StructuralV-nf/issues

===================================================================
Workflow run parameters 
===================================================================

version		: ${params.version}
input		: ${params.input}
reference	: ${params.ref}
annotsvDir	: ${params.annotsvDir}
annotsvMode	: ${params.annotsvMode}
outDir		: ${params.outDir}
workDir		: ${workflow.workDir}
AnnotSV flags	: ${params.extraAnnotsvFlags}

===================================================================
 """

// Help function 
// TODO: once finalised, add all optional and required flags in

def helpMessage() {
log.info"""

Usage:  nextflow run main.nf --input samplesheet.tsv --ref reference.fasta

Required Arguments:

	--input			Full path and name of sample input file (tsv format).

	--ref			Full path and name of reference genome (fasta format).

Optional Arguments:

	--outDir		Full path and name of results directory. 

	--intervals		Full path and name of the intervals file for Manta 
				(bed format).

	--annotsvDir		Full path to the directory housing the prepared
				Annotations_human directory for AnnotSV. 
	
	--minSVsize		Minimum SV size for AnnotSV to annotate (default = 50bp).

	--annotsvMode		Specify full, split, or both for AnnotSV output
				mode (default: both).
	
	--extraAnnotsvFlags	Additionally specify any valid AnnotSV flags. For 
				example: --extraAnnotsvFlags '-SVminSize 50 -vcf 1'. 

""".stripIndent()
}
 
workflow {

if ( params.help == true || params.ref == false || params.input == false ){
	// Invoke the help function above and exit
	helpMessage()
	exit 1

	} else {
	
	// Check inputs file exists
	checkInputs(Channel.fromPath(params.input, checkIfExists: true))
	
	// Split cohort file to collect info for each sample
	input = checkInputs.out
		.splitCsv(header: true, sep:"\t")
		.map { row -> tuple(row.sampleID, file(row.bam), file(row.bai))}

	// Call SVs with Manta  
	manta(input, params.ref, params.ref+'.fai')

	// Rehead manta vcf for merging 
	rehead_manta(manta.out.manta_diploid_convert, manta.out.manta_diploid_convert_tbi)

	// Call SVs with Smoove
	smoove(input, params.ref, params.ref+'.fai')

	// Rehead smoove vcf for merging  
	rehead_smoove(smoove.out.smoove_geno)

	// Run TIDDIT sv
	tiddit_sv(input, params.ref, params.ref+'.fai')
  
	// Rehead TIDDIT vcf for merging
	rehead_tiddit(tiddit_sv.out.tiddit_vcf)

	// Run TIDDIT cov 
	tiddit_cov(input, params.ref, params.ref+'.fai')

	// Collect VCFs for merging
	mergeFile = rehead_tiddit.out.tiddit_VCF
		.concat(rehead_smoove.out.smoove_VCF, rehead_manta.out.manta_VCF)
		.groupTuple() 

	// Run SURVIVOR merge
	survivor_merge(mergeFile)

	// Run SURVIVOR summary
	survivor_summary(survivor_merge.out.mergedVCF)

	// Run AnnotSV (optional)
	if (params.annotsvDir) {
		annotsv(survivor_merge.out.mergedVCF, params.annotsvDir, params.annotsvMode)}
	}}

workflow.onComplete {

summary = """
===================================================================
Workflow execution summary
===================================================================

Duration	: ${workflow.duration}
Success		: ${workflow.success}
workDir		: ${workflow.workDir}
Exit status	: ${workflow.exitStatus}
outDir		: ${params.outDir}

===================================================================

"""

println summary

}
