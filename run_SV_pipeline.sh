#!/bin/bash -l 
#SBATCH --job-name=fix_manta
#SBATCH --account=pawsey0848
#SBATCH --partition=work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mail-user=gavin.monahan@pawsey.org.au
#SBATCH --mail-type=END
#SBATCH --error=%j.%x.err
#SBATCH --output=%j.%x.out

#Load singularity and nextflow modules. The specific version numbers may change over time. 
#Check what modules are available with `module spider <tool_name>`

module load singularity/3.8.6-nompi
module load nextflow/22.04.3
conda activate graphviz

#Run the pipeline
nextflow run /software/projects/pawsey0848/sv/Germline-StructuralV-nf/main.nf --input /scratch/pawsey0848/sv/batch1/samples.tsv --ref /software/projects/pawsey0848/sv/references/hg38_masked/Homo_sapiens_assembly38_masked.fasta --intervals /software/projects/pawsey0848/sv/references/hg38_regions.bed.gz -config /software/projects/pawsey0848/sv/Germline-StructuralV-nf/config/setonix.config --annotsv /software/projects/pawsey0848/sv/references/ --minSVsize 39
