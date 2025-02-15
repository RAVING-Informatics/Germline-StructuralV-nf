#!/bin/bash

#PBS -P er01
#PBS -N pullContainers
#PBS -l walltime=02:00:00
#PBS -l ncpus=1
#PBS -l mem=190GB
#PBS -W umask=022
#PBS -q copyq
#PBS -l wd
#PBS -l storage=scratch/er01+gdata/er01

module load singularity

# TODO create lint test to check this script and nextflow.config have same containers

# specify singularity cache dir (consistent with gadi.config)
export SINGULARITY_CACHE_DIR=/scratch/$PROJECT/$(whoami)/singularity

# pull containers
singularity pull docker://quay.io/biocontainers/bcftools:1.15.1--hfe4b78e_1
singularity pull docker://quay.io/biocontainers/mulled-v2-40295ae41112676b05b649e513fe7000675e9b84:0b4be2c719f99f44df34be7b447b287bb7f86e01-0 
singularity pull docker://brentp/smoove:v0.2.7
singularity pull docker://quay.io/biocontainers/survivor:1.0.7--hd03093a_2
singularity pull docker://sydneyinformaticshub/tiddit:3.5.2
singularity pull docker://sydneyinformaticshub/annotsv:3.2.1
singularity pull docker://quay.io/biocontainers/manta:1.6.0--h9ee0642_2
