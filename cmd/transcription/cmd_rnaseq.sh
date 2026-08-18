#!/bin/bash

module purge
module load singularity/3.5.3
module load nextflow/21.04.1

nextflow -C nextflow_RNAseq.config run RNAseq.nf -resume
