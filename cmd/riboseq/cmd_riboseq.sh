#!/bin/bash

module purge
module load nextflow/21.04.1
module load singularity/3.5.3

nextflow -C nextflow_riboseq.config run riboseq.nf

