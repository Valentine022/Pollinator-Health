#!/bin/bash
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 2      # Request 1 core
#$ -l h_rt=36:0:0 # Request 36 hour runtime
#$ -l h_vmem=8G   # Request 16GB RAM

module load R
Rscript Apis.R
