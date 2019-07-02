## Job script to check raw RNA seq files in fastq (loop) ######


#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 4      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=4G   # Request 1GB RAM

module load fastqc
for f in *.gz; do fastqc $f; done
