######## Use fastq-dump to download files - has to be set up first by adding to PATH ######

#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 4      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=4G   # Request 1GB RAMERR883767

./fastq-dump --gzip ERR883768
./fastq-dump --gzip ERR883769
./fastq-dump --gzip ERR883770
./fastq-dump --gzip ERR883771
./fastq-dump --gzip ERR883772
./fastq-dump --gzip ERR883773
./fastq-dump --gzip ERR883774
./fastq-dump --gzip ERR883775
./fastq-dump --gzip ERR883776
./fastq-dump --gzip ERR883777
./fastq-dump --gzip ERR883778
./fastq-dump --gzip ERR883779
./fastq-dump --gzip ERR883780
./fastq-dump --gzip ERR883781
./fastq-dump --gzip ERR883782
./fastq-dump --gzip ERR883783
./fastq-dump --gzip ERR883784
./fastq-dump --gzip ERR883785
./fastq-dump --gzip ERR883786
./fastq-dump --gzip ERR883787
./fastq-dump --gzip ERR883788
./fastq-dump --gzip ERR883789
./fastq-dump --gzip ERR883790
./fastq-dump --gzip ERR883791
./fastq-dump --gzip ERR883792
./fastq-dump --gzip ERR883793
