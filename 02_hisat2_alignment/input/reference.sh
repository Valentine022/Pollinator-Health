# download .gtf and .dna.toplevel.fa files from EnsemblMetozoa
# set up hisat2 in /bin which is linked in PATH
# https://bioinformatics-core-shared-training.github.io/RNAseq_September_2018/Supplementary_Materials/S1_Getting_raw_reads_from_SRA.html
# useful manual for alignning

# create directories
mkdir input
mkdir results

mkdir hisat2
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip
rm hisat2-2.1.0-Linux_x86_64.zip
export PATH=$PATH:"/data/home/bt18612/bin"

# index reference
hisat2_extract_splice_sites.py -v Bombus_terrestris.Bter_1.0.40.gtf > splice_sites.txt
hisat2_extract_exons.py -v Bombus_terrestris.Bter_1.0.40.gtf > exons.txt
hisat2-build Bombus_terrestris.Bter_1.0.dna.toplevel.fa --ss splice_sites.txt --exon exons.txt Bter.hisat2

# copy files into directory
cp ../../*.fastq .

# align files job script
#!/bin/sh
#$ -cwd           # Set the working directory for the job to the current directory
#$ -pe smp 1      # Request 1 core
#$ -l h_rt=24:0:0 # Request 24 hour runtime
#$ -l h_vmem=1G   # Request 1GB RAM

module load samtools
module load parallel
export PATH="$PATH:/data/home/bt18612/bin"

cat names.txt | parallel -t "hisat2 -x Bter.hisat2 -U {}.fastq > ../results/{}.sam"
for file in ../results/*.sam; do samtools sort $file; done
for file in ../results/*.sam; do samtools index $file; done
