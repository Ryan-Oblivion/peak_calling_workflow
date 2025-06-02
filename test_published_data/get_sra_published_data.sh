#!/usr/bin/bash

#SBATCH --partition=hpc_a10_a
#SBATCH --mem=80GB
#SBATCH --time=6:00:00
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=2


source ~/.bashrc_rj_test.sh

conda activate sra_toolkit_rj

# the sra data came from this study https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239447

# making a directory to get the output dir
mkdir './sra_data'
cd sra_data

prefetch SRR18503706 SRR18503662 SRR18503641 -O '.'

fasterq-dump SRR18503706 SRR18503662 SRR18503641 -O '.'
#fasterq-dump SRR18503662 -O '.'
#fasterq-dump 

# changing the names of the fastq files
mv SRR18503706_1.fastq control_H3k27me3_r1_1.fastq
mv SRR18503706_2.fastq control_H3k27me3_r1_2.fastq

mv SRR18503662_1.fastq control_H3k27me3_r2_1.fastq
mv SRR18503662_2.fastq control_H3k27me3_r2_2.fastq

mv SRR18503641_1.fastq control_H3k27me3_r3_1.fastq
mv SRR18503641_2.fastq control_H3k27me3_r3_2.fastq

gzip *.fastq
