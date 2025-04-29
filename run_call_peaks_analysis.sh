#!/bin/env bash

#SBATCH --mem=20GB
#SBATCH --mail-type=FAIL,END
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --job-name=nextflow_chip
#SBATCH --partition=hpc_a10_a 


source /lustre/fs4/home/rjohnson/.bashrc_rj_test.sh
# source /ru-auth/local/home/rjohnson/.bashrc_rj_test.sh # or use this, should be the same thing 


conda activate nextflow_three

# plot_idr = 0.05 // the default is 0.05 to report in the png plots which peaks passed this threshold
# return_idr = 1  // the default is all peaks will be returned even if the report plots the ones that pass a certain number. 1 as default here will give back all peaks like idr has set

nextflow call_peaks_analysis_pipeline.nf -profile peak_calling_analysis -resume \
--plot_idr 0.1 \
--return_idr 0.1