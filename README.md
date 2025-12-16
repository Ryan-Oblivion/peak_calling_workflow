
# Please look at the README_documentation file




# peak_calling_workflow
# peak_calling_workflow

## commands used to get Hera's broadPeaks
"ls /lustre/fs4/risc_lab/scratch/hcanaj/K562_cutandtag_H1low_alignedbam_012925_hg38/all_aligned_bams_bw_hg38/Peak_Call_cutandtag_hg38/broad/H1low_H3k27me3_*[r1,r2,r3]*_001.trim.st.all.blft.qft.rmdup.sorted.bam_peaks.broadPeak | xargs cp -t ."

"ls /lustre/fs4/risc_lab/scratch/hcanaj/K562_cutandtag_H1low_alignedbam_012925_hg38/all_aligned_bams_bw_hg38/Peak_Call_cutandtag_hg38/broad/Scrm_H3k27me3_*[r1,r2,r3]*_001.trim.st.all.blft.qft.rmdup.sorted.bam_peaks.broadPeak| xargs cp -t ."

## path where I got Hera's bam files

"this is the path where i got hera's bam files from: /lustre/fs4/risc_lab/scratch/hcanaj/K562_cutandtag_H1low_alignedbam_012925_hg38/all_aligned_bams_bw_hg38"

## this is also where to find the duplicate information for each of the replicates
"/lustre/fs4/risc_lab/scratch/hcanaj/K562_cutandtag_H1low_alignedbam_012925_hg38/all_aligned_bams_bw_hg38"
then go into the directory for each replicate and find the duplicate log file "*dups.log"

# use this command to copy them to the current directory you are in
"ls /lustre/fs4/risc_lab/scratch/hcanaj/K562_cutandtag_{H1low,Scrm}_alignedbam_012925_hg38/*_H3k27me3_{r1,r2,r3}/*_dups.log | xargs cp -t ."
