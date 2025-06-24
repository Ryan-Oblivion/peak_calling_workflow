nextflow.enable.dsl=2

// this is the path where i got hera's bam files from: /lustre/fs4/risc_lab/scratch/hcanaj/K562_cutandtag_H1low_alignedbam_012925_hg38/all_aligned_bams_bw_hg38

// the bams and bai files should have 3 replicates if you want it to work with idr
params.control_bams = 'bam_files/H1low_*{bam,bam.bai}'
params.other_control_bams = 'bam_files/dH1_*{bam,bam.bai}'
//params.published_bam = "bam_files/ENCFF915XIL_*{bam,bam.bai}"
//params.published_bam = "bam_files/ENCFF905CZD_*{bam,bam.bai}"

norm_control_bams_index_tuple_ch = Channel.fromFilePairs(params.control_bams)
other_control_bams_index_tuple_ch = Channel.fromFilePairs(params.other_control_bams)
//published_bam_index_tuple_ch = Channel.fromFilePairs(params.published_bam)

control_bams_index_tuple_ch = norm_control_bams_index_tuple_ch.concat(other_control_bams_index_tuple_ch )

//control_bams_index_tuple_ch.view()


params.wt_bams = 'bam_files/Scrm_*{bam,bam.bai}'
wt_bams_index_tuple_ch = Channel.fromFilePairs(params.wt_bams)


params.ref_genome = file('/lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/WholeGenomeFasta/genome.fa')
ref_genome_ch = Channel.value(params.ref_genome)

params.ref_genome_size = file('/lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/WholeGenomeFasta/genome.fa.fai')
ref_genome_size_ch = Channel.value(params.ref_genome_size)

params.gtf_file = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/RNA-seq/rep2/gencode.v38_ERCC.gtf')
gtf_ch = Channel.value(params.gtf_file)

// I need a list of gene regions. got this from irene but also included it in my motif analysis nextflow script
params.wtvs_lowup = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/RNA-seq/rep2/hg38-ERCC-UMI-alignment/DESeq2_results/WTvslowup-genebody.bed')
wtvslowup_genebody_ch = Channel.value(params.wtvs_lowup)

// including the nochange genes just to have something to compare to
params.wtvs_lowdown_nochange = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/RNA-seq/rep2/hg38-ERCC-UMI-alignment/DESeq2_results/WTvslowdown-basemeanmatchnochange-genebody.bed')
wtvslowdown_nochange_ch = Channel.value(params.wtvs_lowdown_nochange)


// I want a multiqc plot of the duplicate levels for the k27 data

params.dups_log = file('./dup_info/*_dups.log')
dups_log_ch = Channel.fromPath(params.dups_log)



params.meth_bed_file = file('/lustre/fs4/home/rjohnson/pipelines/hera_pipeline/results_SE/bl_filt_bed/bed_graphs_deeptools/control_CpG_r1r2r3_fp_filt_filt_coor_sorted_BL_filt_sort2_normalized_cpm.bed')
meth_bed_ch = Channel.value(params.meth_bed_file)


params.up_peaks_file = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/test_published_data/up_regulated_peaks.bed')
up_peaks_ch = Channel.value(params.up_peaks_file)

params.down_peaks_file = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/test_published_data/down_regulated_peaks.bed')
down_peaks_ch = Channel.value(params.down_peaks_file)

params.unchanging_peaks_file = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/test_published_data/unchanging_regulated_peaks.bed')
unchanging_peaks_ch = Channel.value(params.unchanging_peaks_file)

params.master_peaks_file = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/test_published_data/10kb_merged_masterPeak_geo_control_h1low.bed')
master_peaks_ch = Channel.value(params.master_peaks_file)

params.knownGene_bed_file = file('/lustre/fs4/home/rjohnson/downloads/genomes/hg38/genes/knownGene.hg38.bed')
knownGene_ch = Channel.value(params.knownGene_bed_file)


params.proseq_up_genes = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/proseq_data/proseq_up_genes_with_coord.tsv')
proseq_up_gene_ch = Channel.value(params.proseq_up_genes)

params.proseq_down_genes = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/proseq_data/proseq_down_genes_with_coord.tsv')
proseq_down_gene_ch = Channel.value(params.proseq_down_genes)

params.proseq_unchanging_genes = file('/lustre/fs4/home/rjohnson/pipelines/peak_calling_analysis_pipeline/proseq_data/proseq_unchanging_genes_with_coord.tsv')
proseq_unchanging_gene_ch = Channel.value(params.proseq_unchanging_genes)


// now adding the ATAC-seq bigwig signal files to nextflow

params.control_ATAC_bigwig = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/tracks/scr-allmerge.bw')
control_atac_bigwig_ch = Channel.value(params.control_ATAC_bigwig)

params.treatment_ATAC_bigwig = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/ATAC-seq/KA23/tracks/dH1-allmerge.bw')
treatment_atac_bigwig_ch = Channel.value(params.treatment_ATAC_bigwig)


// include {
//     make_alignment_bw_process_control

// }from './modules/peak_analysis_modules.nf'



include {
    mk_bw_call_peaks_workflow;
    plot_histone_data_workflow;
    plot_histone_calledpeak_workflow;
    plot_signal_up_down_peaks_workflow;
    plot_diff_peaks_over_diff_genes_workflow;
    plot_atac_signal_over_diff_peaks_workflow
    //get_diff_peaks_workflow

}from './workflows/call_peaks_workflow.nf'



include {

}from './workflows/find_diff_peaks_workflow.nf'

include {
    meth_enrichment_analysis_workflow



}from './workflows/methylation_enrichment_analysis_workflow.nf'




workflow {

    // now i will put the control bams and wt bams into the peak calling workflow
    // I will also put the reference genome in as the third entry

    mk_bw_call_peaks_workflow(control_bams_index_tuple_ch, wt_bams_index_tuple_ch, ref_genome_ch, ref_genome_size_ch, dups_log_ch )

    // get a channel with the final concat idr peaks
    //final_idr_concat_peaks_ch = mk_bw_call_peaks_workflow.out.concat_idr_peaks
    
    // take the emitted channels from the call peaks workflow
    control_bw_meta_ch = mk_bw_call_peaks_workflow.out.control_meta_bw_ch

    wt_bw_meta_ch = mk_bw_call_peaks_workflow.out.wt_meta_bw_ch

    // getting the cpm normalized bigwigs
    control_meta_cpm_bw = mk_bw_call_peaks_workflow.out.control_meta_cpm_bw_ch
    wt_meta_cpm_bw = mk_bw_call_peaks_workflow.out.wt_meta_cpm_bw_ch

    // the mk_bw_call_peaks_workflow workflow also emits the broad peaks that were called so I will grab those
    all_broadpeaks_ch = mk_bw_call_peaks_workflow.out.broadpeaks_ch

    // adding a workflow here to get the differential peaks

    //get_diff_peaks_workflow(all_broadpeaks_ch)


    // here I want to make a workflow that plots the chromatin features at a list of annotated sites (genes)
    // Also looking at the peaks around tss

    ////////// for testing do not run the plotting yet it takes too long //////////////////////////

    // NOT DOING THIS WORKFLOW ANYMORE. IT IS POINTLES TO SHOW THIS SIGNAL WITHOUT THE SCRM AND H1LOW ON THE SAME PLOT
    //plot_histone_data_workflow(control_bw_meta_ch, wt_bw_meta_ch, wtvslowup_genebody_ch, wtvslowdown_nochange_ch)

    // I will just make another workflow for plotting each histone-rep bigwig with its corresponding histone-rep broadPeak

    // NOT DOING THIS EITHER BECASUE THE SINGAL OVER THE LOOSE PARAMETER BROADPEAKS FROM MACS IS NOT USEFUL INFO BEFORE THE IDR AND MERGING FILTERING 
    //plot_histone_calledpeak_workflow(control_bw_meta_ch, wt_bw_meta_ch, all_broadpeaks_ch)



    // now putting the methylation bed file in the workflow with the up genes and down genes

    meth_enrichment_analysis_workflow(meth_bed_ch, wtvslowup_genebody_ch, wtvslowdown_nochange_ch, ref_genome_size_ch)



    // WILL HAVE TO FIGURE OUT A NAMING FOR THE PEAK FILES SO I CAN AUTOMATE THIS (IF HISTONE IN BIGWIG IS SAME AS HISTONE_LABEL IN PEAKS THEN RUN. SOMETHING LIKE THAT)
    // now lets plot histone marks signal (cpm normalized bigwig) over the up and down regulated peaks

    plot_signal_up_down_peaks_workflow(control_meta_cpm_bw, wt_meta_cpm_bw, up_peaks_ch, down_peaks_ch)

    // emitting the combined bigwig meta channel grouped by both histone and replicate
    // it has this format tuple(condition, histone, replicate, file_name, file)
    combined_bigwig_meta_2grouped_ch = plot_signal_up_down_peaks_workflow.out.experiment_group_meta_cpm_ch


    // now lets see how to get differential peaks to be plot over the up and down genes TSS plus 5kb
    // I will also plot the signal over genes tss plus 20kb in this workflow, because I already made the 20kb changed from 5kb
    // will put the combined bigwig channel here emitted from the other workflow
    plot_diff_peaks_over_diff_genes_workflow(up_peaks_ch, down_peaks_ch, unchanging_peaks_ch, master_peaks_ch, wtvslowup_genebody_ch, wtvslowdown_nochange_ch, gtf_ch, ref_genome_size_ch, knownGene_ch, proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch, combined_bigwig_meta_2grouped_ch)

    
    // now to add a workflow to plot the atac-seq signal over the cut&run peaks

    plot_atac_signal_over_diff_peaks_workflow(control_atac_bigwig_ch, treatment_atac_bigwig_ch, up_peaks_ch, down_peaks_ch, unchanging_peaks_ch )
    




}