nextflow.enable.dsl=2


params.control_bams = 'bam_files/H1low_*{bam,bam.bai}'
control_bams_index_tuple_ch = Channel.fromFilePairs(params.control_bams)

//control_bams_index_tuple_ch.view()


params.wt_bams = 'bam_files/Scrm_*{bam,bam.bai}'
wt_bams_index_tuple_ch = Channel.fromFilePairs(params.wt_bams)


params.ref_genome = file('/lustre/fs4/risc_lab/store/risc_data/downloaded/hg38/genome/Sequence/WholeGenomeFasta/genome.fa')
ref_genome_ch = Channel.value(params.ref_genome)

// I need a list of gene regions. got this from irene but also included it in my motif analysis nextflow script
params.wtvs_lowup = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/RNA-seq/rep2/hg38-ERCC-UMI-alignment/DESeq2_results/WTvslowup-genebody.bed')
wtvslowup_genebody_ch = Channel.value(params.wtvs_lowup)

// including the nochange genes just to have something to compare to
params.wtvs_lowdown_nochange = file('/lustre/fs4/risc_lab/scratch/iduba/linker-histone/RNA-seq/rep2/hg38-ERCC-UMI-alignment/DESeq2_results/WTvslowdown-basemeanmatchnochange-genebody.bed')
wtvslowdown_nochange_ch = Channel.value(params.wtvs_lowdown_nochange)

// include {
//     make_alignment_bw_process_control

// }from './modules/peak_analysis_modules.nf'



include {
    mk_bw_call_peaks_workflow;
    plot_histone_data_workflow;
    plot_histone_calledpeak_workflow;
    get_diff_peaks_workflow

}from './workflows/call_peaks_workflow.nf'



include {

}from './workflows/find_diff_peaks_workflow.nf'




workflow {

    // now i will put the control bams and wt bams into the peak calling workflow
    // I will also put the reference genome in as the third entry

    mk_bw_call_peaks_workflow(control_bams_index_tuple_ch, wt_bams_index_tuple_ch, ref_genome_ch )

    // take the emitted channels from the call peaks workflow
    control_bw_meta_ch = mk_bw_call_peaks_workflow.out.control_meta_bw_ch

    wt_bw_meta_ch = mk_bw_call_peaks_workflow.out.wt_meta_bw_ch

    // the mk_bw_call_peaks_workflow workflow also emits the broad peaks that were called so I will grab those
    all_broadpeaks_ch = mk_bw_call_peaks_workflow.out.broadpeaks_ch

    // adding a workflow here to get the differential peaks

    get_diff_peaks_workflow(all_broadpeaks_ch)


    // here I want to make a workflow that plots the chromatin features at a list of annotated sites (genes)
    // Also looking at the peaks around tss

    plot_histone_data_workflow(control_bw_meta_ch, wt_bw_meta_ch, wtvslowup_genebody_ch, wtvslowdown_nochange_ch)

    // I will just make another workflow for plotting each histone-rep bigwig with its corresponding histone-rep broadPeak

    plot_histone_calledpeak_workflow(control_bw_meta_ch, wt_bw_meta_ch, all_broadpeaks_ch)









}