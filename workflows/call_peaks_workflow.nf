



include {
    make_alignment_bw_process_control;
    make_alignment_bw_process_wt;
    plot_histone_at_genes_process;
    merge_peaks_bedtools_process;
    macs2_call_peaks_process_both;
    plot_histones_at_peaks_process;
    find_idr_in_replicates_process;
    multiqc_process;
    mk_bedgraph_process;
    seacr_peakcalls_process;
    sicer2_peakcall_process;
    mk_bed_for_sicer2_process;
    get_pval_bedgraph;
    kenttools_get_bigwig;
    find_diff_peaks_R_process;
    plot_at_up_down_peaks_process;
    signal_over_gene_tss_process;
    bedtools_stranded_create_process;
    merge_concat_peaks_process;
    diff_peaks_intersect_diff_genes_process;
    atac_signal_over_peaks_process;
    get_CpG_islands_in_peaks_process;
    plot_over_diff_cpg_regions_process;
    atac_enrich_counts_process;
    r_atac_enrich_plot_process
    // macs2_call_peaks_process_wt
    

}from '../modules/peak_analysis_modules.nf'




workflow mk_bw_call_peaks_workflow {


    take:
    control_bams_index_tuple_ch
    wt_bams_index_tuple_ch
    ref_genome_ch
    ref_genome_size_ch
    dups_log_ch

    


    main:

    // first I need to separate the control and wt by their histone marks

    control_bams_index_tuple_ch
        .map { unique_name, bam_index_tuple -> 
        
        bam_path = bam_index_tuple[0]
        bai_path = bam_index_tuple[1]

        basename = bam_path.baseName
        tokens = basename.tokenize("_")
        
        condition_label = tokens[0]

        histone_label = tokens[1]

        replicate_label = tokens[2]

        bio_label = tokens[3]

        bam_file_name = bam_path.name

        meta_name = "${tokens[0]}_${tokens[1]}_${tokens[2]}_${tokens[3]}"



        tuple(condition_label, histone_label, replicate_label, meta_name, bam_file_name, bam_path, bai_path)//.join(','))


        }
        .groupTuple(by:1, sort: true) // how the grouping would look  control h3k27me3 [[H1low, H1low, H1low], H3k27me3, [r2, r3, r1], [H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam, H1low_H3k27me3_r3_S27_001.trim.st.all.blft.qft.rmdup.sorted.bam, H1low_H3k27me3_r1_S25_001.trim.st.all.blft.qft.rmdup.sorted.bam],
        // .multiMap{row ->

        // histone_key = row[1] // hoping the second element in the grouping list is the histone key. look above at a snippet of how grouping looked

        // result = [:]
        // if (histone_key.contains('H3k27me3')) { result.h3k27me3 = row} 

        // if (histone_key.contains('H3k9me3')) { result.h3k9me3 = row}

        // return result

        // // h3k27me3: histone_key.contains('H3k27me3') ? row : null
        // // h3k9me3: histone_key.contains('H3k9me3') ? row : null

        // }
        .set{control_bam_meta_ch} // how one of the channels look [H1low,H3k27me3,r2,H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam,/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/bin/H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam,/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/bin/H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam.bai]
        //.view()
    //control_bam_meta_ch.view { v -> "control $v"} // not using multimap. here is how it looks.  [[H1low, H1low, H1low], H3k27me3, [r2, r3, r1], [H1low_H3k27me3_r2_S26_001.trim.st.all.blft.qft.rmdup.sorted.bam, H1low_H3k27me3_r3_S27_001.trim.st.all.blft.qft.rmdup.sorted.bam, H1low_H3k27me3_r1_S25_001.trim.st.all.blft.qft.rmdup.sorted.bam],

    //control_bam_meta_ch.view()

    wt_bams_index_tuple_ch
        .map { unique_name, bam_index_tuple -> 
        
        bam_path = bam_index_tuple[0]
        bai_path = bam_index_tuple[1]

        basename = bam_path.baseName
        tokens = basename.tokenize("_")
        
        condition_label = tokens[0]

        histone_label = tokens[1]

        replicate_label = tokens[2]

        bio_label = tokens[3]

        bam_file_name = bam_path.name

        meta_name = "${tokens[0]}_${tokens[1]}_${tokens[2]}_${tokens[3]}"

        

        tuple(condition_label, histone_label, replicate_label, meta_name, bam_file_name, bam_path, bai_path)//.join(', '))


        }
        .groupTuple(by:1, sort: true) // how the grouping would look [Scrm, Scrm, Scrm], H3k27me3, [r3, r2, r1], [Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.sorted.bam, Scrm_H3k27me3_r2_S2_001.trim.st.all.blft.qft.rmdup.sorted.bam, Scrm_H3k27me3_r1_S1_001.trim.st.all.blft.qft.rmdup.sorted.bam],
        //.view()
        // .multiMap{row ->

        // histone_key = row[1] // hoping the second element in the grouping list is the histone key. look above at a snippet of how grouping looked

        // result = [:]
        // if (histone_key.contains('H3k27me3')) { result.h3k27me3 = row} 

        // if (histone_key.contains('H3k9me3')) { result.h3k9me3 = row}

        // return result

        // //h3k27me3: histone_key.contains('H3k27me3') ? row : null
        // //h3k9me3: histone_key.contains('H3k9me3') ? row : null

        // }
        .set{wt_bam_meta_ch} // how one of the wt channels look with out grouping, but i ended up grouping so look above [Scrm,H3k27me3,r3,Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.sorted.bam,/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/bin/Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.sorted.bam,/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/bin/Scrm_H3k27me3_r3_S3_001.trim.st.all.blft.qft.rmdup.sorted.bam.bai]
    // wt_bam_meta_ch.h3k27me3.view { v -> "h3k27me3 $v"}
    // wt_bam_meta_ch.h3k9me3.view { v -> "h3k9me3 $v"}


    // first I need to make a process to make bigwig files out of each bam in the histone marks groups.
    // these bigwig files will be made from the bam files 
    // all_control_bams = control_bam_meta_ch.h3k27me3.concat(control_bam_meta_ch.h3k9me3)
    // all_control_bams.view()
    
    //control_bam_meta_ch.view()

    make_alignment_bw_process_control(control_bam_meta_ch)

    make_alignment_bw_process_wt(wt_bam_meta_ch)


    control_meta_bw_ch = make_alignment_bw_process_control.out.bigwig_meta_ch

    wt_meta_bw_ch = make_alignment_bw_process_wt.out.bigwig_meta_ch

    control_meta_cpm_bw_ch = make_alignment_bw_process_control.out.cpm_bigwig_meta_ch
    wt_meta_cpm_bw_ch = make_alignment_bw_process_wt.out.cpm_bigwig_meta_ch

    // I want to view how the meta bigwig output channel looks.
    //control_meta_bw_ch.view()

    // now that i have the meta channels where I grouped them by the histone marks, I can put it into a process to call peaks on all the files for that histone mark and another process will spawn calling peaks for the other histone marks in parallel

    // i think if i concat the control_bam_meta_ch and the wt_bam_meta_ch I can parallelize the process
    wt_bam_meta_ch
        .concat(control_bam_meta_ch)
        .transpose()
        //.view()
        .set{concat_wt_control_bam_meta_ch}
    //concat_wt_control_bam_meta_ch.view()
    macs2_call_peaks_process_both(concat_wt_control_bam_meta_ch, ref_genome_ch) // might need ref_genome
    //macs2_call_peaks_process_wt()

    // i want to get the ppois files from the macs2 process
    ppois_files_ch = macs2_call_peaks_process_both.out.ppois_macs2_file

    get_pval_bedgraph(ppois_files_ch, ref_genome_size_ch)

    pval_bedgraph_ch = get_pval_bedgraph.out.pvalue_bedgraph_file

    chrom_size_ch = get_pval_bedgraph.out.chrom_size_file

    kenttools_get_bigwig(pval_bedgraph_ch, ref_genome_size_ch)
    //chrom_size_ch


    // now getting the channel for the broadpeaks and will have to emit it from the workflow and put into a new workflow to plot the bigwig signal onto the called broad peaks

    broadpeaks_ch = macs2_call_peaks_process_both.out.broadpeaks
    broadpeaks_ch_backup = macs2_call_peaks_process_both.out.broadpeaks
    //broadpeaks_ch.view() // now got the peaks so time to emit this channel.

    // okay, another thing to do is to use bedtools to merge peaks by 1kb, 2kb, and 5kb to see how they look

    merge_peaks_bedtools_process(broadpeaks_ch)




    broadpeaks_ch_backup
        .map{ peakpath -> 
        
        basename = peakpath.baseName
        file_name = peakpath.name

        tokens = basename.tokenize("_")

        condition = tokens[0]
        histone = tokens[1]
        replicate = tokens[2]
        bio_rep = tokens[3]

        grouping_key = "${condition}_${histone}"

        tuple(grouping_key, condition, histone, replicate, bio_rep, file_name, basename, peakpath)
        
        }
        .groupTuple(by:0, sort:true)
        //.view()
        .set{broadpeak_gtuple_meta_ch}

    // split the broadpeaks so i have them grouped by their condition label


    // now i want to get the idr peaks per each replicate combination
    //broadpeak_gtuple_meta_ch.view()
    find_idr_in_replicates_process(broadpeak_gtuple_meta_ch)

    // these were the peaks without 10kb merged
    concat_idr_peaks = find_idr_in_replicates_process.out.final_concat_peaks

    
    // now with the concat peaks, I need to put them in R and get the list of up peaks and down peaks

    concat_idr_peaks
        .map {file -> 
        
        basename = file.baseName

        file_name = file.name

        // concat_IDR_H1low_H3k27me3_r1_vs_r2_vs_r3
        tokens = basename.tokenize("_")

        condition = tokens[2]
        histone = tokens[3]

        tuple(histone, condition, file_name, file)

        
        
        }
        .groupTuple(by:0, sort: true)
        //.view()
        // H1low is first and scrambled is next but if the file is named something else you should not hard code which is which
        // example: [H3k27me3, [conditions], [concat_IDR_H1low_H3k27me3_r1_vs_r2_vs_r3_0.4_pairs.broadPeak, concat_IDR_Scrm_H3k27me3_r1_vs_r2_vs_r3_0.4_pairs.broadPeak], [/lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/work/32/d1583cfe217bf08b629e50f5d2c843/concat_IDR_H1low_H3k27me3_r1_vs_r2_vs_r3_0.4_pairs.broadPeak, /lustre/fs4/risc_lab/scratch/rjohnson/pipelines/peak_calling_analysis_pipeline/work/72/7a61681f524cbe7de6d547ad41a86d/concat_IDR_Scrm_H3k27me3_r1_vs_r2_vs_r3_0.4_pairs.broadPeak]]
        .set{group_concat_idr_peaks_ch}

    // now I should merge these concat peaks but keep them separate and output them in the same meta channel as above

    merge_concat_peaks_process(group_concat_idr_peaks_ch)

    first_10kb_peakfile = merge_concat_peaks_process.out.first_10kb_merged_peak
    second_10kb_peakfile = merge_concat_peaks_process.out.second_10kb_merged_peak
    //first_10kb_peakfile.view()

    both_10kb_peakfiles = first_10kb_peakfile.concat(second_10kb_peakfile)

    both_10kb_peakfiles
        .map { file -> 
        
        file_basename = file.baseName
        file_name = file.name

        tokens = file_basename.tokenize("_")

        condition_label = tokens[2]
        exper_type = tokens[3]

        tuple(condition_label, exper_type, file_name, file)

        }
        .groupTuple(by:1, sort:true)
        .set{group_10kb_concat_idr_peaks_ch}


    //group_10kb_concat_idr_peaks_ch = merge_concat_peaks_process.out.merged_10kb_concat_peaks
    //group_10kb_concat_idr_peaks_ch.view()


    // now I also need the bam files
    // load all bams in but use the histone mark to get the correct bams
    // the bam meta ch looks like this tuple(condition_label, histone_label, replicate_label, meta_name, bam_file_name, bam_path, bai_path)
    
     control_bams_index_tuple_ch
        .concat(wt_bams_index_tuple_ch)
        .map { key, tuple ->
    
        
            bam = tuple[0]
            bai = tuple[1]
            basename = bam.baseName
            file_name = bam.name

            tokens = basename.tokenize("_")

            condition = tokens[0]
            histone = tokens[1]

            bam
            //tuple(histone, condition, file_name, bam)
        
        }
        .collect()
        // just get all bams
        // . map { bam ->
        
        // basename = bam.baseName
        // file_name = bam.name

        // tokens = basename.tokenize("_")

        // condition = tokens[0]
        // histone = tokens[1]

        // tuple(histone, condition, file_name, bam)
        
        // }
        // .groupTuple(by:0, sort:true)
        //.view()
        //.set{meta_bam_histone_group_tuple_ch}
        .set{all_bams_paths}
    
    //all_bams_paths.view()
    //meta_bam_histone_group_tuple_ch.view()
    
    //group_concat_idr_peaks_ch.view()
    
    // filtering channels 
    // h3k27me3_idr_peaks_ch = group_concat_idr_peaks_ch.filter { it[0] == 'H3k27me3' }
    // h3k27me3_bams_ch = meta_bam_histone_group_tuple_ch.filter { it[0] == 'H3k27me3' }

    // h3k27me3_idr_peaks_ch.view()
    // h3k27me3_bams_ch.view()

    // using the 10kb merged idr peaks
    find_diff_peaks_R_process(group_10kb_concat_idr_peaks_ch, all_bams_paths)

    diff_peaks_tuple = find_diff_peaks_R_process.out.diff_peaks_ch
    //diff_peaks_tuple.view()

    // try keeping the peaks in separate channels
    // then when I put this into a process or workflow, just check to see if the histones match
    master_peaks_list_ch = find_diff_peaks_R_process.out.master_peak_emit.collect()
    up_peaks_list_ch = find_diff_peaks_R_process.out.up_peaks_emit.collect()
    down_peaks_list_ch = find_diff_peaks_R_process.out.down_peaks_emit.collect()
    unchanging_peaks_list_ch = find_diff_peaks_R_process.out.unchanging_peaks_emit.collect()

    //up_peaks_list_ch
        //.collect()
        //.view()



    if (params.make_html_report = true) {
        // running multiqc on the duplicate log files

        // group by histone 
        dups_log_ch
            .flatten()
            .map { file ->
            
            file_basename = file.baseName
            tokens = file_basename.tokenize("_")
            condition = tokens[0]
            histone = tokens[1]
            
            //grouping_key  = "${condition}_${histone}"
            tuple( condition, histone, file)
            }
            .groupTuple(by:1)
            .set{dups_log_meta_ch}
        multiqc_process(dups_log_meta_ch)
    }


    // I should try seacr now //////////////
    //control_bam_meta_ch.view()
    //wt_bam_meta_ch.view()
    
    // possibly concate the two above so i can run process for getting bedgraph files in parallel
    control_bams_index_tuple_ch
        .concat(wt_bams_index_tuple_ch)
        //.view()
        .set{combined_bam_index_tuple_ch}

    //combined_bam_index_tuple_ch.view()
    // it takes the reference genome also
    mk_bedgraph_process(combined_bam_index_tuple_ch, ref_genome_size_ch)

    bedgraphs_for_seacr_ch = mk_bedgraph_process.out.bedgraph_for_seacr
    
    // then the seacr process 
    //seacr_peakcalls_process(bedgraphs_for_seacr_ch)

    // now make a idr process for seacr peaks
    // well seacr does not give back peaks that when merged with idr will total more that 20
    // and idr needs to have a post merge of 20 peaks.
    
    /* seacr_peaks = seacr_peakcalls_process.out.seacr_peaks

    seacr_peaks
        .map{ peakpath -> 
        
        basename = peakpath.baseName
        file_name = peakpath.name

        

        tokens = basename.tokenize("_")

        condition = tokens[0]
        histone = tokens[1]
        replicate = tokens[2]
        bio_rep = tokens[3]

        grouping_key = "${condition}_${histone}"

        tuple(grouping_key, condition, histone, replicate, bio_rep, file_name, basename, peakpath)
        
        }
        .groupTuple(by:0, sort:true)
        //.view()
        .set{seacrpeaks_gtuple_meta_ch}



    seacr_idr_process(seacrpeaks_gtuple_meta_ch)
    */


    // now making a process to get sicer2 peaks
    // I can use the bam files because I added bedtools in the sicer2 conda environment I made

    //combined_bam_index_tuple_ch.view()
    //sicer2_peakcall_process(combined_bam_index_tuple_ch)

    // ill just make a process to create bed files from the bams and then input it into sicer2 env that is by itself
    // mk_bed_for_sicer2_process(combined_bam_index_tuple_ch)

    // bed_for_sicer2_ch = mk_bed_for_sicer2_process.out.bed_for_sicer2

    //bed_for_sicer2_ch
        //.filter ( ~/.*H1low.*/)
        //.view()
        //.set{hlow_bed_for_sicer2_ch}
    
    //sicer2_peakcall_process(hlow_bed_for_sicer2_ch)

    emit:
    control_meta_bw_ch
    wt_meta_bw_ch
    broadpeaks_ch
    control_meta_cpm_bw_ch
    wt_meta_cpm_bw_ch
    group_10kb_concat_idr_peaks_ch
    master_peaks_list_ch
    up_peaks_list_ch
    down_peaks_list_ch
    unchanging_peaks_list_ch
    //diff_peaks_tuple
    //concat_idr_peaks

}


/*
workflow get_diff_peaks_workflow {


    take:

    all_broadpeaks_ch



    main:

    //all_broadpeaks_ch.view{file -> "this is the broadpeak: $file"}




}*/

workflow plot_histone_data_workflow {


    take:
    control_bw_meta_ch
    wt_bw_meta_ch
    wtvslowup_genebody_ch
    wtvslowdown_nochange_ch




    main:

    control_bw_meta_ch
        .map {condition_label, histone_label, replicate_label, bigwig_path ->

        bigwig_name = bigwig_path.name  // this is vectorized. it got the file name for all three in the list removing the paths
        //bigwig_name2 = bigwig_path[[1]].name 
        //bigwig_name3 = bigwig_path[[2]].name 

        file_basename = bigwig_path.baseName

        // tokens = file_basename.tokenize("_")

        // condition = tokens[0]
        // histone = tokens[1]
        // replicate = tokens[2]
        // hopefully i can end with this keeping k9 and k27 separated still

        //tuple(condition_label, histone_label, replicate_label, bigwig_name, bigwig_path)
        //tuple(bigwig_name,bigwig_path)
        bigwig_path

        }
        //.sort{ a, b -> a[2] <=> b[2]}
        //.transpose()
        //.groupTuple(by:1, sort:true)
        //.view()
        .set{control_bw_meta2_ch}

    wt_bw_meta_ch
        .map {condition_label, histone_label, replicate_label, bigwig_path ->

        bigwig_name = bigwig_path.name  // this is vectorized. it got the file name for all three in the list removing the paths
        //bigwig_name2 = bigwig_path[[1]].name 
        //bigwig_name3 = bigwig_path[[2]].name 

        file_basename = bigwig_path.baseName

        // tokens = file_basename.tokenize("_")

        // condition = tokens[0]
        // histone = tokens[1]
        // replicate = tokens[2]
        // hopefully i can end with this keeping k9 and k27 separated still

        //tuple(condition_label, histone_label, replicate_label, bigwig_name, bigwig_path)
        //tuple(bigwig_name,bigwig_path)
        bigwig_path

        }
        //.sort{ a, b -> a[2] <=> b[2]}
        //.transpose()
        //.groupTuple(by:1, sort:true)
        //.view()
        .set{wt_bw_meta2_ch}
        
    wt_bw_meta2_ch
        .concat(control_bw_meta2_ch)
        .flatten()
        .map{ file ->

        file_basename = file.baseName
        file_name = file.name
        tokens = file_basename.tokenize("_")

        condition = tokens[0]
        histone = tokens[1]
        replicate = tokens[2]

        tuple(condition, histone, replicate, file_name, file)


        }
        .groupTuple(by:1)
        //.view()
        .set {experiment_group_meta_ch}
    // making a process that takes both control and wt and makes a plot. I need to give it a gene list also

    plot_histone_at_genes_process(experiment_group_meta_ch, wtvslowup_genebody_ch, wtvslowdown_nochange_ch)





}


workflow plot_histone_calledpeak_workflow {



    take:
    control_bw_meta_ch2
    wt_bw_meta_ch2
    all_broadpeaks_ch2






    main:

    // i need to make a meta channel for the broadPeaks where I have in a list the broadPeaks for each condition and histone
    
    control_bw_meta_ch2
        .concat(wt_bw_meta_ch2)
        .map { condition, histone, replicate, bigwig_paths ->

        
        bigwig_paths
        
        }
        .flatten()
        .concat(all_broadpeaks_ch2)
        .map { paths ->

        file_basenames = paths.baseName

        file_name = paths.name
        
        tokens = file_basenames.tokenize("_")

        condition = tokens[0]

        histone = tokens[1]

        replicate = tokens[2]

        // i need to get the bio rep

        bio_rep = tokens[3]

        grouping_key = "${condition}_${histone}_${replicate}_${bio_rep}"

        tuple(grouping_key, condition, histone, replicate, file_name, paths)

        }
        .groupTuple(by:0)
        //.view()
        .set{meta_bw_peak_ch}
    //meta_bw_peak_ch.view()
    plot_histones_at_peaks_process(meta_bw_peak_ch)
    
    
    // all_broadpeaks_ch2
    //     .map {peak_paths ->

    //     peak_basename = peak_paths.baseName

    //     peak_filename = peak_paths.name

    //     tokens = peak_basename.tokenize("_")

    //     condition_label = tokens[0]
    //     histone_label = tokens[1]
    //     replicate_label = tokens[2]

    //     grouping_key = "${condition_label}_${histone_label}"

    //     tuple( condition_label, histone_label, replicate_label, peak_filename, peak_paths)


    //     }
    //     .groupTuple(by:1)
    //     // .multiMap {tuple ->

    //     // histone = tuple[1]

    //     // k27me3: histone=="H3k27me3"
    //     // k9me3: histone=="H3k9me3"

    //     // }
    //     .filter { tuple ->
        
    //     //tuple(condition, histone, replicate, peak_filename, peak_basename, peak_paths)
    //     tuple[1]=="H3k27me3"
        

    //     }
    //     //.transpose()
    //     //.view()
    //     .set{broadpeak_k27_ch}
    
    // // now getting the k9 channel

    // all_broadpeaks_ch2
    //     .map {peak_paths ->

    //     peak_basename = peak_paths.baseName

    //     peak_filename = peak_paths.name

    //     tokens = peak_basename.tokenize("_")

    //     condition_label = tokens[0]
    //     histone_label = tokens[1]
    //     replicate_label = tokens[2]

    //     grouping_key = "${condition_label}_${histone_label}"

    //     tuple( condition_label, histone_label, replicate_label, peak_filename, peak_paths)


    //     }
    //     .groupTuple(by:1)
    //     // .multiMap {tuple ->

    //     // histone = tuple[1]

    //     // k27me3: histone=="H3k27me3"
    //     // k9me3: histone=="H3k9me3"

    //     // }
    //     .filter { tuple ->
        
    //     //tuple(condition, histone, replicate, peak_filename, peak_basename, peak_paths)
    //     tuple[1]=="H3k9me3"
        

    //     }
    //     //.transpose()
    //     //.view()
    //     .set{broadpeak_k9_ch} // how it looks [[Scrm, Scrm, Scrm, H1low, H1low, H1low], H3k9me3, [r1, r3, r2, r1, r3, r2], 


    //     //.set{broadpeak_meta_ch}

    

    // control_bw_meta_ch2
    //     .map {condition_label, histone_label, replicate_label, bigwig_path ->

    //     bigwig_name = bigwig_path.name  // this is vectorized. it got the file name for all three in the list removing the paths
        

    //     file_basename = bigwig_path.baseName

        
    //     bigwig_path

    //     }
        
    //     .set{control_bw_meta3_ch}

    // wt_bw_meta_ch2
    //     .map {condition_label, histone_label, replicate_label, bigwig_path ->

    //     bigwig_name = bigwig_path.name  // this is vectorized. it got the file name for all three in the list removing the paths
         

    //     file_basename = bigwig_path.baseName

    //     bigwig_path

    //     }
        
    //     .set{wt_bw_meta3_ch}
        
    // wt_bw_meta3_ch
    //     .concat(control_bw_meta3_ch)
    //     //.concat(broadpeak_meta_ch)
    //     //.view()
    //     .flatten()
    //     .map{ file ->

    //     file_basename = file.baseName
    //     file_name = file.name
    //     tokens = file_basename.tokenize("_")

    //     condition = tokens[0]
    //     histone = tokens[1]
    //     replicate = tokens[2]

    //     grouping_key = "${condition}_${histone}_${replicate}"

    //     tuple(grouping_key, condition, histone, replicate, file_name, file_basename, file)


    //     }
    //     .groupTuple(by:2)
    //     .map{grouping_key, condition, histone, replicate, file_name, file_basename, file ->

    //     tuple(condition, histone, replicate, file_name, file)

    //     }
    //     //.view()
    //     // .groupTuple(by:0)
    //     // //.transpose() 
    //     // .map { condition, histone, replicate, file_name, file_basename, file ->

    //     // tuple(file_name, file)

    //     // }
    //     // .view()
    //     .set {experiment_group_bigwigs_meta_ch}

    // // i need to separate k27 and k9 bigwigs

    // experiment_group_bigwigs_meta_ch
    //     .filter{ tuple ->

    //     histone = tuple[1]

    //     histone == "H3k27me3"
    //     }
    //     .transpose()
    //     //.view()
    //     .set{bigwig_k27_meta_ch}

    // // now for the bigwig k9

    // experiment_group_bigwigs_meta_ch
    //     .filter{ tuple ->

    //     histone = tuple[1]

    //     histone == "H3k9me3"
    //     }
    //     .transpose()
    //     //.view()
    //     .set{bigwig_k9_meta_ch}


    // // if i now concat the two channels i can hopefully be sure the first will be k27 and the second will be k9 if i put it that way

    // bigwig_k27_meta_ch
    //     .concat(bigwig_k9_meta_ch)
    //     .set{concat_bw_k27_k9_ch}

    // broadpeak_k27_ch.transpose()
    //     .concat(broadpeak_k9_ch.transpose())
    //     .set{concat_peak_k27_k9_ch}

    //concat_bw_k27_k9_ch.view()
    //concat_peak_k27_k9_ch.view()


    // see how this looks and how i can get the broad peaks to be in here

    //experiment_group_meta_ch.view()

    //plot_histones_at_peaks_process(concat_bw_k27_k9_ch, concat_peak_k27_k9_ch)


}

// this will be very similar to the plotting histone signal at genes workflow
workflow plot_signal_up_down_peaks_workflow {




    take:

    control_meta_cpm_bw
    wt_meta_cpm_bw
    up_peaks_ch
    down_peaks_ch
    bisulfate_bigwig_ch
    master_peaks_ch
    unchanging_peaks_list_true
    cpg_island_unmasked_ch



    main:

    // need to combine the two bigwig channels

    control_meta_cpm_bw
        .map {condition_label, histone_label, replicate_label, bigwig_path ->

        bigwig_name = bigwig_path.name  // this is vectorized. it got the file name for all three in the list removing the paths
        //bigwig_name2 = bigwig_path[[1]].name 
        //bigwig_name3 = bigwig_path[[2]].name 

        file_basename = bigwig_path.baseName

        // tokens = file_basename.tokenize("_")

        // condition = tokens[0]
        // histone = tokens[1]
        // replicate = tokens[2]
        // hopefully i can end with this keeping k9 and k27 separated still

        //tuple(condition_label, histone_label, replicate_label, bigwig_name, bigwig_path)
        //tuple(bigwig_name,bigwig_path)
        bigwig_path

        }
        //.sort{ a, b -> a[2] <=> b[2]}
        //.transpose()
        //.groupTuple(by:1, sort:true)
        //.view()
        .set{control_bw_cpm_meta2_ch}

    wt_meta_cpm_bw
        .map {condition_label, histone_label, replicate_label, bigwig_path ->

        bigwig_name = bigwig_path.name  // this is vectorized. it got the file name for all three in the list removing the paths
        //bigwig_name2 = bigwig_path[[1]].name 
        //bigwig_name3 = bigwig_path[[2]].name 

        file_basename = bigwig_path.baseName

        // tokens = file_basename.tokenize("_")

        // condition = tokens[0]
        // histone = tokens[1]
        // replicate = tokens[2]
        // hopefully i can end with this keeping k9 and k27 separated still

        //tuple(condition_label, histone_label, replicate_label, bigwig_name, bigwig_path)
        //tuple(bigwig_name,bigwig_path)
        bigwig_path

        }
        //.sort{ a, b -> a[2] <=> b[2]}
        //.transpose()
        //.groupTuple(by:1, sort:true)
        //.view()
        .set{wt_bw_cpm_meta2_ch}
        
    wt_bw_cpm_meta2_ch
        .concat(control_bw_cpm_meta2_ch)
        .flatten()
        .map{ file ->

        file_basename = file.baseName
        file_name = file.name
        tokens = file_basename.tokenize("_")

        condition = tokens[0]
        histone = tokens[1]
        replicate = tokens[2]

        tuple(condition, histone, replicate, file_name, file)


        }
        .groupTuple(by: [1,2])  // here i can try groupting tuple by 2 fields or columns. By the histone and then by the replicate !!!
        //.view()
        .set {experiment_group_meta_cpm_ch}

    // new group cpm ch for joining
    wt_bw_cpm_meta2_ch
        .concat(control_bw_cpm_meta2_ch)
        .flatten()
        .map{ file ->

        file_basename = file.baseName
        file_name = file.name
        tokens = file_basename.tokenize("_")

        condition = tokens[0]
        histone = tokens[1]
        replicate = tokens[2]

        tuple(condition, histone, replicate, file_name, file)


        }
        .groupTuple(by: [1,2])  // here i can try groupting tuple by 2 fields or columns. By the histone and then by the replicate !!!
        //.view()
        .set {experiment_group_meta_to_join_ch}

    // try to make a meta map of up down masterpeaks
    up_peaks_ch
        .concat(down_peaks_ch)
        .concat(unchanging_peaks_list_true)
        .concat(master_peaks_ch)
        .flatten()
        .map {peaks -> 
        
        peak_basename = peaks.baseName

        tokens = peak_basename.tokenize("_")

        diff_peak = tokens[0] // this will be up, down or master
        exper_type = tokens[1] // this will be k27me3 or k9me3 or anytype of experiment

        tuple(diff_peak, exper_type, peaks)


        }
        .groupTuple(by:1)
        //.view()
        .set{peaks_group_meta_ch}

    experiment_group_meta_to_join_ch
        //.map { histone, condition, replicate, file_name, file ->
        //tuple(histone, tuple(histone, condition, replicate, file_name, file))

        //}
        .combine(peaks_group_meta_ch, by:1)

        //.map{histone, meta_tuple, peaks ->

        //def (expr_type, conditions, replicate, file_name, files) = meta_tuple

        //tuple(expr_type, conditions, replicate, file_name, files, peaks)

        //}

        .view()
        .set{exper_rep_bigwig_peak_group}
    // now using .join to group the two channels

    // now to plot the signal at the peaks
    
    // this below is the version that used the geo data with the h1low to make masterpeaks
    //plot_at_up_down_peaks_process(experiment_group_meta_cpm_ch, up_peaks_ch, down_peaks_ch, bisulfate_bigwig_ch, master_peaks_ch, cpg_island_unmasked_ch)

    ////////// uncomment when using only master peaks generated from data only in this pipeline //////////
    plot_at_up_down_peaks_process(exper_rep_bigwig_peak_group, bisulfate_bigwig_ch, cpg_island_unmasked_ch)
    // here i can try groupting tuple by 2 fields or columns. By the histone and then by the replicate !!!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////



    emit:

    experiment_group_meta_cpm_ch

    // this is the new combined bigwig meta ch that is grouped by histone and replicate but also has the peaks files there also
    exper_rep_bigwig_peak_group

}

workflow plot_diff_peaks_over_diff_genes_workflow {



    take:
    up_peaks_ch
    down_peaks_ch
    unchanging_peaks_ch
    master_peaks_ch
    up_genes_ch
    down_genes_ch
    nochanging_genes_ch
    gtf_ch
    ref_genome_size_ch
    knownGene_ch
    proseq_up_gene_ch
    proseq_down_gene_ch
    proseq_unchanging_gene_ch
    combined_bigwig_meta_2grouped_ch // this will have the peak type and the peaks paths if using the version that only uses data generated in the pipeline
    


    main:

    // first lets get the the strandedness of the up and down genes then use bedtools slop to get 5kb based on strand

    bedtools_stranded_create_process(up_genes_ch, down_genes_ch, nochanging_genes_ch, gtf_ch, ref_genome_size_ch) 

    // changed from 5kb to 20kb
    up_genes_with_20kb = bedtools_stranded_create_process.out.up_genes_20kb_stranded
    down_genes_with_20kb = bedtools_stranded_create_process.out.down_genes_20kb_stranded
    nochange_genes_with_20kb = bedtools_stranded_create_process.out.nochange_genes_20kb_stranded


    // I might just concat the up and down peak files to get a large diff peak file then plot over up and down genes
    // the combined_bigwig_meta_2grouped_ch is grouped by histone and replicate and has this format tuple(condition, histone, replicate, file_name, file)
    // this will have independent channels where in each histone, you have a single replicate with its treatment bigwig and its control bigwig. this compares the control and treatment in the same replicate and same histone
    
    //signal_over_gene_tss_process(up_peaks_ch, down_peaks_ch, up_genes_with_20kb, down_genes_with_20kb, nochange_genes_with_20kb, up_genes_ch, down_genes_ch, combined_bigwig_meta_2grouped_ch) //version for outside data
    signal_over_gene_tss_process( up_genes_with_20kb, down_genes_with_20kb, nochange_genes_with_20kb, up_genes_ch, down_genes_ch, combined_bigwig_meta_2grouped_ch)

    // now make a process to find which peaks intersect the up and down genes
    
    //not using for now
    //diff_peaks_intersect_diff_genes_process(up_peaks_ch, down_peaks_ch, unchanging_peaks_ch, master_peaks_ch, up_genes_with_20kb, down_genes_with_20kb, nochange_genes_with_20kb, knownGene_ch, proseq_up_gene_ch, proseq_down_gene_ch, proseq_unchanging_gene_ch)



}


workflow plot_atac_signal_over_diff_peaks_workflow {


    take:

    control_atac_bigwig

    treatment_atac_bigwig

    up_peaks

    down_peaks

    unchanging_peaks

    // new additions
    down_atac_peaks_ch
    up_atac_peaks_ch
    combined_bigwig_meta_2grouped_ch
    cpg_island_unmasked_ch

    main:

    // now make a process that will plot this information

    // with the histone meta bigwig channels, just merge the three replicates into one for each condition in the process using cat

    // version using peaks from outside pipeline
    // atac_signal_over_peaks_process(control_atac_bigwig, treatment_atac_bigwig, up_peaks, down_peaks, unchanging_peaks,  down_atac_peaks_ch, up_atac_peaks_ch, combined_bigwig_meta_2grouped_ch, cpg_island_unmasked_ch)

    // version using peaks from the pipeline
    atac_signal_over_peaks_process(control_atac_bigwig, treatment_atac_bigwig, down_atac_peaks_ch, up_atac_peaks_ch, combined_bigwig_meta_2grouped_ch, cpg_island_unmasked_ch)



    //emit:
}

workflow find_then_plot_cpgIslands_in_peaks_workflow {


    take:

    up_peaks_ch
    down_peaks_ch
    unchanging_peaks_ch
    cpg_island_unmasked_ch
    combined_bigwig_meta_2grouped_ch




    main:

    // now first lets find which CpG islands overlap with the different peaks


    get_CpG_islands_in_peaks_process(up_peaks_ch, down_peaks_ch, unchanging_peaks_ch, cpg_island_unmasked_ch)

    cpg_up_ch = get_CpG_islands_in_peaks_process.out.cpg_up_regions
    cpg_down_ch = get_CpG_islands_in_peaks_process.out.cpg_down_regions
    cpg_unchanging_ch = get_CpG_islands_in_peaks_process.out.cpg_unchanging_regions

    // now to plot the signal over these up, down, unchanging cpg regions
    plot_over_diff_cpg_regions_process(cpg_up_ch, cpg_down_ch, cpg_unchanging_ch, combined_bigwig_meta_2grouped_ch)


}


workflow get_roadmap_histone_enrichment_workflow {


    take:
    roadmap_broad_histones
    roadmap_narrow_histones
    control_atac_bigwig_ch
    treatment_atac_bigwig_ch



    main:

    roadmap_broad_histones
        .map {file ->
        
        basename = file.baseName
        name = file.name

        tuple(name, basename, file)


        }
        //.view()
        .set{roadmap_broad_histone_meta_ch}

    roadmap_narrow_histones
        .map {file ->
        
        basename = file.baseName
        name = file.name

        tuple(name, basename, file)


        }
        //.view()
        .set{roadmap_narrowhistone_meta_ch}
    // now i should combine the two atac bigwig files then flatten them so they are run in parallel

    control_atac_bigwig_ch
        .concat(treatment_atac_bigwig_ch)
        //.flatten()
        // .map{ file ->
        
        // basename = file.baseName

        // tokens = basename.tokenize("_")

        // condition_label = tokens[0]
        // experiment_label = tokens[1]
        // replicate_label = tokens[2]
        // tuple(condition_label, experiment_label, replicate_label)
        // }
        .view()
        .set{atac_bigwig_cat}
    
    // now make a deeptools process that uses multibigwigsummary to get the counts

    atac_enrich_counts_process(roadmap_broad_histone_meta_ch, roadmap_narrowhistone_meta_ch, atac_bigwig_cat)

    enrich_counts_tab_ch = atac_enrich_counts_process.out.raw_enrichment_counts.collect()

    enrich_counts_tab_ch
        .flatten()
        .map{ file ->
        // i need to get the tokens and find which used broad and narrow then group them by that, so i have scrm with h1low broad, and scrm with h1low narrow
        basename = file.baseName

        filename = file.name

        tokens = basename.tokenize("_")

        condition = tokens[0]
        peak_type = tokens[1]

        tuple(peak_type, basename, filename, file)

        }
        .groupTuple(by:0)
        .view(tuple -> "This is the meta channel for enrichment counts tab file grouped by peak type: $tuple")
        .set{enrich_counts_tab_meta_ch}
        
    // now I want to make a process that will take the output tab files from above and make the enrichemnt in R

    r_atac_enrich_plot_process(enrich_counts_tab_meta_ch)


    //emit:
}