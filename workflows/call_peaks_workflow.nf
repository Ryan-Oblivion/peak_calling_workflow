




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
    mk_bed_for_sicer2_process
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
    seacr_peakcalls_process(bedgraphs_for_seacr_ch)

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