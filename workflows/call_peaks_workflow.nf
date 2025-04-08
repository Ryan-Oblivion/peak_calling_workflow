




include {
    make_alignment_bw_process_control;
    make_alignment_bw_process_wt;
    plot_histone_at_genes_process;
    macs2_call_peaks_process_both
    // macs2_call_peaks_process_wt

}from '../modules/peak_analysis_modules.nf'




workflow mk_bw_call_peaks_workflow {


    take:
    control_bams_index_tuple_ch
    wt_bams_index_tuple_ch
    ref_genome_ch

    emit:
    control_meta_bw_ch
    wt_meta_bw_ch


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


    wt_bams_index_tuple_ch
        .map { unique_name, bam_index_tuple -> 
        
        bam_path = bam_index_tuple[0]
        bai_path = bam_index_tuple[1]

        basename = bam_path.baseName
        tokens = basename.tokenize("_")
        
        condition_label = tokens[0]

        histone_label = tokens[1]

        replicate_label = tokens[2]

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

    macs2_call_peaks_process_both(concat_wt_control_bam_meta_ch, ref_genome_ch) // might need ref_genome
    //macs2_call_peaks_process_wt()


}

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