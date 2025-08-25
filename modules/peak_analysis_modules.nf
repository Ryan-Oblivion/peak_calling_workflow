

process make_alignment_bw_process_control {

    debug true

    label 'normal_small_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'

    publishDir "./bigwigs/${histone_label}_bigwigs/", mode: 'copy', pattern: "*${histone_label}*.bigwig"


    input:
    
    // not using multimap so i will have two instances of the process.

    tuple val(condition_label), val(histone_label), val(replicate_label), val(meta_name), val(bam_file_name), path(bam_path), path(bai_path)
    
    // if i were using multimap then i would do this below
    //tuple val(condition_label_k27), val(histone_label_k27), val(replicate_label_k27), val(bam_file_name_k27), path(bam_path_k27), path(bai_path_k27)

    //tuple val(condition_label_k9), val(histone_label_k9), val(replicate_label_k9), val(bam_file_name_k9), path(bam_path_k9), path(bai_path_k9)


    output:

    tuple val("${condition_label}"), val("${histone_label}"), val("${replicate_label}"), path("*_raw*.bigwig"), emit: bigwig_meta_ch

    // the cpm normalized bigwigs 
    tuple val("${condition_label}"), val("${histone_label}"), val("${replicate_label}"), path("*_cpm*.bigwig"), emit: cpm_bigwig_meta_ch

    


    script:

    // there are three  technical replicates for each label above, except histone_label 

    // rep1_file = "${bam_file_name[0]}"
    // rep2_file = "${bam_file_name[1]}"
    // rep3_file = "${bam_file_name[2]}"

    // rep1_meta_name = "${meta_name[0]}"
    // rep2_meta_name = "${meta_name[1]}"
    // rep3_meta_name = "${meta_name[2]}"


    // rep1_k27_file = "${bam_file_name_k27[2]}"
    // rep2_k27_file = "${bam_file_name_k27[1]}"
    // rep3_k27_file = "${bam_file_name_k27[0]}"

    // rep1_k9_file = "${bam_file_name_k9[2]}"
    // rep2_k9_file = "${bam_file_name_k9[1]}"
    // rep3_k9_file = "${bam_file_name_k9[0]}"

    // doing this because it works in the for loop and now if I am given more data, it can handel it
    bam_file_list = bam_file_name.join(" ")




    """
    #!/usr/bin/env bash

    ###### Deeptools parameters ##########




    #######################################

    echo "will this work in the for loop? ${bam_file_list}"

    # making a bigwig file for each bam file
    #for bam in "\${rep1_file}" "\${rep2_file}" "\${rep3_file}"; do

    # need to find a way for this to just take how ever many bam files are put here
    # so making a list from the channel
    for bam in ${bam_file_list}; do
        
        # strip the bam file name
        bigwig_out_name="\$(basename \$bam .bam)_raw.bigwig"
        
        bamCoverage \
        --bam \$bam \
        --outFileName \$bigwig_out_name \
        --outFileFormat "bigwig" \
        --scaleFactor 1 \
        --extendReads


        # now also making the cpm bigwig

        cpm_bigwig_out_name="\$(basename \$bam .bam)_cpm_extendreads_normalized.bigwig"
        
        bamCoverage \
        --bam \$bam \
        --outFileName \$cpm_bigwig_out_name \
        --outFileFormat "bigwig" \
        --scaleFactor 1 \
        --normalizeUsing CPM \
        --extendReads

    done

    echo "for loop finished"



    """
}

// exact same process as above but for wt not control

process make_alignment_bw_process_wt {

    debug true

    label 'normal_small_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'

    publishDir "./bigwigs/${histone_label}_bigwigs/", mode: 'copy', pattern: "*${histone_label}*.bigwig"


    input:
    
    // not using multimap so i will have two instances of the process.

    tuple val(condition_label), val(histone_label), val(replicate_label), val(meta_name), val(bam_file_name), path(bam_path), path(bai_path)
    
    // if i were using multimap then i would do this below
    //tuple val(condition_label_k27), val(histone_label_k27), val(replicate_label_k27), val(bam_file_name_k27), path(bam_path_k27), path(bai_path_k27)

    //tuple val(condition_label_k9), val(histone_label_k9), val(replicate_label_k9), val(bam_file_name_k9), path(bam_path_k9), path(bai_path_k9)


    output:

    tuple val("${condition_label}"), val("${histone_label}"), val("${replicate_label}"), path("*_raw*.bigwig"), emit: bigwig_meta_ch

    // the cpm normalized bigwigs 
    tuple val("${condition_label}"), val("${histone_label}"), val("${replicate_label}"), path("*_cpm*.bigwig"), emit: cpm_bigwig_meta_ch


    script:

    // there are three  technical replicates for each label above, except histone_label 

    // rep1_file = "${bam_file_name[0]}"
    // rep2_file = "${bam_file_name[1]}"
    // rep3_file = "${bam_file_name[2]}"


    // rep1_k27_file = "${bam_file_name_k27[2]}"
    // rep2_k27_file = "${bam_file_name_k27[1]}"
    // rep3_k27_file = "${bam_file_name_k27[0]}"

    // rep1_k9_file = "${bam_file_name_k9[2]}"
    // rep2_k9_file = "${bam_file_name_k9[1]}"
    // rep3_k9_file = "${bam_file_name_k9[0]}"

    // doing this because it works in the for loop and now if I am given more data, it can handel it
    bam_file_list = bam_file_name.join(" ")



    """
    #!/usr/bin/env bash

    ###### Deeptools parameters ##########




    #######################################

    echo "will this work in the for loop? ${bam_file_list}"

    # making a bigwig file for each bam file
    #for bam in "\${rep1_file}" "\${rep2_file}" "\${rep3_file}"; do

    # need to find a way for this to just take how ever many bam files are put here
    # so making a list from the channel
    for bam in ${bam_file_list}; do
        
        # strip the bam file name
        bigwig_out_name="\$(basename \$bam .bam)_raw.bigwig"
        
        bamCoverage \
        --bam \$bam \
        --outFileName \$bigwig_out_name \
        --outFileFormat "bigwig" \
        --scaleFactor 1 \
        --extendReads


        # now also making the cpm bigwig

        cpm_bigwig_out_name="\$(basename \$bam .bam)_cpm_extendreads_normalized.bigwig"
        
        bamCoverage \
        --bam \$bam \
        --outFileName \$cpm_bigwig_out_name \
        --outFileFormat "bigwig" \
        --scaleFactor 1 \
        --normalizeUsing CPM \
        --extendReads

    done





    """
}


process plot_histone_at_genes_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools-3.5.6_rj'

    label 'normal_big_resources'
    
    //label 'super_big_resources'

    publishDir "./heatmaps", mode: 'copy', pattern: '*'


    input:

    debug true

    tuple val(condition_label), val(histone_label), val(replicate_label), val(bw_names), path(bigwig_filepath)

    //tuple val(wt_condition_label), val(wt_histone_label), val(wt_replicate_label), val(wt_bw_names), path(wt_bigwig_file)

    path(wtvslowup_genebody)

    path(wtvslowdown_nochange_genebody)



    output:

    path("${png_heatmap_upgenes}"), emit: upgene_histone_heatmap

    path("${png_heatmap_down_unchanging_genes}"), emit: down_unchanging_gene_histone_heatmap

    path("${png_heatmap_both}"), emit: both_histone_heatmap
    //path("${png_wt_heatmap}"), emit: gene_histone_heatmap_wt




    script:

    name_list = []

    num_files = bw_names.size()
    
    for (int i = 0; i < num_files; i++) {
        true_basename = "${bw_names[i]}".replaceFirst(/\..*/, '')
        name_list << true_basename

    }

    //newName = bw_names.replaceFirst(/\..*/, '')
    //list_control_bw_names = "${control_bw_names.toList()}"

    out_matrix_scores_upgenes = "matrix_gene_${histone_label}_lowup_genebody.mat.gz"
    //out_matrix_scores_wt = "matrix_gene_${histone_label}_lowup_genebody.mat.gz"

    out_matrix_scores_down_unchange_genes = "matrix_gene_${histone_label}_lowdown_unchanging_genebody.mat.gz"

    out_matrix_scores_both = "matrix_gene_${histone_label}_up_and_down_unchanging_genebody.mat.gz"

    png_heatmap_upgenes = "${histone_label}_histone_features_at_lowup_genebody.png"
    png_heatmap_down_unchanging_genes = "${histone_label}_histone_features_at_lowdown_unchanging_genebody.png"
    png_heatmap_both = "${histone_label}_histone_features_at_up_and_down_unchanging_genebody.png"

    //png_wt_heatmap = "${wt_histone_label}_histone_features_at_lowup_genebody.png"

    // making output name for genes file that will have no zeros

    up_genes_nozero = "${wtvslowup_genebody.baseName}_noZerolength.bed"
    down_unchanging_genes_nozero = "${wtvslowdown_nochange_genebody.baseName}_noZerolength.bed"

    """
    #!/usr/bin/env bash

    ########### deeptools params ##########



    #######################################

    #echo ' this is the list of names: "\${name_list}"'

    # fix the genebody file
    awk  '\$2!=\$3 {print \$0}' "${wtvslowup_genebody}" > "${up_genes_nozero}"

    computeMatrix scale-regions -S ${bw_names.join(' ')} \
    -R "${up_genes_nozero}" \
    --beforeRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 3000 \
    --skipZeros \
    --quiet \
    --numberOfProcessors "max" \
    -o "${out_matrix_scores_upgenes}"

    
    plotHeatmap -m "${out_matrix_scores_upgenes}" \
    -out "${png_heatmap_upgenes}" \
    --samplesLabel ${name_list.join(' ')} \
    --labelRotation 30 \
    --heatmapWidth 8 \
    --dpi 300 \
    --sortUsing sum



    # now to do this with the other genebody file 'down unchanging'

    awk  '\$2!=\$3 {print \$0}' "${wtvslowdown_nochange_genebody}" > "${down_unchanging_genes_nozero}"

    computeMatrix scale-regions -S ${bw_names.join(' ')} \
    -R "${down_unchanging_genes_nozero}" \
    --beforeRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 3000 \
    --skipZeros \
    --quiet \
    --numberOfProcessors "max" \
    -o "${out_matrix_scores_down_unchange_genes}"

    
    plotHeatmap -m "${out_matrix_scores_down_unchange_genes}" \
    -out "${png_heatmap_down_unchanging_genes}" \
    --samplesLabel ${name_list.join(' ')} \
    --labelRotation 30 \
    --heatmapWidth 8 \
    --dpi 300 \
    --sortUsing sum


    echo ' the up genes and down unchanging genes plot has completed. Now starting on the plot for both up and down genes together'


    # now plotting both together

    computeMatrix scale-regions -S ${bw_names.join(' ')} \
    -R "${up_genes_nozero}" "${down_unchanging_genes_nozero}" \
    --beforeRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 3000 \
    --skipZeros \
    --quiet \
    --numberOfProcessors "max" \
    -o "${out_matrix_scores_both}"

    plotHeatmap -m "${out_matrix_scores_both}" \
    -out "${png_heatmap_both}" \
    --samplesLabel ${name_list.join(' ')} \
    --labelRotation 30 \
    --heatmapWidth 8 \
    --dpi 300 \
    --sortUsing sum







    """
}




process macs2_call_peaks_process_both {

    //debug true
    errorStrategy 'ignore'

    label 'normal_big_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/macs2_rj'

    publishDir "./peak_files/${condition_label}", mode: 'copy', pattern: '*'



    input:

    tuple val(condition_label), val(histone_label), val(replicate_label), val(meta_name), val(bam_file_name), path(bam_path), path(bai_path)

    path(ref_genome)
    // condition_label has the same condition multiple times (3 in this run)
    // the replicate_label had multiple also (3) but not the same (r1,r2,r3)


    output:
    path("*broad*"), emit: broadpeaks
    path("*"), emit: macs2_files
    path("*ppois*"), emit: ppois_macs2_file



    script:

    //num_files = bam_file_name.size()
    
    old_peak_files = "${bam_file_name}_old_macs2_stats"
    

        
    """
    #!/usr/bin/env bash

    ##### macs2 params ######
    # use this code to find params used
    # macs2 callpeak --help

    # changed --fe-cutoff from 2 to 1 which is default. to see if it will give results in the bio rep s10 file
    # that didn't change anything

    # I want to call peaks less stringently so we can bring a bit more noise into the peaks that were called and then idr will find the best ones
    # so changing --fe-cutoff from 2 to 1, but i will remove it altogether
    # changing -qvalue from '0.05' to '0.01' probably not using it anymore
    # using --pvalue 1e-3 
    #########################

    #echo "these are the file names: \${bam_file_name}"

    #bam_list=( \${bam_file_name.join(' ')} )

        

    #bam_basename=\$(basename \${bam_list[i]})


    macs2 callpeak \
    --treatment ${bam_file_name} \
    --format "BAM" \
    --gsize "hs" \
    --keep-dup '1' \
    --outdir . \
    --name ${bam_file_name} \
    --bdg  \
    --trackline \
    --pvalue '1e-3' \
    --broad \
    --cutoff-analysis \
    --nolambda

    # first compute the sval score then use it in macs2 bdgcmp
    chipReads=\$(cat "${bam_file_name}_peaks.broadPeak"| wc -l | awk '{printf "%f", \$1/1000000}')
    
    # since we do not have a control we do not have to get the control reads and compare it with the chipReads(C&T reads) to find which is the lower reads and use that as sval
    # so we will just use the chipReads as sval
    sval=\$(echo "\$chipReads")

    echo "sval when creating pvalue bigwig from macs2 is \$sval"


    # now we need to create the signal pvalue bigwig file. (according to the encode chip-seq pipeline google doc: https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit?tab=t.0)

    macs2 bdgcmp \
    -t *_treat_pileup.bdg \
    -c *_control_lambda.bdg \
    -o ${bam_file_name}_ppois.bdg \
    -m ppois \
    -S \$sval




    # this is me trying to get the old peak files that produce the correct signal in peaks
    #macs2 callpeak \
    --treatment \${bam_file_name} \
    --format "BAM" \
    --gsize "hs" \
    --keep-dup '1' \
    --outdir . \
    --name \${old_peak_files} \
    --bdg  \
    --trackline \
    --qvalue '0.05' \
    --broad \
    --fe-cutoff 2 \
    --cutoff-analysis
    
    # use bedtools merge, in another process, on all the broadpeak files  and the option -b for merging at 1kb, 2kb, 5kb. generate different merged files for the same data to compare

    # I need to remove nolambda because scrm r1 only is calling 7 peaks and the IDR tool needs at least 20 peaks to run
    # this is according to the IDR error output


    """
    

    
}

process get_pval_bedgraph {
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'
    label 'normal_big_resources'
    publishDir "./peak_files/pval_signal_files", mode: 'copy', pattern: '*'


    input:

    path(ppois_file)

    path(genome_chr_size)


    output:

    path("${bedgraph_file_out}"), emit: pvalue_bedgraph_file
    path("genome.bed"), emit: chrom_size_file


    script:

    bedgraph_file_out = "${ppois_file.baseName}.pval.signal.bedgraph"

    """
    #!/usr/bin/env bash

    # making the proper genome.bed file from the chromosome size file generated. the genome.fa.fai
    awk 'BEGIN{OFS="\t"}{ print \$1, 0, \$2 }' "${genome_chr_size}" > genome.bed

    bedtools slop \
    -i ${ppois_file} \
    -g ${genome_chr_size} \
    -b 0 \
    | bedtools intersect \
        -a stdin \
        -b genome.bed \
    | sort -k1,1 -k2,2n \
    > ${bedgraph_file_out}





    """
}

process kenttools_get_bigwig {
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/uscs_utils_rj'
    label 'normal_big_resources'
    publishDir "./peak_files/pval_signal_files", mode: 'copy', pattern: '*'


    input:
    path(bedgraph_files)
    path(chrom_size_file)

    output:
    path("${out_pval_signal_bigwig}"), emit: pval_signal_bigwig_files

    script:

    out_pval_signal_bigwig = "${bedgraph_files.baseName}.bigwig"

    """
    #!/usr/bin/env bash

    bedGraphToBigWig \
    ${bedgraph_files} \
    ${chrom_size_file} \
    ${out_pval_signal_bigwig}





    """

}

process merge_peaks_bedtools_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'

    label 'normal_big_resources'

    publishDir "merged_broadpeaks/", mode: 'copy', pattern:'*'


    input:
    path(peakpath)


    output:
    path("*_merged.bed"), emit: all_merged_broadpeaks



    script:

    peak_basename = peakpath.baseName
    peak_filename = peakpath.name

    merged_1kb = "${peak_basename}_1kb_merged.bed"
    merged_2kb = "${peak_basename}_2kb_merged.bed"
    merged_5kb = "${peak_basename}_5kb_merged.bed"

    """
    #!/usr/bin/env bash

    ###### bedtools merge parameters to use ##########


    ##################################################

    # doing 1kb first
    bedtools merge \
    -i ${peak_filename} \
    -d 1000 \
    -c 1 \
    -o count \
    > ${merged_1kb}

    # doing 2kb first
    bedtools merge \
    -i ${peak_filename} \
    -d 2000 \
    -c 1 \
    -o count \
    > ${merged_2kb}

    # doing 5kb first
    bedtools merge \
    -i ${peak_filename} \
    -d 5000 \
    -c 1 \
    -o count \
    > ${merged_5kb}




    """






}

process merge_concat_peaks_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'

    publishDir "idr_results/merged_concat_peaks", mode: 'copy', pattern: "*.bed"

    label 'normal_big_resources'

    errorStrategy 'ignore'

    debug true


    input:

    tuple val(histone), val(condition), val(concat_peaks_names), path(concat_peaks_paths)



    output:

    //tuple val(histone), val(condition), val(list_of_10kb_merged_names), val(list_of_10kb_merged_paths) , emit: merged_10kb_concat_peaks

    path("${output_10kb_merged_first_mPeak}"), emit: first_10kb_merged_peak

    path("${output_10kb_merged_second_mPeak}"), emit: second_10kb_merged_peak

    path("${output_30kb_merged_first_mPeak}"), emit: first_30kb_merged_peak

    path("${output_30kb_merged_second_mPeak}"), emit: second_30kb_merged_peak

    path("${output_100kb_merged_first_mPeak}"), emit: first_100kb_merged_peak

    path("${output_100kb_merged_second_mPeak}"), emit: second_100kb_merged_peak

    path("${concat_10kb_merged_both_mPeak}"), emit: concat_10kb_merged_peak

    path("${concat_30kb_merged_both_mPeak}"), emit: concat_30kb_merged_peak

    path("${concat_100kb_merged_both_mPeak}"), emit: concat_100kb_merged_peak



    script:



    first_mPeak = concat_peaks_names[0]
    second_mPeak = concat_peaks_names[1]

    println(first_mPeak)
    println(second_mPeak)

    
    
    output_10kb_merged_first_mPeak = first_mPeak.replace(/.broadPeak/, "_10kb_merged.bed")
    output_10kb_merged_second_mPeak = second_mPeak.replace(/.broadPeak/, "_10kb_merged.bed")
    concat_10kb_merged_both_mPeak = "concat_master_bothConditions_${histone}_10kb_merged.bed"

    output_30kb_merged_first_mPeak = first_mPeak.replace(/.broadPeak/, "_30kb_merged.bed")
    output_30kb_merged_second_mPeak = second_mPeak.replace(/.broadPeak/, "_30kb_merged.bed")
    concat_30kb_merged_both_mPeak = "concat_master_bothConditions_${histone}_30kb_merged.bed"

    output_100kb_merged_first_mPeak = first_mPeak.replace(/.broadPeak/, "_100kb_merged.bed")
    output_100kb_merged_second_mPeak = second_mPeak.replace(/.broadPeak/, "_100kb_merged.bed")
    concat_100kb_merged_both_mPeak = "concat_master_bothConditions_${histone}_100kb_merged.bed"


    //list_of_10kb_merged_names =  [output_10kb_merged_first_mPeak, output_10kb_merged_second_mPeak ]
    //list_of_10kb_merged_paths =  [output_10kb_merged_first_mPeak, output_10kb_merged_second_mPeak ]

    """
    #!/usr/bin/env bash

    ###### bedtools merge parameters to use ##########
    # the -c parameter: allows you to choose a column and determine what will happen when combined with -o parameter
    # the -o parameter: you can get the counts of how many regions were merged to create that one region. "count" will report the the sum of merged intervals

    ##################################################

    # have to sort the broadPeak files first
    sort -k1,1 -k2,2n ${first_mPeak} > sorted_first_peak.broadPeak

    sort -k1,1 -k2,2n ${second_mPeak} > sorted_second_peak.broadPeak

    

    # merging the first concat peaks by 10kb
    bedtools merge \
    -i sorted_first_peak.broadPeak \
    -d 10000 \
    > ${output_10kb_merged_first_mPeak}


    # merging the second concat peaks by 10kb
    bedtools merge \
    -i sorted_second_peak.broadPeak \
    -d 10000 \
    > ${output_10kb_merged_second_mPeak}

    # concatenating 10kb
    cat ${output_10kb_merged_first_mPeak} ${output_10kb_merged_second_mPeak} > ${concat_10kb_merged_both_mPeak}


    # merging the first concat peaks by 30kb
    bedtools merge \
    -i sorted_first_peak.broadPeak \
    -d 30000 \
    > ${output_30kb_merged_first_mPeak}


    # merging the second concat peaks by 30kb
    bedtools merge \
    -i sorted_second_peak.broadPeak \
    -d 30000 \
    > ${output_30kb_merged_second_mPeak}

    # concatenating 30kb
    cat ${output_30kb_merged_first_mPeak} ${output_30kb_merged_second_mPeak} > ${concat_30kb_merged_both_mPeak}


    # merging the first concat peaks by 100kb
    bedtools merge \
    -i sorted_first_peak.broadPeak \
    -d 100000 \
    > ${output_100kb_merged_first_mPeak}


    # merging the second concat peaks by 100kb
    bedtools merge \
    -i sorted_second_peak.broadPeak \
    -d 100000 \
    > ${output_100kb_merged_second_mPeak}

    # concatenating 100kb
    cat ${output_100kb_merged_first_mPeak} ${output_100kb_merged_second_mPeak} > ${concat_100kb_merged_both_mPeak}


    """
}


process find_idr_in_replicates_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/idr-2.0_rj'
    label 'super_big_resources'
    publishDir "idr_results/${histone[0]}/${condition[0]}", mode: 'copy', pattern:'*'

    // because of some experiment types (like histone H3k36me2) not having enough peaks called in all replicates, we cannot have this process ending the entire pipeline just becasue of one experiment type
    // so i will add ignore errorStrategy to this process also
    errorStrategy 'ignore'


    input:

    tuple val(grouping_key), val(condition), val(histone), val(replicate), val(bio_rep), val(file_name), val(basename), path(peakpath)
    // for the file_name there will be three replicates, all in order. this will be the case for each instance this process is called when parallelized
    // idr takes only 2 replicates at a time
    // we will do rep 1 and 2, then rep 2 and 3, then rep 1/2 and 2/3
    // this process will only work for if we have three replicates in our data


    output:
    // now i need to output the final idr broadpeak and then get all the others for this experiment histone mark and combine them
    // I can create another process that will filter these using the blacklist and bedtools intersect but lets look at it as is for now. 
    path("*broadPeak"), emit: idr_peaks

    path("*.png"), emit:idr_pngs

    path("${concat_peaks_name}"), emit: final_concat_peaks

    script:

    // the conditions is always the same but it has three in the val above. ex: [Hlow, Hlow, Hlow]
    condition_label = condition[0] // so i can just use one of them to rebuild a file name if needed

    // same with the histone ex: [H3k27me3, H3k27me3, H3k27me3]
    histone_label = histone[0]

    // the replicates are in order but not the same ex: [r1, r2, r3]
    // but i can still store the correct thing
    rep_label1 = replicate[0]
    rep_label2 = replicate[1]
    rep_label3 = replicate[2]

    // same with bio_reps, they are in order and correspond properly with the order of rep_labels
    bio_label1 = bio_rep[0]
    bio_label2 = bio_rep[1]
    bio_label3 = bio_rep[2]

    // these will be the broadpeak file names, the paths are already in the process directory so i can use the correct name in order. the paths will not be in order
    peak1 = file_name[0]
    peak2 = file_name[1]
    peak3 = file_name[2]

    // making the sorted filename
    sort_peak1 = "sorted_${peak1}.gz"
    sort_peak2 = "sorted_${peak2}.gz"
    sort_peak3 = "sorted_${peak3}.gz"

    // now to create the output names for checking each replicate

    idr_out_1_2_name = "IDR_${condition_label}_${histone_label}_${rep_label1}_vs_${rep_label2}_${params.return_idr}.broadPeak"
    idr_out_2_3_name = "IDR_${condition_label}_${histone_label}_${rep_label2}_vs_${rep_label3}_${params.return_idr}.broadPeak"
    idr_out_1_3_name = "IDR_${condition_label}_${histone_label}_${rep_label1}_vs_${rep_label3}_${params.return_idr}.broadPeak"
    idr_out_final_name = "IDR_${condition_label}_${histone_label}_${rep_label1}_vs_${rep_label2}_vs_${rep_label3}_${params.return_idr}.broadPeak"
    idr_out_final_merged_name = "IDR_${condition_label}_${histone_label}_${rep_label1}_vs_${rep_label2}_vs_${rep_label3}_2_merged_${params.return_idr}.broadPeak"

    concat_peaks_name = "concat_IDR_${condition_label}_${histone_label}_${rep_label1}_vs_${rep_label2}_vs_${rep_label3}_${params.return_idr}_pairs.broadPeak"

    // name for pooled broad peak files 
    broadpeak_pool = "broadpeak_pool.broadPeak"
    sort_broadpeak_pool = "sort_${broadpeak_pool}.gz"

    // now creating the file name for the idr's of all the reps

    // idr_final_sorted_file = "${condition_label}_${histone_label}_${rep_label1}_vs_${rep_label2}_vs_${rep_label3}_IDR_${params.return_idr}.broadPeak.gz"

    """
    #!/usr/bin/env bash

    ###### idr parameters ########
    # plot_idr = 0.05 // the default is 0.05 to report in the png plots which peaks passed this threshold
    # return_idr = 1  // the default is all peaks will be returned even if the report plots the ones that pass a certain number. 1 as default here will give back all peaks like idr has set

    #--peak-list is not provided
    #Peaks are grouped by overlap and then merged. The merged peak aggregate value is determined by --peak-merge-method.

    #Peaks that don't overlap another peak in every other replicate are not included unless --use-nonoverlapping-peaks is set.
    ##############################

    # i need to pool the 3 reps so i can have the pooled file for input in idr
    cat ${peak1} ${peak2} ${peak3} > ${broadpeak_pool}

    # i have to gzip the file
    #gzip -nc \${peak1} > "\${peak1}.gz"
    #gzip -nc \${peak2} > "\${peak2}.gz"
    #gzip -nc \${peak3} > "\${peak3}.gz"


    # now sorting all the peaks by their pvalue column
    sort -k8,8nr "${peak1}" | gzip > ${sort_peak1}
    sort -k8,8nr "${peak2}" | gzip > ${sort_peak2}
    sort -k8,8nr "${peak3}"| gzip > ${sort_peak3}
    sort -k8,8nr ${broadpeak_pool} | gzip > "${sort_broadpeak_pool}"

    # now i have everything to put into idr
    # i think its best to return all peaks that passed the merging criteria for rep1 vs rep2 and rep2 vs rep3
    # then when i do rep1_2 vs rep2_3 I only return the peaks that have the 0.05 IDR threshold passed; these peaks will also be the ones that pass the merging threshold.
    # you can change this idr by changing the value in the --return_idr parameter when running nextflow

    # i will make a if then statement where i will choose to run the idr if the length of the peak file is larger than 21 lines
    peak1_length=\$(less ${sort_peak1} | wc -l )
    peak2_length=\$(less ${sort_peak2} | wc -l )
    peak3_length=\$(less ${sort_peak3} | wc -l )

    # rep 1 vs 2
    
    if ((\$peak1_length > 21 && \$peak2_length > 21)); then
        idr --samples ${sort_peak1} ${sort_peak2} \
        --input-file-type broadPeak \
        --output-file ${idr_out_1_2_name} \
        --rank p.value \
        --idr-threshold ${params.return_idr} \
        --soft-idr-threshold ${params.plot_idr} \
        --plot \
        --use-best-multisummit-IDR
    else
        touch ${idr_out_1_2_name}
        touch "place_holder.png"
    fi 

    # now rep 2 vs 3

    if ((\$peak2_length > 21 && \$peak3_length > 21)); then
        idr --samples ${sort_peak2} ${sort_peak3} \
        --input-file-type broadPeak \
        --output-file ${idr_out_2_3_name} \
        --rank p.value \
        --idr-threshold ${params.return_idr} \
        --soft-idr-threshold ${params.plot_idr} \
        --plot \
        --use-best-multisummit-IDR
    else
        touch ${idr_out_2_3_name}
        touch "place_holder.png"
    fi 

    # now doing 1 vs 3

    if ((\$peak1_length > 21 && \$peak3_length > 21)); then
        idr --samples ${sort_peak1} ${sort_peak3} \
        --input-file-type broadPeak \
        --output-file ${idr_out_1_3_name} \
        --rank p.value \
        --idr-threshold ${params.return_idr} \
        --soft-idr-threshold ${params.plot_idr} \
        --plot \
        --use-best-multisummit-IDR
    else
        touch ${idr_out_1_3_name}
        touch "place_holder.png"
    fi 

    # maybe we dont look at the final output and just select the output of peak pairs that has the most peaks that pass the 0.05 threshold
    # this is according to section 4d of the encode 3 pipeline 'https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit?tab=t.0#heading=h.9ecc41kilcvq' 
    # I will instead load all the files into R and use only the one from each condition that has the max number of peaks called
    # then I will merge the two files; the max from hlow and the max from scrm to get the masterpeak
    
    #idr_1_2_count=\$(less \${idr_out_1_2_name} | wc -l)
    #idr_1_3_count=\$(less \${idr_out_1_3_name} | wc -l)
    #idr_2_3_count=\$(less \${idr_out_2_3_name} | wc -l)


        


    # now the final output
    # this would be merging according to idr standards, keeping all the peaks from the 2 of the 3 pairs
    
    if ((\$peak1_length > 21 && \$peak2_length > 21 && \$peak3_length > 21)); then
        idr --samples ${idr_out_1_2_name} ${idr_out_2_3_name} \
        --input-file-type broadPeak \
        --output-file ${idr_out_final_merged_name} \
        --rank p.value \
        --plot \
        --use-best-multisummit-IDR
    fi


    # just concatenate them and see
    
    cat ${idr_out_1_2_name} ${idr_out_2_3_name} ${idr_out_1_3_name} > ${concat_peaks_name}
    



    #if ((\$peak1_length > 21 && \$peak2_length > 21 && \$peak3_length > 21)); then
        #idr --samples \${idr_out_1_2_name} \${idr_out_2_3_name} \
        --input-file-type broadPeak \
        --output-file \${idr_out_final_name} \
        --idr-threshold \${params.return_idr} \
        --soft-idr-threshold \${params.plot_idr} \
        --rank p.value \
        --plot \
        --use-best-multisummit-IDR
    #fi


    # will have to ask johanna what this next bit of code does, but i can convert it to work here in nextflow as i did above
    # it checks to see if column 12 is greater than or equal to the idr_thresh_transformed, if so print all the columns

    # dont need this since it might be using the wrong column and IDR already has a parameter to get only peaks passing the threshold
    
    # IDR_THRESH_TRANSFORMED=\$(awk -v p=0.05 'BEGIN{print -log(p)/log(10)}')

    # can change this to print 0 instead of listing all of them in. that will print all the columns
    
    # if you read the documentation, broad peak outputs do not have a summit columns (10,18,22 for 2 replicates), so that means i am actually using column 11 not 12 for broad peaks not narrow peaks
    #awk 'BEGIN{OFS="\t"} \$11>='"\${IDR_THRESH_TRANSFORMED}"' {print \$0}' \${idr_out_final_name} | \
    sort | uniq | sort -k7n,7n | gzip -nc > \${idr_final_sorted_file}



    """



}



process find_diff_peaks_R_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/r_language'

    publishDir "./nextflow_R_script_outputs/${idr_histone}/", mode: 'copy', pattern: '*'

    label 'normal_big_resources'

    //debug true

    input:

    tuple val(idr_condition), val(idr_histone), val(merge_dist), val(idr_peak_name), path(idr_peak_path)

    // these are all experiment bams.
    // i need to use the histone label, and condition label to get the correct bams in this process
    path(bam_path_list)
    


    
    output:

    path("*.{png,pdf,bed}"), emit: all_r_plots

    path("${master_peak_export_out}"), emit: master_peak_emit

    tuple val("${full_condition}"), val("${idr_histone}"), path("${up_peaks_out}"), path("${down_peaks_out}"), path("${unchanging_peaks_out}"), emit: diff_peaks_ch 

    path("${up_peaks_out}"), emit: up_peaks_emit
    path("${down_peaks_out}"), emit: down_peaks_emit
    path("${unchanging_peaks_out}"), emit: unchanging_peaks_emit



    script:

    //full_condition = "${idr_condition[0]}vs${idr_condition[1]}"
    full_condition = "${idr_condition}"

    //first_idr_peak = idr_peak_name[0]
    //second_idr_peak = idr_peak_name[1]

    bam_name_list = bam_path_list.collect { "\"${it.getName()}\"" }.join(',')

    //peaks_list = [first_idr_peak, second_idr_peak]

    // need the output file name for masterpeak export

    master_peak_export_out = "masterpeak_${idr_histone}_${merge_dist}_maxgap_${idr_condition}_${idr_histone}.bed"

    up_peaks_out = "up_${idr_histone}_regulated_peaks.bed"
    down_peaks_out = "down_${idr_histone}_regulated_peaks.bed"
    unchanging_peaks_out = "unchanging_${idr_histone}_regulated_peaks.bed"


    """
    #!/usr/bin/env Rscript

    bam_list = c(${bam_name_list})

    print(bam_list)

    #peaklist = list("\${first_idr_peak}", "\${second_idr_peak}")
    #peaklist = \${peaks_list}

    #print(peaklist)

    #best_hlow_idr
    #best_scrm_idr

    library(DESeq2)
    library(GenomicRanges)
    library(chromVAR)
    library(tidyr)
    library(EnhancedVolcano)
    library(readr)

    print(dir())
    # making the master peak genomic ranges object
    
    #mPeak = GRanges()

    #for (peakfile in c("./\${first_idr_peak}", "./\${second_idr_peak}")) {
    
        #peaktable = read.table(peakfile, header = FALSE, sep = "\t")
    
        #gr_object = GRanges(seqnames = peaktable\$V1, IRanges(start = peaktable\$V2, end = peaktable\$V3), strand = "*" )
    
        #gr_object = rtracklayer::import(peakfile, format = "BED")

        #mPeak = append(mPeak, gr_object)
    #}

    #peaktable = read.table("\${idr_peak_name}", header = FALSE)
    #mPeak = GRanges(seqnames = peaktable\$V1, IRanges(start = peaktable\$V2, end = peaktable\$V3), strand = "*" )

    mPeak = rtracklayer::import("${idr_peak_name}", format = "BED")

    # making sure there are no redundant peaks
    # this will allow me to resizd all of the regions by 50kb on each side, so any regions that are less than 100kb away from eachother will be merged by reduce because they overlap
    
    
    #masterPeak2 = reduce(mPeak)

    # using resize to extnd regions
    extended_master_peak = resize(mPeak, width = width(mPeak)+100000, fix = "center")

    # before trimming i need to add the sequence lengths of the correct genome
    # will have to find a way to automate this step for other genomes

    library(BSgenome.Hsapiens.UCSC.hg38)
    #seq_lengths <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)
    #seqlengths(extended_master_peak) <- seq_lengths[seqlevels(extended_master_peak)]
    #trimmed_master_peak = trim(extended_master_peak)
    #masterPeak_beta = reduce(trimmed_master_peak)

    masterPeak_beta = reduce(extended_master_peak)
    # now to resize again to remove the artificial 50kb from the start and end of the newly reduced peak
    
    new_peak_size = width(masterPeak_beta)-100000
    masterPeak = resize(masterPeak_beta,  width = new_peak_size , fix = "center")
    
    #print(masterPeak)

    # now I want to keep the standard chromosomes
    #seqnames_to_keep =masterPeak@seqnames@values[1:23]

    #masterPeak = masterPeak[seqnames(masterPeak) %in% seqnames_to_keep]

    # hoping to export my GRanges object master peaks to a bed file
    #write_tsv("\${master_peak_export_out}")
    rtracklayer::export.bed(masterPeak, con = "${master_peak_export_out}")

    

    list_of_mBams = c()
    
    for (bam in bam_list) {
  
        tokens = strsplit(bam, split = "_")[[1]]
        #print(tokens)
        
        histone = tokens[2]
        #print(histone)

        #list_of_mBams = list()

        if ("${idr_histone}" == histone) {
            #print(bam)
            list_of_mBams = c(list_of_mBams, bam)

        }
    }

    print(list_of_mBams)
   
    
    # now doing the next section to make the count matrix and the condition and experiment design

    # making the matrix 

    countsMatrix = matrix(data = NA, length(masterPeak), length(list_of_mBams) )


    # I want to get the chr start end of the peaks and put them in the matrix as column names so I know which peaks I am looking at from the deseq2 results in the differential peaks.
    seqnames(masterPeak)

    # i should put the unique ids from the master peak object in the matrix as row names.
    unique_masterPeak_ids = paste0(as.character(seqnames(masterPeak)), ":" ,start(masterPeak), "-", end(masterPeak))

    rownames(countsMatrix) = unique_masterPeak_ids

    #getting the list of bam base names to add to the matrix column names
    list_bam_basenames = list()


    # for deseq2 i need the condition design. so hlow and scrm
    # then i will find a way to tally how many of each are there so i can automate the rep count
    condition_design = list()

    type_design = list()

    for (x in c(1:length(list_of_mBams))) {
    
    path_bam = list_of_mBams[[x]]
    print(path_bam)
    bam_basename = basename(path_bam)
    bam_tokens = strsplit(basename(path_bam), split = "_")[[1]]
    
    # labeling the important tokens so it is easier to keep track of
    condition = bam_tokens[1]
    histone = bam_tokens[2]
    replicate = bam_tokens[3]
    
    
    # for later parts I need the condition and to know how many times to repeat them
    condition_design = append(condition_design, paste(condition,histone,sep="_"))
    
    
    # also get the replicate too for type design
    type_design = append(type_design, replicate)
    #type_design = append(type_design, paste(histone,replicate, sep="_"))
    
    
    
    # using chromVAR getCounts function
    fragment_counts = getCounts(path_bam, masterPeak, paired= TRUE, by_rg = FALSE, format = "bam")
    
    
    # putting the fragment counts in the column labeled by the bam name
    countsMatrix[,x] = counts(fragment_counts)[,1]
    
    list_bam_basenames = append(list_bam_basenames, bam_basename)
    
    }

    colnames(countsMatrix) = list_bam_basenames

    ## this part i would have to tell the use to name their conditions as control and treatment so it can make the design treatment vs control

    #first removing the low count fragments. only keeping the rows that are above 5 count in total

    keep_Rows = which(rowSums(countsMatrix) > 5)

    filt_countmatrix = countsMatrix[keep_Rows,]

    # now to get the condition_design
    condition_counts = table(unlist(condition_design))

    # this gives back the names and the counts give back the counts for each of the names
    condition_names = names(condition_counts)
    condition_num = as.numeric(condition_counts)




    # now i can put both lists in to get back the experiment design
    condition_factor = factor(rep(condition_names, times=condition_num))


    # I want it so the treatment is second in the list here so it is first in the experimental design in deseq2. This way the experimental design will have the diff results going up or down in the treatment vs control

    # for treatment I need to make nextflow take the users input so they can relevel this section
    #treatment = "H1low_H3k27me3"
    treatment = "${params.treatment_name}_${idr_histone}"

    if (levels(condition_factor)[1] == treatment ) {
    condition_factor = relevel(condition_factor, levels(condition_factor)[2])
    }else {
    condition_factor
    }
    print(condition_factor)

    # repeating the above to have another column with type (replicates)
    type_counts = table(unlist(type_design))
    type_names = names(type_counts)
    type_num = as.numeric(type_counts)

    #type_factor = factor(rep(type_names, times=type_counts))
    type_factor = factor(rep(type_names, times=type_num[1]))

    # I want to get the idr threshold and use that as input for the file names and other things

    #peak_file = basename(peaklist[[1]])
    peak_file = basename("./${idr_peak_name}")

    idr_used = strsplit(peak_file, split = "_")[[1]][10]
    idr_used

    # now for deseq2 workflow


    # now to do the normal deseq2 workflow

    dds = DESeqDataSetFromMatrix(countData = filt_countmatrix,
                                colData = DataFrame(condition_factor, type_factor),
                                design = ~ condition_factor)


    # using the function on our data
    DDS = DESeq(dds)

    norm_DDS = counts(DDS, normalized = TRUE) # normalization with respect to the sequencing depth

    # adding _norm onto the column names in the normalized matrix
    colnames(norm_DDS) = paste0(colnames(norm_DDS), "_norm")


    # provides independent filtering using the mean of normalized counts
    res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")


    # this is looking at the differences between the 3 deseq analyzed options
    countMatDiff = cbind(filt_countmatrix, norm_DDS, res)

    head(countMatDiff)




    # getting the results name and addding to the coef we want to shrink
    experiment_design_name = resultsNames(DDS)[2]

    # useful for visualization and ranking of genes or in this case peaks
    resLFC = lfcShrink(DDS, coef= resultsNames(DDS)[2], type = "apeglm")


    # finding the up and down regulated counts that pass the threshold

    up_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange >= 1) ,]
    total_up_reg = length(up_reg[,1])
    #total_up_reg
    down_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange <= -1) ,]
    total_down_reg = length(down_reg[,1])
    total_down_reg

    unchanging_reg = resLFC[which(resLFC\$padj < 0.05 & resLFC\$log2FoldChange <1 & resLFC\$log2FoldChange > -1) ,]
    total_unchanging_reg = length(unchanging_reg[,1])
    total_unchanging_reg

    resLFC\$Label = ifelse(resLFC\$padj > 0.05, "not significant", ifelse(abs(resLFC\$log2FoldChange) >=1, "|LFC| > 1 & padj < 0.05", "padj < 0.05" ))
    ma_plot_labeled = ggplot(resLFC, aes(x = baseMean, y = log2FoldChange, color=Label))+
    labs(caption = paste("Up reg (green +1 LFC & padj <0.05) = ", total_up_reg, "\n", "Down reg (green -1 LFC & padj < 0.05) = ", total_down_reg), title = experiment_design_name, subtitle = paste("This means the first condition has these points more than in the second condition. The idr threshold for masterPeaks is",idr_used, sep = " "))+
    theme(plot.caption = element_text(size = 25, face = "bold"), axis.text.x = element_text(size = 25), axis.text.y = element_text(size = 25), axis.title.x = element_text(size = 25), axis.title.y = element_text(size = 25), legend.text = element_text(size = 25))+
    scale_colour_manual(values=c("red", "darkgrey", "pink"))+
    scale_x_continuous(trans='log10')+
    ylim(c(min(-max(resLFC\$log2FoldChange),min(resLFC\$log2FoldChange)), max(max(resLFC\$log2FoldChange),-min(resLFC\$log2FoldChange))))+
    geom_point(data = resLFC, aes(x = baseMean, y = log2FoldChange), size = 3)

    print(ma_plot_labeled)


    pdf(file = paste(experiment_design_name,"IDR", idr_used,"MA_plot_our_data_counts.pdf", sep = "_"), width = 12, height = 12)

    print(ma_plot_labeled)
    dev.off()

    png(filename = paste(experiment_design_name,"IDR", idr_used,"MA_plot_our_data_counts.png", sep = "_"), width = 900, height = 900, antialias = "subpixel")

    print(ma_plot_labeled)
    dev.off()


    volcano_plot_removed_reps = EnhancedVolcano(resLFC,
                    lab = rownames(resLFC),
                    title = experiment_design_name, subtitle = paste("This means the first condition has these points more than in the second condition. The idr threshold for masterPeaks is",idr_used, sep = " "),
                    x = 'log2FoldChange', FCcutoff = 1,
                    y = 'padj', pCutoff = 0.05, labSize = 0 # making this 0 so it doesn't show the region names in the volcano plot 
                    
                    )

    print(volcano_plot_removed_reps)

    png(filename = paste(experiment_design_name,"IDR",idr_used,"volcano_plot_our_data_counts.png", sep = "_"), width = 900, height = 900, antialias = "subpixel")
    print(volcano_plot_removed_reps)
    dev.off()

    pdf(file = paste(experiment_design_name,"IDR",idr_used,"volcano_plot_our_data_counts.pdf", sep = "_"), width = 12, height = 12)
    print(volcano_plot_removed_reps)
    dev.off()
    



    # using rlog over vst for transformation

    rld = rlog(DDS, blind=FALSE)

    #it is clear that the rlog transformation inherently accounts for differences in sequencing depth *
    head(assay(rld), 5)


    library(ggplot2)

    pcaData = plotPCA(rld, intgroup=c("condition_factor","type_factor"), returnData = TRUE )
    pcaData

    percentVar <- round(100 * attr(pcaData, "percentVar"))


    pca_plot_rlog = ggplot(pcaData, aes(PC1, PC2, color=condition_factor, shape=type_factor)) +
    ggtitle(paste("PCA plot using Rlog transform IDR", idr_used, sep = " ")) +
    theme(plot.caption = element_text(size = 30, face = "bold"), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15))+
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
    #pca_plot
    name_for_rlog_png_file =paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used, "rlog.png", sep = "_")

    png(filename = name_for_rlog_png_file, width = 900, height = 900, antialias = "subpixel")
    print(pca_plot_rlog)
    dev.off()

    name_for_rlog_pdf_file =paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used, "rlog.pdf", sep = "_")

    pdf(file = name_for_rlog_pdf_file, width = 10, height = 10)
    print(pca_plot_rlog)
    dev.off()


    # testing with vst
    #vsd_t = vst(DDS, blind = FALSE)

    # for histone mark k36me2 vst fails so i have to use the direct function
    vsd_t <- tryCatch({
        vst(DDS, blind = FALSE)
    }, error = function(e) {
        message("vst() failed, using varianceStabilizingTransformation() instead.")
        varianceStabilizingTransformation(DDS, blind = FALSE)
    })

    head(assay(vsd_t), 5)

    pcaData2 = plotPCA(vsd_t, intgroup=c("condition_factor","type_factor"), returnData = TRUE )
    pcaData2

    percentVar <- round(100 * attr(pcaData2, "percentVar"))
    pca_plot_vst = ggplot(pcaData2, aes(PC1, PC2, color=condition_factor, shape=type_factor)) +
    ggtitle(paste("PCA plot using VST transform IDR", idr_used, sep = " ")) +
    theme(plot.caption = element_text(size = 30, face = "bold"), axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15), legend.text = element_text(size = 15))+
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()


    #name_for_vst_png_file =paste(experiment_design_name,"PCA_plot_IDR", idr_used, "vst.png", sep = "_")

    pdf(file = paste(experiment_design_name,"PCA_plot_IDR_our_data_counts",idr_used,"vst_all_reps.pdf", sep = "_"), width = 10, height = 10)
    print(pca_plot_vst)
    dev.off()

    png(file = paste(experiment_design_name,"PCA_plot_IDR_our_data_counts", idr_used,"vst_all_reps.png", sep = "_"), width = 900, height = 900)
    print(pca_plot_vst)
    dev.off()

    pca_plot_rlog

    pca_plot_vst


    rtracklayer::export.bed(row.names(up_reg), con = "${up_peaks_out}")

    rtracklayer::export.bed(row.names(down_reg), con = "${down_peaks_out}")

    # now exporting the unchanging peaks

    rtracklayer::export.bed(row.names(unchanging_reg), con = "${unchanging_peaks_out}")

    # now hoping to get the peak lengths histone


    ###### if any errors happen here then dont do anything ######
    tryCatch({
    up_peaks = read.table(file = './${up_peaks_out}', header = FALSE, sep = "\t")

    up_peak_lengths = up_peaks\$V3 - up_peaks\$V2
    print(max(up_peak_lengths))

    png(filename = "hist_up_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
    hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
    axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
    dev.off()

    # just to view it here, not needed in nextflow here.
    #hist.default(up_peak_lengths, xaxt = "n", breaks = 1e2)
    #axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))

    down_peaks = read.table(file = './${down_peaks_out}', header = FALSE, sep = "\t")

    down_peak_lengths = down_peaks\$V3 - down_peaks\$V2
    print(max(down_peak_lengths))

    png(filename = "hist_down_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
    hist.default(down_peak_lengths, xaxt = "n", breaks = 1e2)
    axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
    dev.off()


    unchanging_peaks = read.table(file = './${unchanging_peaks_out}', header = FALSE, sep = "\t")

    unchanging_peak_lengths = unchanging_peaks\$V3 - unchanging_peaks\$V2
    print(max(unchanging_peak_lengths))

    png(filename = "hist_unchanging_regulated_${idr_histone}_peak_lengths.png", height = 1000, width = 1000)
    hist.default(unchanging_peak_lengths, xaxt = "n", breaks = 1e2)
    axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
    dev.off()

    ########## now violin plots for peak lengths ################

    up_peak_df = DataFrame(peak_lengths_${idr_histone} = up_peak_lengths, category = "up_peaks_lengths")
    down_peak_df = DataFrame(peak_lengths_${idr_histone} = down_peak_lengths, category = "down_peaks_lengths")
    unchanging_peak_df = DataFrame(peak_lengths_${idr_histone} = unchanging_peak_lengths, category = "unchanging_peaks_lengths")


    df_all_peak_lengths = rbind(up_peak_df, down_peak_df, unchanging_peak_df)

    #df_all_peak_lengths

    df_all_peak_lengths_gg =  ggplot(df_all_peak_lengths, aes( category, peak_lengths_${idr_histone}))
    
    all_peak_lengths_violin = df_all_peak_lengths_gg+
        geom_violin()+
        scale_y_continuous(labels = scales::label_number())
    ggsave("${idr_histone}_up_down_unchanging_peak_length_violin_plot.png", plot = all_peak_lengths_violin)
    
    
    }, error = function(x) {
    
    message("some of the peak files had no lenght so plotting is pointless")
    })

    """
}


/*

# making the matrix 

    countsMatrix = matrix(data = NA, length(masterPeak), length(list_of_mBams) )


    # I want to get the chr start end of the peaks and put them in the matrix as column names so I know which peaks I am looking at from the deseq2 results in the differential peaks.
    seqnames(masterPeak)

    # i should put the unique ids from the master peak object in the matrix as row names.
    unique_masterPeak_ids = paste0(as.character(seqnames(masterPeak)), ":" ,start(masterPeak), "-", end(masterPeak))

    rownames(countsMatrix) = unique_masterPeak_ids

    #getting the list of bam base names to add to the matrix column names
    list_bam_basenames = list()


    # for deseq2 i need the condition design. so hlow and scrm
    # then i will find a way to tally how many of each are there so i can automate the rep count
    condition_design = list()

    type_design = list()

    length_of_bams = length(list_of_mBams)

    for (x in c(1:length_of_bams) ) {
    
        path_bam = list_of_mBams[x]
        print(path_bam)
        #bam_basename = basename(path_bam)
        bam_tokens = strsplit(path_bam, split = "_")[[1]]
        
        # labeling the important tokens so it is easier to keep track of
        condition = bam_tokens[1]
        histone = bam_tokens[2]
        replicate = bam_tokens[3]
        
        
        # for later parts I need the condition and to know how many times to repeat them
        condition_design = append(condition_design, paste(condition,histone,sep="_"))
        
        
        # also get the replicate too for type design
        type_design = append(type_design, replicate)
        #type_design = append(type_design, paste(histone,replicate, sep="_"))
        
        
        
        # using chromVAR getCounts function
        fragment_counts = getCounts(path_bam, masterPeak, paired= TRUE, by_rg = FALSE, format = "bam")
        
        
        # putting the fragment counts in the column labeled by the bam name
        countsMatrix[,x] = counts(fragment_counts)[,1]
        
        list_bam_basenames = append(list_bam_basenames, bam_basename)
        
    }

    colnames(countsMatrix) = list_bam_basenames 
    print(countsMatrix)



*/

process mk_bedgraph_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'
    label 'normal_big_resources'

    publishDir "./bedgraph_for_seacr/", mode: 'copy', pattern: "*"

    debug true 

    input:

    // the bam is the first index and the bai is the second index
    tuple val(file_basename), path(bam_index_tuple) 

    path(ref_genome_fai)


    output:

    path("${bedgraph_out_name}"), emit: bedgraph_for_seacr


    script:

    bam_files = bam_index_tuple[0]
    bai_files = bam_index_tuple[1]

    // out file names

    genome_size_name = "${ref_genome_fai}.genome"

    out_bed_name = "${file_basename}.bed"
    cleaned_out_bed_name = "${file_basename}.clean.bed"
    fragments_out_bed_name = "${file_basename}.fragments.bed"

    bedgraph_out_name = "${file_basename}.fragments.bedgraph"


    """
    #!/usr/bin/env bash

    # get the correct 2 fields for ref genome size from the fai file

    cut -f1,2 ${ref_genome_fai} > ${genome_size_name}

    # following the steps from SRACR github to get the bedgraph input

    bedtools bamtobed -bedpe -i ${bam_files} > ${out_bed_name}
    awk '\$1==\$4 && \$6-\$2 < 1000 {print \$0}' ${out_bed_name} > ${cleaned_out_bed_name}
    cut -f 1,2,6 ${cleaned_out_bed_name} | sort -k1,1 -k2,2n -k3,3n > ${fragments_out_bed_name}

    bedtools genomecov -bg -i ${fragments_out_bed_name} -g ${genome_size_name} > ${bedgraph_out_name}

    echo "\$(ll  ${bedgraph_out_name})"
    """
}


// before this i need to make a proces for making the bedgraph files
process seacr_peakcalls_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/seacr_1.3_rj'
    label 'normal_big_resources'

    publishDir "./seacr_peaks/", mode: 'copy', pattern: "*"



    input:

    path(bedgraph_files)



    output:

    path("*.stringent.bed"), emit: seacr_peaks



    script:

    file_basename = bedgraph_files.baseName

    // seacr will just add this ".stringent.bed" on to the end of the base name
    stringent_non_seacr_out_name = "${file_basename}.non"
    stringent_norm_seacr_out_name = "${file_basename}.norm"

    """
    #!/usr/bin/env bash

    # trying non vs norm
    SEACR_1.3.sh ${bedgraph_files} 0.1 non stringent ${stringent_non_seacr_out_name}

    # both resulted in the same peaks being called
    # SEACR_1.3.sh \${bedgraph_files} 0.1 norm stringent \${stringent_norm_seacr_out_name}


    """
}

// making the bed file for sicer2 

process mk_bed_for_sicer2_process {

    label 'normal_big_resources'

    // i am able to put bam files directly as input becasue i added bedtools to the conda environment
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'

    //publishDir "./sicer2_peaks/", mode: 'copy', pattern:'*'


    input:

    tuple val(tuple_key), path(bam_bai_list)


    output:

    path("*.bed"), emit: bed_for_sicer2



    script:

    bam_file = bam_bai_list[0]
    bai_file = bam_bai_list[1]

    bam_basename = bam_file.baseName

    bed_out_name = "${bam_basename}.bed"


    """
    #!/usr/bin/env bash

    bedtools bamtobed \
    -i ${bam_file} \
    > ${bed_out_name}
    





    """

}

process sicer2_peakcall_process {

    label 'normal_big_resources'

    // i am able to put bam files directly as input becasue i added bedtools to the conda environment
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/sicer2_rj'

    publishDir "./sicer2_peaks/", mode: 'copy', pattern:'*'


    input:

    //tuple val(tuple_key), path(bam_bai_file)

    path(bed_file)


    output:

    // this output is cgisland because i am using recognicer instead of sicer
    path("*.cgisland"), emit: sicer2_cgislands

    path("*normalized.wig"), emit: sicer2_wig

    path("*islandfiltered.bed"), emit: sicer2_peak_file


    script:

    //bam_file = bam_bai_file[0]
    //bai_file = bam_bai_file[1]

    """
    #!/usr/bin/env bash

    ###### sicer2 parameters ##########

    # only for sicer, i am using recognicer for broad peaks. --e_value : this requires user input if no control is set ( i am only using a treatment) default 1000
    
    # --step_size : the number of windows in one graining unit. default is 3
    # --step_score : the minimum number of positive elements in the graining unit to call the unit positive. Default value is 2.
    
    ###################################

    
    recognicer \
    -t ${bed_file} \
    --species hg38 \
    --step_size 3 \
    --step_score 2 \
    --significant_reads

    #sicer \
    -t \${bam_file} \
    --e_value 1000 \
    --species hg38 \
    --significant_reads








    """
}

/*
process seacr_idr_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/idr-2.0_rj'
    label 'super_big_resources'
    publishDir "./seacr_peaks/seacr_idr_results/${histone[0]}/${condition[0]}", mode: 'copy', pattern:'*'


    input:

    tuple val(grouping_key), val(condition), val(histone), val(replicate), val(bio_rep), val(file_name), val(basename), path(peakpath)
    // for the file_name there will be three replicates, all in order. this will be the case for each instance this process is called when parallelized
    // idr takes only 2 replicates at a time
    // we will do rep 1 and 2, then rep 2 and 3, then rep 1/2 and 2/3
    // this process will only work for if we have three replicates in our data


    output:
    // now i need to output the final idr broadpeak and then get all the others for this experiment histone mark and combine them
    // I can create another process that will filter these using the blacklist and bedtools intersect but lets look at it as is for now. 
    path("*broadPeak"), emit: idr_peaks

    path("*.png"), emit:idr_pngs



    script:

    // the conditions is always the same but it has three in the val above. ex: [Hlow, Hlow, Hlow]
    condition_label = condition[0] // so i can just use one of them to rebuild a file name if needed

    // same with the histone ex: [H3k27me3, H3k27me3, H3k27me3]
    histone_label = histone[0]

    // the replicates are in order but not the same ex: [r1, r2, r3]
    // but i can still store the correct thing
    rep_label1 = replicate[0]
    rep_label2 = replicate[1]
    rep_label3 = replicate[2]

    // same with bio_reps, they are in order and correspond properly with the order of rep_labels
    bio_label1 = bio_rep[0]
    bio_label2 = bio_rep[1]
    bio_label3 = bio_rep[2]

    // these will be the broadpeak file names, the paths are already in the process directory so i can use the correct name in order. the paths will not be in order
    peak1 = file_name[0]
    peak2 = file_name[1]
    peak3 = file_name[2]

    // making the sorted filename
    sort_peak1 = "sorted_${peak1}.gz"
    sort_peak2 = "sorted_${peak2}.gz"
    sort_peak3 = "sorted_${peak3}.gz"

    // now to create the output names for checking each replicate

    idr_out_1_2_name = "IDR_${condition_label}_${histone_label}_${rep_label1}_vs_${rep_label2}_${params.return_idr}.broadPeak"
    idr_out_2_3_name = "IDR_${condition_label}_${histone_label}_${rep_label2}_vs_${rep_label3}_${params.return_idr}.broadPeak"
    idr_out_1_3_name = "IDR_${condition_label}_${histone_label}_${rep_label1}_vs_${rep_label3}_${params.return_idr}.broadPeak"
    idr_out_final_name = "IDR_${condition_label}_${histone_label}_${rep_label1}_vs_${rep_label2}_vs_${rep_label3}_${params.return_idr}.broadPeak"
    idr_out_final_merged_name = "IDR_${condition_label}_${histone_label}_${rep_label1}_vs_${rep_label2}_vs_${rep_label3}_2_merged_${params.return_idr}.broadPeak"

    concat_peaks_name = "concat_IDR_${condition_label}_${histone_label}_${rep_label1}_vs_${rep_label2}_vs_${rep_label3}_${params.return_idr}_pairs.broadPeak"

    // name for pooled broad peak files 
    broadpeak_pool = "broadpeak_pool.broadPeak"
    sort_broadpeak_pool = "sort_${broadpeak_pool}.gz"

    // now creating the file name for the idr's of all the reps

    // idr_final_sorted_file = "${condition_label}_${histone_label}_${rep_label1}_vs_${rep_label2}_vs_${rep_label3}_IDR_${params.return_idr}.broadPeak.gz"

    """
    #!/usr/bin/env bash

    ###### idr parameters ########
    # plot_idr = 0.05 // the default is 0.05 to report in the png plots which peaks passed this threshold
    # return_idr = 1  // the default is all peaks will be returned even if the report plots the ones that pass a certain number. 1 as default here will give back all peaks like idr has set

    #--peak-list is not provided
    #Peaks are grouped by overlap and then merged. The merged peak aggregate value is determined by --peak-merge-method.

    #Peaks that don't overlap another peak in every other replicate are not included unless --use-nonoverlapping-peaks is set.
    ##############################

    # i need to pool the 3 reps so i can have the pooled file for input in idr
    cat ${peak1} ${peak2} ${peak3} > ${broadpeak_pool}

    # i have to gzip the file
    #gzip -nc \${peak1} > "\${peak1}.gz"
    #gzip -nc \${peak2} > "\${peak2}.gz"
    #gzip -nc \${peak3} > "\${peak3}.gz"


    # now sorting all the peaks by their pvalue column
    sort -k8,8nr "${peak1}" | gzip > ${sort_peak1}
    sort -k8,8nr "${peak2}" | gzip > ${sort_peak2}
    sort -k8,8nr "${peak3}"| gzip > ${sort_peak3}
    sort -k8,8nr ${broadpeak_pool} | gzip > "${sort_broadpeak_pool}"

    # now i have everything to put into idr
    # i think its best to return all peaks that passed the merging criteria for rep1 vs rep2 and rep2 vs rep3
    # then when i do rep1_2 vs rep2_3 I only return the peaks that have the 0.05 IDR threshold passed; these peaks will also be the ones that pass the merging threshold.
    # you can change this idr by changing the value in the --return_idr parameter when running nextflow

    # i will make a if then statement where i will choose to run the idr if the length of the peak file is larger than 21 lines
    peak1_length=\$(less ${sort_peak1} | wc -l )
    peak2_length=\$(less ${sort_peak2} | wc -l )
    peak3_length=\$(less ${sort_peak3} | wc -l )

    # rep 1 vs 2
    
    if ((\$peak1_length > 21 && \$peak2_length > 21)); then
        idr --samples ${sort_peak1} ${sort_peak2} \
        --input-file-type broadPeak \
        --output-file ${idr_out_1_2_name} \
        --rank p.value \
        --idr-threshold ${params.return_idr} \
        --soft-idr-threshold ${params.plot_idr} \
        --plot \
        --use-best-multisummit-IDR
    else
        touch ${idr_out_1_2_name}
    fi 

    # now rep 2 vs 3

    if ((\$peak2_length > 21 && \$peak3_length > 21)); then
        idr --samples ${sort_peak2} ${sort_peak3} \
        --input-file-type broadPeak \
        --output-file ${idr_out_2_3_name} \
        --rank p.value \
        --idr-threshold ${params.return_idr} \
        --soft-idr-threshold ${params.plot_idr} \
        --plot \
        --use-best-multisummit-IDR
    else
        touch ${idr_out_2_3_name}
    fi 

    # now doing 1 vs 3

    if ((\$peak1_length > 21 && \$peak3_length > 21)); then
        idr --samples ${sort_peak1} ${sort_peak3} \
        --input-file-type broadPeak \
        --output-file ${idr_out_1_3_name} \
        --rank p.value \
        --idr-threshold ${params.return_idr} \
        --soft-idr-threshold ${params.plot_idr} \
        --plot \
        --use-best-multisummit-IDR
    else
        touch ${idr_out_1_3_name}
    fi 

    # maybe we dont look at the final output and just select the output of peak pairs that has the most peaks that pass the 0.05 threshold
    # this is according to section 4d of the encode 3 pipeline 'https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit?tab=t.0#heading=h.9ecc41kilcvq' 
    # I will instead load all the files into R and use only the one from each condition that has the max number of peaks called
    # then I will merge the two files; the max from hlow and the max from scrm to get the masterpeak
    
    #idr_1_2_count=\$(less \${idr_out_1_2_name} | wc -l)
    #idr_1_3_count=\$(less \${idr_out_1_3_name} | wc -l)
    #idr_2_3_count=\$(less \${idr_out_2_3_name} | wc -l)


        


    # now the final output
    # this would be merging according to idr standards, keeping all the peaks from the 2 of the 3 pairs
    
    if ((\$peak1_length > 21 && \$peak2_length > 21 && \$peak3_length > 21)); then
        idr --samples ${idr_out_1_2_name} ${idr_out_2_3_name} \
        --input-file-type broadPeak \
        --output-file ${idr_out_final_merged_name} \
        --rank p.value \
        --plot \
        --use-best-multisummit-IDR
    fi


    # just concatenate them and see
    
    cat ${idr_out_1_2_name} ${idr_out_2_3_name} ${idr_out_1_3_name} > ${concat_peaks_name}
    



    #if ((\$peak1_length > 21 && \$peak2_length > 21 && \$peak3_length > 21)); then
        #idr --samples \${idr_out_1_2_name} \${idr_out_2_3_name} \
        --input-file-type broadPeak \
        --output-file \${idr_out_final_name} \
        --idr-threshold \${params.return_idr} \
        --soft-idr-threshold \${params.plot_idr} \
        --rank p.value \
        --plot \
        --use-best-multisummit-IDR
    #fi


    # will have to ask johanna what this next bit of code does, but i can convert it to work here in nextflow as i did above
    # it checks to see if column 12 is greater than or equal to the idr_thresh_transformed, if so print all the columns

    # dont need this since it might be using the wrong column and IDR already has a parameter to get only peaks passing the threshold
    
    # IDR_THRESH_TRANSFORMED=\$(awk -v p=0.05 'BEGIN{print -log(p)/log(10)}')

    # can change this to print 0 instead of listing all of them in. that will print all the columns
    
    # if you read the documentation, broad peak outputs do not have a summit columns (10,18,22 for 2 replicates), so that means i am actually using column 11 not 12 for broad peaks not narrow peaks
    #awk 'BEGIN{OFS="\t"} \$11>='"\${IDR_THRESH_TRANSFORMED}"' {print \$0}' \${idr_out_final_name} | \
    sort | uniq | sort -k7n,7n | gzip -nc > \${idr_final_sorted_file}



    """



}
*/


process multiqc_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/multiqc_rj'
    label 'normal_small_resources'
    publishDir "./dup_info", mode: 'copy', pattern: '*'

    input:
    tuple val(condition) , val(histone), path(log_files)


    output:

    path("*"), emit: multiqc_dup_info

    script:

    // log_basename = log_files[0].baseName
    // tokens = log_basename.tokenize("_")
    // condition = tokens[0]
    // histone = tokens[1]


    """
    #!/usr/bin/env bash

    multiqc . 
    
    mv multiqc_report.html multiqc_${histone}_report.html



    """


}

process plot_histones_at_peaks_process {

    //debug true

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools-3.5.6_rj'

    label 'normal_big_resources'

    //label 'super_big_resources'

    publishDir "./heatmaps", mode: 'copy', pattern: '*'

    input:
    tuple val(grouping_key), val(condition_label), val(histone_label), val(replicate_label), val(bw_peak_names), path(bw_peak_filepath)

    //tuple val(peak_condition_label), val(peak_histone_label), val(peak_replicate_label), val(peak_names), path(peak_filepath)



    output:

    path("${heatmap_out_name}"), emit: histone_peak_heatmap
    path("${profile_out_name}"), emit: profile_peak_heatmap




    script:

    bigwig_file_name = bw_peak_names[0]
    peak_file_name = bw_peak_names[1]



    out_matrix_name = "matrix_${grouping_key}.mat.gz"

    heatmap_out_name = "${grouping_key}_bigwig_signal_over_${grouping_key}_peaks_heatmap.png"


    profile_out_name = "${grouping_key}_bigwig_signal_over_${grouping_key}_peaks_profile.png"

    
    true_bw_name = "${bigwig_file_name}".replaceFirst(/\..*/, '')
    
    
    //name_list = []

    //num_names = bw_names.size()

    // for (int i=0; i< num_names; i++) {
        
        
    //     bw_true_name = "${bw_names[i]}".replaceFirst(/\..*/, '')

    //     name_list << bw_true_name
    // }

    """
    #!/usr/bin/env bash

    ###### deeptools parameters ####



    #################################

    # this will help me debug and make sure the bigwig and peak channels were aligned

    echo "the bigwig file for S flag is: ${bigwig_file_name}, and the peak file for R flag is: ${peak_file_name}"

    computeMatrix scale-regions \
    -S ${bigwig_file_name} \
    -R ${peak_file_name}\
    --outFileName "${out_matrix_name}" \
    --skipZeros \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --numberOfProcessors max


    echo "true bigwig name ${true_bw_name}"
    plotHeatmap \
    -m "${out_matrix_name}" \
    --outFileName "${heatmap_out_name}" \
    --sortUsing max \
    --heatmapWidth 15 \
    --heatmapHeight 20 \
    --dpi 300 \
    --labelRotation 30 \
    --samplesLabel ${true_bw_name} \
    --regionsLabel "${true_bw_name} called peaks"

    plotProfile \
    -m "${out_matrix_name}" \
    --outFileName "${profile_out_name}" \
    --plotWidth 8 \
    --plotHeight 10 \
    --labelRotation 30 \
    --dpi 300 \
    --samplesLabel ${true_bw_name}






    """
}





// this process will be very similar to the plot histone at genes process, but just at the peaks
// this will not be good becasue it uses both histone marks over the peaks that come from k27
// I have to do the version above instead
process plot_at_up_down_peaks_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools-3.5.6_rj'

    label 'normal_big_resources'
    
    //label 'super_big_resources'

    publishDir "./heatmaps/histone_signal_over_up_down_peaks", mode: 'copy', pattern: '*'


    input:

    debug true

    errorStrategy 'ignore'

    tuple  val(histone_label), val(condition_label),  val(replicate_label), val(bw_names), path(bigwig_filepath), val(peak_type), path(peak_filepath)

    // this is the version below that will use the geo control data, uncomment it and also the up peaks, down_peaks, and master_peaks
    // comment this below and up, down, master peaks if using the above that has all of it inside
    //tuple  val(histone_label), val(condition_label),  val(replicate_label), val(bw_names), path(bigwig_filepath)
    

    //path(up_peaks) //comment out if using data generated in pipeline for master peaks

    //path(down_peaks)  //comment out if using data generated in pipeline for master peaks

    path(bisulfate_cpg_bigwig)

    //path(master_peaks)   //comment out if using data generated in pipeline for master peaks

    path(cpg_island_unmasked_bed)



    output:

    // path("${png_heatmap_uppeaks}"), emit: uppeaks_histone_heatmap

    // path("${png_heatmap_down_peaks}"), emit: down_peaks_histone_heatmap

    path("${png_heatmap_both}"), emit: both_histone_heatmap
    //path("${png_wt_heatmap}"), emit: gene_histone_heatmap_wt

    path("*.png"), emit: all_png_files
    path("*.svg"), emit: all_svg_files



    script:

    // here i need to make sure I choose only the peaks that match the experiment type
    // these should be matched already but now to get all four peaks up, down, unchanging, and master peaks
    up_peaks = peak_filepath[0]
    down_peaks = peak_filepath[1]
    unchanging_peaks = peak_filepath[2]
    master_peaks = peak_filepath[3]

    // i need to get the bigwig name from the bisulfate bigwig file path


    bisulfate_bigwig_name = "${bisulfate_cpg_bigwig.name}"

    true_basename_bisulfate = "${bisulfate_bigwig_name}".replaceFirst(/\._fp*/, '')

    name_list = []

    num_files = bw_names.size()
    
    for (int i = 0; i < num_files; i++) {
        true_basename = "${bw_names[i]}".replaceFirst(/\..*/, '')
        name_list << true_basename

    }

    //newName = bw_names.replaceFirst(/\..*/, '')
    //list_control_bw_names = "${control_bw_names.toList()}"

    // out_matrix_scores_uppeaks = "matrix_peaks_${histone_label}_up_peaks.mat.gz"
    // //out_matrix_scores_wt = "matrix_gene_${histone_label}_lowup_genebody.mat.gz"

    // out_matrix_scores_down_peaks = "matrix_peaks_${histone_label}_down_peaks.mat.gz"

    out_matrix_scores_both = "matrix_peaks_${histone_label}_up_and_down_peaks.mat.gz"

    // png_heatmap_uppeaks = "${histone_label}_${replicate_label}_histone_features_at_up_peaks.png"
    // png_profile_up_peaks = "${histone_label}_${replicate_label}_histone_features_at_up_peaks_profile.png"
    // png_heatmap_down_peaks = "${histone_label}_${replicate_label}_histone_features_at_down_peaks.png"
    // png_profile_down_peaks = "${histone_label}_${replicate_label}_histone_features_at_down_peaks_profile.png"
    png_heatmap_both = "${histone_label}_${replicate_label}_histone_features_at_up_down_and_master_peaks_cpg_islands.png"
    png_profile_both_peaks = "${histone_label}_${replicate_label}_histone_features_at_up_down_and_master_peaks_cpg_islands_profile.png" 

    // exporting to svg to edit text size manually

    svg_heatmap_both = "${histone_label}_${replicate_label}_histone_features_at_up_down_and_master_peaks_cpg_islands.svg"
    svg_profile_both_peaks = "${histone_label}_${replicate_label}_histone_features_at_up_down_and_master_peaks_cpg_islands_profile.svg" 


    //png_wt_heatmap = "${wt_histone_label}_histone_features_at_lowup_genebody.png"

    // making output name for genes file that will have no zeros

    master_split_basename = "${master_peaks.baseName}".split('_')
    short_master_peak_name = "${master_split_basename[0]}_${master_split_basename[1]}"

    up_split_basename = "${up_peaks.baseName}".split('_')
    short_up_peak_name = "${up_split_basename[0]}_${up_split_basename[1]}_${up_split_basename[3]}"

    down_split_basename = "${down_peaks.baseName}".split('_')
    short_down_peak_name = "${down_split_basename[0]}_${down_split_basename[1]}_${down_split_basename[3]}"


    //up_peaks_nozero = "${up_peaks.baseName}_noZerolength.bed"
    up_peaks_nozero = "${short_up_peak_name}_noZerolength.bed"

    //down_peaks_nozero = "${down_peaks.baseName}_noZerolength.bed"
    down_peaks_nozero = "${short_down_peak_name}_noZerolength.bed"

    //master_peaks_nozero ="${master_peaks.baseName}_noZerolength.bed"
    master_peaks_nozero ="${short_master_peak_name}_noZerolength.bed"
    

    // i need to just get the matrix for plotting signal over cpg islands and then the output heatmap and profile plot names

    //out_matrix_cpg_islands = "matrix_peaks_${histone_label}_unmasked_cpg_regions.mat.gz"

    //cpg_png_heatmap

    """
    #!/usr/bin/env bash

    ########### deeptools params ##########



    #######################################

    #echo ' this is the list of names: "\${name_list}"'

    # fix the genebody file
    awk  '\$2!=\$3 {print \$0}' "${up_peaks}" > "${up_peaks_nozero}"

    #computeMatrix reference-point -S \${bw_names.join(' ')} \
    -R "\${up_peaks_nozero}" \
    --referencePoint center \
    --beforeRegionStartLength 50000 \
    --afterRegionStartLength 50000 \
    --skipZeros \
    --quiet \
    --binSize 200 \
    --numberOfProcessors "max" \
    -o "\${out_matrix_scores_uppeaks}"

    
    #plotHeatmap -m "\${out_matrix_scores_uppeaks}" \
    -out "\${png_heatmap_uppeaks}" \
    --colorMap RdBu_r \
    --samplesLabel \${name_list.join(' ')} \
    --labelRotation 30 \
    --heatmapWidth 8 \
    --sortUsing sum

    #plotProfile -m "\${out_matrix_scores_uppeaks}" \
    -out "\${png_profile_up_peaks}" \
    --plotHeight 10 \
    --plotWidth 10 \
    --perGroup \
    --samplesLabel \${name_list.join(' ')} \
    --plotTitle "signal over up peaks"


    # now to do this with the other genebody file 'down unchanging'

    awk  '\$2!=\$3 {print \$0}' "${down_peaks}" > "${down_peaks_nozero}"

    awk  '\$2!=\$3 {print \$0}' "${master_peaks}" > "${master_peaks_nozero}"

    #computeMatrix reference-point -S \${bw_names.join(' ')} \
    -R "\${down_peaks_nozero}" \
    --referencePoint center \
    --beforeRegionStartLength 20000 \
    --afterRegionStartLength 20000 \
    --skipZeros \
    --quiet \
    --binSize 200 \
    --numberOfProcessors "max" \
    -o "\${out_matrix_scores_down_peaks}"

    
    #plotHeatmap -m "\${out_matrix_scores_down_peaks}" \
    -out "\${png_heatmap_down_peaks}" \
    --colorMap RdBu_r \
    --samplesLabel \${name_list.join(' ')} \
    --labelRotation 30 \
    --heatmapWidth 8 \
    --sortUsing sum

    # now I want to use plot profile to make a plot per bed file(the genes) instead of per bigwig file (the signal), this will put all the signal bigwig in the same plot

    #plotProfile -m "\${out_matrix_scores_down_peaks}" \
    -out "\${png_profile_down_peaks}" \
    --plotHeight 10 \
    --plotWidth 10 \
    --perGroup \
    --samplesLabel \${name_list.join(' ')} \
    --plotTitle "signal over down peaks"


    echo ' the up genes and down unchanging genes plot has completed. Now starting on the plot for both up and down genes together'


    # now plotting both together

    # I THINK THE BEST THING IS TO ONLY HAVE BOTH UP AND DOWN TOGETHER

    computeMatrix reference-point -S ${bw_names.join(' ')} ${bisulfate_bigwig_name} \
    -R "${up_peaks_nozero}" "${down_peaks_nozero}" "${master_peaks_nozero}" "${cpg_island_unmasked_bed}" \
    --referencePoint center \
    --beforeRegionStartLength 50000 \
    --afterRegionStartLength 50000 \
    --skipZeros \
    --quiet \
    --binSize 1000 \
    --numberOfProcessors "max" \
    -o "${out_matrix_scores_both}"

    plotHeatmap -m "${out_matrix_scores_both}" \
    -out "${png_heatmap_both}" \
    --colorMap 'Reds' \
    --zMin 0 \
    --zMax 3 \
    --samplesLabel ${name_list.join(' ')} 'CpG_site_bigwig_signal' \
    --labelRotation 30 \
    --sortUsing sum \
    --perGroup \
    --heatmapWidth 8 \
    --heatmapHeight 20 \
    --dpi 300 
    #--plotTitle "Bigwig Signal and CpG signal Over Up and Down Peaks"

    plotProfile -m "${out_matrix_scores_both}" \
    -out "${png_profile_both_peaks}" \
    --plotHeight 20 \
    --plotWidth 20 \
    --perGroup \
    --dpi 300 \
    --samplesLabel ${name_list.join(' ')} 'CpG_site_bigwig_signal' \
    --plotTitle "Bigwig Signal and CpG signal over both peaks"

    // svg

    plotHeatmap -m "${out_matrix_scores_both}" \
    -out "${svg_heatmap_both}" \
    --colorMap 'Reds' \
    --zMin 0 \
    --zMax 3 \
    --samplesLabel ${name_list.join(' ')} 'CpG_site_bigwig_signal' \
    --labelRotation 30 \
    --sortUsing sum \
    --perGroup \
    --heatmapWidth 8 \
    --heatmapHeight 20 \
    --dpi 300 
    #--plotTitle "Bigwig Signal and CpG signal Over Up and Down Peaks"

    plotProfile -m "${out_matrix_scores_both}" \
    -out "${svg_profile_both_peaks}" \
    --plotHeight 20 \
    --plotWidth 20 \
    --perGroup \
    --dpi 300 \
    --samplesLabel ${name_list.join(' ')} 'CpG_site_bigwig_signal' \
    --plotTitle "Bigwig Signal and CpG signal over both peaks"





    """


}

// just putting this process close here to plot atac signal over cut&run peaks

process atac_signal_over_peaks_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools-3.5.6_rj'

    label 'normal_big_resources'
    
    //label 'super_big_resources'

    errorStrategy 'ignore'

    publishDir "./heatmaps/atac_signal_over_peaks/", mode: 'copy', pattern: '*'


    input:
    path(control_atac_bigwig)
    path(treatment_atac_bigwig) 
    
    // i need to remember to uncomment this if i am using the peaks from outside
    // path(up_peaks)
    // path(down_peaks)
    // path(unchanging_peaks)

    // new additions
    path(down_atac_peaks)
    path(up_atac_peaks)

    // have to change below to reflect the grouping with bigwig and peaks
    //tuple val(condition_label), val(histone_label), val(replicate_label), val(bw_names), path(bigwig_filepath)
    
    // new version with peaks
    tuple  val(histone_label), val(condition_label),  val(replicate_label), val(bw_names), path(bigwig_filepath), val(peak_type), path(peak_filepath)

    // adding the cpg island regions to see if they are being more accessible
    path(cpg_island_unmasked_regions)




    output:

    path("*.{png,svg}"), emit: atac_png_plots


    script:

    up_peaks = peak_filepath[0]
    down_peaks = peak_filepath[1]
    unchanging_peaks = peak_filepath[2]
    master_peaks = peak_filepath[3]

    out_matrix_scores_3 = "matrix_atac_bw_signal_over_${histone_label}_peaks.mat.gz"

    png_heatmap_3 = "atac_bigwig_signal_features_at_all_${histone_label}_peaks_cpg_regions_heatmap.png"
    svg_heatmap_3 = "atac_bigwig_signal_features_at_all_${histone_label}_peaks_cpg_regions_heatmap.svg"
    png_profile_3 = "atac_bigwig_signal_features_at_all_${histone_label}_peaks_cpg_regions_profile.png" 


    // fixing the names over the plots and adding it for unchanging also

    unchanging_split_basename = "${unchanging_peaks.baseName}".split('_')
    short_unchanging_peak_name = "${unchanging_split_basename[0]}_${unchanging_split_basename[1]}_${unchanging_split_basename[3]}"

    up_split_basename = "${up_peaks.baseName}".split('_')
    short_up_peak_name = "${up_split_basename[0]}_${up_split_basename[1]}_${up_split_basename[3]}"

    down_split_basename = "${down_peaks.baseName}".split('_')
    short_down_peak_name = "${down_split_basename[0]}_${down_split_basename[1]}_${down_split_basename[3]}"


    //up_peaks_nozero = "${up_peaks.baseName}_noZerolength.bed"
    up_peaks_nozero = "${short_up_peak_name}_noZerolength.bed"

    //down_peaks_nozero = "${down_peaks.baseName}_noZerolength.bed"
    down_peaks_nozero = "${short_down_peak_name}_noZerolength.bed"

    //master_peaks_nozero ="${master_peaks.baseName}_noZerolength.bed"
    unchanging_peaks_nozero ="${short_unchanging_peak_name}_noZerolength.bed"


    //up_peaks_nozero = "${up_peaks.baseName}_noZerolength.bed"
    //down_peaks_nozero = "${down_peaks.baseName}_noZerolength.bed"
    //unchanging_peaks_nozero = "${unchanging_peaks.baseName}_noZerolength.bed"

    out_matrix_scores_atac_peaks = "matrix_atac_bw_signal_over_atac_peaks.mat.gz"
    out_matrix_scores_expr_signal_atac_peaks = "matrix_expr_bw_signal_over_atac_peaks.mat.gz"

    svg_heatmap_atac_signal_over_atac_peaks = "atac_bigwig_signal_features_at_atac_peaks_heatmap.svg"
    png_heatmap_atac_signal_over_atac_peaks = "atac_bigwig_signal_features_at_atac_peaks_heatmap.png"
    
    svg_heatmap_expr_signal_over_atac_peaks = "${histone_label}_${replicate_label}_bigwig_signal_features_at_atac_peaks_heatmap.svg" 
    png_heatmap_expr_signal_over_atac_peaks = "${histone_label}_${replicate_label}_bigwig_signal_features_at_atac_peaks_heatmap.png" 

    
    name_list = []

    num_files = bw_names.size()
    
    for (int i = 0; i < num_files; i++) {
        true_basename = "${bw_names[i]}".replaceFirst(/\..*/, '')
        name_list << true_basename

    }



    """
    #!/usr/bin/env bash

    awk  '\$2!=\$3 {print \$0}' "${up_peaks}" > "${up_peaks_nozero}"
    awk  '\$2!=\$3 {print \$0}' "${down_peaks}" > "${down_peaks_nozero}"
    awk  '\$2!=\$3 {print \$0}' "${unchanging_peaks}" > "${unchanging_peaks_nozero}"


    computeMatrix reference-point -S ${control_atac_bigwig} ${treatment_atac_bigwig} \
    -R "${up_peaks_nozero}" "${down_peaks_nozero}" "${unchanging_peaks_nozero}" "${cpg_island_unmasked_regions}" \
    --referencePoint center \
    --beforeRegionStartLength 50000 \
    --afterRegionStartLength 50000 \
    --skipZeros \
    --quiet \
    --binSize 200 \
    --numberOfProcessors "max" \
    -o "${out_matrix_scores_3}"

    plotHeatmap -m "${out_matrix_scores_3}" \
    -out "${png_heatmap_3}" \
    --colorMap 'Reds' \
    --zMin 0 \
    --zMax "auto" \
    --samplesLabel ${control_atac_bigwig} ${treatment_atac_bigwig} \
    --labelRotation 30 \
    --sortUsing sum \
    --perGroup \
    --heatmapWidth 9 \
    --heatmapHeight 15 \
    --dpi 300 
    #--plotTitle "ATAC Bigwig Signal Over Up and Down and Unchanging Peaks"

    # making the svg heatmap first

    plotHeatmap -m "${out_matrix_scores_3}" \
    -out "${svg_heatmap_3}" \
    --colorMap 'Reds' \
    --zMin 0 \
    --zMax "auto" \
    --samplesLabel ${control_atac_bigwig} ${treatment_atac_bigwig} \
    --labelRotation 30 \
    --sortUsing sum \
    --perGroup \
    --heatmapWidth 9 \
    --heatmapHeight 15 \
    --dpi 300 
    #--plotTitle "ATAC Bigwig Signal Over Up and Down and Unchanging Peaks"

    plotProfile -m "${out_matrix_scores_3}" \
    -out "${png_profile_3}" \
    --plotHeight 10 \
    --plotWidth 10 \
    --perGroup \
    --dpi 300 \
    --samplesLabel ${control_atac_bigwig} ${treatment_atac_bigwig} \
    --plotTitle "ATAC signal over up down unchanging peaks"


    # now they want the atac peaks to be used and all signal plotted over them ##############

    computeMatrix reference-point -S ${control_atac_bigwig} ${treatment_atac_bigwig} \
    -R "${up_atac_peaks}" "${down_atac_peaks}"  \
    --referencePoint center \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --skipZeros \
    --quiet \
    --binSize 200 \
    --numberOfProcessors "max" \
    -o "${out_matrix_scores_atac_peaks}"

    plotHeatmap -m "${out_matrix_scores_atac_peaks}" \
    -out "${svg_heatmap_atac_signal_over_atac_peaks}" \
    --colorMap 'Reds' \
    --zMin 0 \
    --zMax "auto" \
    --regionsLabel "${up_atac_peaks.name}_ATAC_peaks" "${down_atac_peaks.name}_ATAC_peaks"  \
    --samplesLabel ${control_atac_bigwig} ${treatment_atac_bigwig} \
    --labelRotation 30 \
    --sortUsing sum \
    --perGroup \
    --heatmapWidth 9 \
    --heatmapHeight 15 \
    --dpi 300 
    #--plotTitle "ATAC Bigwig Signal Over atac Peaks"

    plotHeatmap -m "${out_matrix_scores_atac_peaks}" \
    -out "${png_heatmap_atac_signal_over_atac_peaks}" \
    --colorMap 'Reds' \
    --zMin 0 \
    --zMax "auto" \
    --regionsLabel "${up_atac_peaks.name}_ATAC_peaks" "${down_atac_peaks.name}_ATAC_peaks"  \
    --samplesLabel ${control_atac_bigwig} ${treatment_atac_bigwig} \
    --labelRotation 30 \
    --sortUsing sum \
    --perGroup \
    --heatmapWidth 9 \
    --heatmapHeight 15 \
    --dpi 300 
    #--plotTitle "ATAC Bigwig Signal Over atac Peaks"


    # doing experiment signal over atac peaks now #####################

    computeMatrix reference-point -S ${bw_names.join(' ')} \
    -R "${up_atac_peaks}" "${down_atac_peaks}"  \
    --referencePoint center \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --skipZeros \
    --quiet \
    --binSize 200 \
    --numberOfProcessors "max" \
    -o "${out_matrix_scores_expr_signal_atac_peaks}"

    plotHeatmap -m "${out_matrix_scores_expr_signal_atac_peaks}" \
    -out "${svg_heatmap_expr_signal_over_atac_peaks}" \
    --colorMap 'Reds' \
    --zMin 0 \
    --zMax "auto" \
    --regionsLabel "${up_atac_peaks.name}_ATAC_peaks" "${down_atac_peaks.name}_ATAC_peaks"  \
    --samplesLabel ${name_list.join(' ')} \
    --labelRotation 30 \
    --sortUsing sum \
    --perGroup \
    --heatmapWidth 9 \
    --heatmapHeight 15 \
    --dpi 300 
    #--plotTitle "${histone_label} Bigwig Signal Over atac Peaks"

    plotHeatmap -m "${out_matrix_scores_expr_signal_atac_peaks}" \
    -out "${png_heatmap_expr_signal_over_atac_peaks}" \
    --colorMap 'Reds' \
    --zMin 0 \
    --zMax "auto" \
    --regionsLabel "${up_atac_peaks.name}_ATAC_peaks" "${down_atac_peaks.name}_ATAC_peaks"  \
    --samplesLabel ${name_list.join(' ')} \
    --labelRotation 30 \
    --sortUsing sum \
    --perGroup \
    --heatmapWidth 9 \
    --heatmapHeight 15 \
    --dpi 300 
    #--plotTitle "${histone_label} Bigwig Signal Over atac Peaks"

    


    """
}

process bedtools_stranded_create_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'
    label 'normal_big_resources'
    publishDir "./diff_genes/", mode: 'copy', pattern: '*'


    input:
    path(up_genes)

    path(down_genes)

    path(nochange_genes)

    path(gtf_file)

    path(genome_size)



    output:

    path("${up_genes_20kb_slop}"), emit: up_genes_20kb_stranded
    path("${down_genes_20kb_slop}"), emit: down_genes_20kb_stranded
    path("${nochange_genes_20kb_slop}"), emit: nochange_genes_20kb_stranded




    script:

    stranded_up_genes = "${up_genes.baseName}_stranded.bed"
    stranded_down_genes = "${down_genes.baseName}_stranded.bed"
    stranded_nochange_genes = "${nochange_genes.baseName}_stranded.bed"

    // changing from 5kb to 20kb
    up_genes_20kb_slop = "${up_genes.baseName}_stranded_20kbSlop.bed"
    down_genes_20kb_slop = "${down_genes.baseName}_stranded_20kbSlop.bed"
    nochange_genes_20kb_slop = "${nochange_genes.baseName}_stranded_20kbSlop.bed"

    """
    #!/usr/bin/env bash

    # getting the strandedness onto the up genes

    awk 'FNR==NR{a[\$4];next} \$3=="gene" {match(\$0,/gene_id \\"([^\\"]+)\\"/,m); gsub(/\\..*/,\"\",m[1]); if(m[1] in a) print \$1 "\\t" \$4 "\\t" \$5 "\\t" m[1] "\\t." "\\t" \$7}' ${up_genes} ${gtf_file} > ${stranded_up_genes}


    # getting the strandedness onto the down genes
    awk 'FNR==NR{a[\$4];next} \$3=="gene" {match(\$0,/gene_id \\"([^\\"]+)\\"/,m); gsub(/\\..*/,\"\",m[1]); if(m[1] in a) print \$1 "\\t" \$4 "\\t" \$5 "\\t" m[1] "\\t." "\\t" \$7}'  ${down_genes} ${gtf_file} > ${stranded_down_genes}

    # getting strandedness on the nochange genes
    awk 'FNR==NR{a[\$4];next} \$3=="gene" {match(\$0,/gene_id \\"([^\\"]+)\\"/,m); gsub(/\\..*/,\"\",m[1]); if(m[1] in a) print \$1 "\\t" \$4 "\\t" \$5 "\\t" m[1] "\\t." "\\t" \$7}'  ${nochange_genes} ${gtf_file} > ${stranded_nochange_genes}



    # getting +5kb onto the genes based on strandedness
    # using the -s for strandedness and -l to add the 5000pb on the start coordinate but on the nevative strand it will add it to the end coordinate
    # have to use -b with -s or both -l and -r with -s 

    
    # changing from 5kb to 20kb

    bedtools slop \
    -i ${stranded_up_genes} \
    -g ${genome_size} \
    -s \
    -b 20000 \
    > ${up_genes_20kb_slop}


    # now doing the same for down genes 

    bedtools slop \
    -i ${stranded_down_genes} \
    -g ${genome_size} \
    -s \
    -b 20000 \
    > ${down_genes_20kb_slop}


    # now doing the same for nochange genes 

    bedtools slop \
    -i ${stranded_nochange_genes} \
    -g ${genome_size} \
    -s \
    -b 20000 \
    > ${nochange_genes_20kb_slop}

    """
}


process signal_over_gene_tss_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools-3.5.6_rj'

    label 'normal_big_resources'

    publishDir "heatmaps/signal_over_tss", mode: 'copy', pattern: '*'



    input:
    //path(up_peaks)
    //path(down_peaks)

    path(up_proseq_genes) //path(up_genes_20kb)

    path(down_proseq_genes)//path(down_genes_20kb)

    path(unchanging_proseq_genes) //path(nochange_genes_20kb)

    path(up_genes_norm)
    path(down_genes_norm)

    //tuple val(condition_label), val(histone_label), val(replicate_label), val(bw_names), path(bigwig_filepath)  // this is the version that uses data from outside the pipeline. uncomment if needed to change

    tuple val(histone_label), val(condition_label), val(replicate_label), val(bw_names), path(bigwig_filepath), val(peak_type), path(peak_filepath)


    output:

    path("*.{png,svg,pdf}"), emit: signal_over_gene_tss_heatmaps



    script:

    name_list = []

    num_files = bw_names.size()
    
    for (int i = 0; i < num_files; i++) {
        true_basename = "${bw_names[i]}".replaceFirst(/\..*/, '')
        name_list << true_basename

    }

    //concat_diff_peaks = "diff_peaks.bed"

    //up_peaks_nozero = "${up_peaks.baseName}_noZerolength.bed"
    //down_peaks_nozero = "${down_peaks.baseName}_noZerolength.bed"

    up_genes_tss_out_matrix_scores = "matrix_${histone_label}_${replicate_label}_signal_up_gene_tss.mat.gz"
    down_genes_tss_out_matrix_scores = "matrix_${histone_label}_${replicate_label}_signal_down_gene_tss.mat.gz"

    both_genes_tss_out_matrix_scores = "matrix_${histone_label}_${replicate_label}_signal_both_gene_tss.mat.gz"
    both_genes_tss_out_heatmap = "${histone_label}_${replicate_label}_signal_at_both_gene_tss_20kb.png"
    both_genes_tss_out_heatmap_svg = "${histone_label}_${replicate_label}_signal_at_both_gene_tss_20kb.svg"

    heatmap_only_pdf = "${histone_label}_${replicate_label}_signal_at_both_gene_tss_20kb_heatmap_only.pdf"


    """
    #!/usr/bin/env bash

    #cat \${up_peaks} \${down_peaks} > \${concat_diff_peaks}

    # fix the genebody file
    #awk  '\$2!=\$3 {print \$0}' "\${up_peaks}" > "\${up_peaks_nozero}"

    # now to do this with the other genebody file 'down unchanging'

    #awk  '\$2!=\$3 {print \$0}' "\${down_peaks}" > "\${down_peaks_nozero}"

    #cat \${up_peaks_nozero} \${down_peaks_nozero} > \${concat_diff_peaks}

    # now plotting 
    # because I put the TSS out 20kb, I will make the region around that new tss 20kb so I dont get anything in the gene body

    computeMatrix reference-point \
    -S ${bw_names.join(" ")} \
    -R ${up_proseq_genes} ${down_proseq_genes} ${unchanging_proseq_genes} \
    --referencePoint "TSS" \
    --beforeRegionStartLength 20000 \
    --afterRegionStartLength 20000 \
    --binSize 200 \
    --numberOfProcessors "max" \
    --outFileName ${both_genes_tss_out_matrix_scores}
    

    plotHeatmap \
    --matrixFile ${both_genes_tss_out_matrix_scores} \
    --outFileName ${both_genes_tss_out_heatmap} \
    --colorMap 'Reds' \
    --zMin 0 \
    --zMax "auto" \
    --samplesLabel ${name_list.join(' ')} \
    --labelRotation 30 \
    --heatmapWidth 9 \
    --heatmapHeight 15 \
    --dpi 300 \
    --sortUsing sum \
    --perGroup 

    # making the svg version
    plotHeatmap \
    --matrixFile ${both_genes_tss_out_matrix_scores} \
    --outFileName ${both_genes_tss_out_heatmap_svg} \
    --colorMap 'Reds' \
    --zMin 0 \
    --zMax "auto" \
    --samplesLabel ${name_list.join(' ')} \
    --labelRotation 30 \
    --heatmapWidth 9 \
    --heatmapHeight 15 \
    --dpi 300 \
    --sortUsing sum \
    --perGroup 

    # plotting the heatmap only
    plotHeatmap \
    --matrixFile ${both_genes_tss_out_matrix_scores} \
    --outFileName ${heatmap_only_pdf} \
    --colorMap 'Reds' \
    --whatToShow 'heatmap and colorbar' \
    --zMin 0 \
    --zMax "auto" \
    --samplesLabel ${name_list.join(' ')} \
    --labelRotation 30 \
    --heatmapWidth 9 \
    --heatmapHeight 15 \
    --dpi 300 \
    --sortUsing sum \
    --perGroup 








    """
}

process diff_peaks_intersect_diff_genes_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'
    label 'normal_big_resources'
    publishDir "./intersect_peaks_in_genes", mode: 'copy', pattern: '*'



    input:

    path(up_peaks)
    path(down_peaks)
    path(unchanging_peaks)
    path(master_peaks)
    path(up_20kb_genes)
    path(down_20kb_genes)
    path(nochange_20kb_genes)
    path(knownGene_bed)
    path(proseq_up_gene_ch)
    path(proseq_down_gene_ch)
    path(proseq_unchanging_gene_ch)



    output:

    tuple path("${out_up_genes_diff_peaks}"), path("${out_down_genes_diff_peaks}"), emit: diff_peaks_in_diff_genes_ch

    //tuple path("${out_closest_gene_to_up_peaks}"), path("${out_closest_gene_to_down_peaks}"), emit: closest_genes_to_peaks_ch

    path("*.bed"), emit: bed_files_ch

    script:

    

    sorted_knownGene_file = "${knownGene_bed.baseName}.sorted.bed"

    sorted_up_genes = "${up_20kb_genes.baseName}.sorted.bed"

    sorted_down_genes = "${down_20kb_genes.baseName}.sorted.bed"

    sorted_nochange_genes = "${nochange_20kb_genes.baseName}.sorted.bed"

    out_up_genes_diff_peaks = "diff_peaks_in_up_genes_20kb.bed"

    out_down_genes_diff_peaks = "diff_peaks_in_down_genes_20kb.bed"

    out_nochange_genes_diff_peaks = "diff_peaks_in_nochange_genes_20kb.bed"

    out_up_genes_master_peaks = "master_peaks_in_up_genes_20kb.bed"

    out_down_genes_master_peaks = "master_peaks_in_down_genes_20kb.bed"

    out_nochange_genes_master_peaks = "master_peaks_in_nochange_genes_20kb.bed"

    //out_closest_gene_to_up_peaks = "closest_gene_in_up_peaks.bed"

    //out_closest_gene_to_down_peaks = "closest_gene_in_down_peaks.bed"

    closest_up_proseq_gene_to_up_peaks = "closest_up_proseq_gene_up_peaks.bed"
    closest_up_proseq_gene_to_down_peaks = "closest_up_proseq_gene_down_peaks.bed"
    closest_up_proseq_gene_to_unchanging_peaks = "closest_up_proseq_gene_unchanging_peaks.bed"
    //closest_proseq_gene_to_unchanging_peaks = "closest_proseq_gene_unchanging_peaks.bed"

    // now i also want to use the normal up genes and find which of those are close to up or down peaks
    //closest_norm_up_genes_to_up_peaks = "closest_norm_up_genes_to_up_peaks.bed"
    //closest_norm_up_genes_to_down_peaks  = "closest_norm_up_genes_to_down_peaks.bed"
    //closest_norm_up_genes_to_unchanging_peaks  = "closest_norm_up_genes_to_unchanging_peaks.bed"


    // have to do this for the diff proseq genes
    closest_diff_genes_to_up_peaks = "closest_diff_genes_to_up_peaks.bed"
    closest_diff_genes_to_down_peaks = "closest_diff_genes_to_down_peaks.bed"
    closest_diff_genes_to_nochange_peaks = "closest_diff_genes_to_nochange_peaks.bed"

    closest_proseq_diff_genes_to_up_peaks = "closest_proseq_diff_genes_to_up_peaks.bed"
    closest_proseq_diff_genes_to_down_peaks = "closest_proseq_diff_genes_to_down_peaks.bed"
    closest_proseq_diff_genes_to_nochange_peaks = "closest_proseq_diff_genes_to_nochange_peaks.bed"


    """
    #!/usr/bin/env bash

    # sort the knownGene

    sort -k1,1 -k2,2n ${knownGene_bed} > ${sorted_knownGene_file}

    # sorting the files 

    sort -k1,1 -k2,2n ${up_20kb_genes} > ${sorted_up_genes}

    sort -k1,1 -k2,2n ${down_20kb_genes} > ${sorted_down_genes}

    sort -k1,1 -k2,2n ${nochange_20kb_genes} > ${sorted_nochange_genes}


    # looking at diff peaks in up genes first

    bedtools intersect \
    -a ${sorted_up_genes} \
    -b ${up_peaks} ${down_peaks} \
    -wa \
    -wb \
    -sorted \
    -filenames \
    > ${out_up_genes_diff_peaks}

    # now to get the diff peaks in down genes 
    bedtools intersect \
    -a ${sorted_down_genes} \
    -b ${up_peaks} ${down_peaks} \
    -wa \
    -wb \
    -sorted \
    -filenames \
    > ${out_down_genes_diff_peaks}

    # now to get the diff peaks in nochange genes 
    bedtools intersect \
    -a ${sorted_nochange_genes} \
    -b ${up_peaks} ${down_peaks} \
    -wa \
    -wb \
    -sorted \
    -filenames \
    > ${out_nochange_genes_diff_peaks}

    # now two separate files for looking at which master peaks are in up genes or down genes

    bedtools intersect \
    -a ${sorted_up_genes} \
    -b ${master_peaks}\
    -wa \
    -wb \
    -sorted \
    -filenames \
    > ${out_up_genes_master_peaks}

    bedtools intersect \
    -a ${sorted_down_genes} \
    -b ${master_peaks}\
    -wa \
    -wb \
    -sorted \
    -filenames \
    > ${out_down_genes_master_peaks}

    bedtools intersect \
    -a ${sorted_nochange_genes} \
    -b ${master_peaks}\
    -wa \
    -wb \
    -sorted \
    -filenames \
    > ${out_nochange_genes_master_peaks}

    # should just do the bedtools closest here

    # looking for genes that are closest to up peaks, then the genes that are closest to down peaks

    ############ params for bedtools closest ################

    # -D parameter: report the distance away from a that the closest b is found also, and tells if upstream or downstream (negative being upstream when using ref parameter).
    # -filenames : use the file names when reporting which gene came from where
    # -mdb : find the closest interval among all the files, instead of reporting the closest interval from all the files
    ########################################################

    # PROBABLY NOT ANALYZING THE CLOSEST PEAK IN A FULL KNOWNGENE FILE ANYMORE
    # first up peaks
    #bedtools closest \
    -a \${up_peaks} \
    -b \${sorted_knownGene_file} \
    -D ref \
    > \${out_closest_gene_to_up_peaks}



    # now down peaks
    #bedtools closest \
    -a \${down_peaks} \
    -b \${sorted_knownGene_file} \
    -D ref \
    > \${out_closest_gene_to_down_peaks}


    # using proseq genes
    # changing to finding the closest diff gene to the peaks among all diff gene files 

    # up proseq genes to up peaks
    #bedtools closest \
    -a \${up_peaks} \
    -b \${proseq_up_gene_ch} \
    -D ref \
    > \${closest_up_proseq_gene_to_up_peaks}

    # up proseq genes to down peaks
    #bedtools closest \
    -a \${down_peaks} \
    -b \${proseq_up_gene_ch} \
    -D ref \
    > \${closest_up_proseq_gene_to_down_peaks}

    # up proseq genes to unchanging peaks
    #bedtools closest \
    -a \${unchanging_peaks} \
    -b \${proseq_up_gene_ch} \
    -D ref \
    > \${closest_up_proseq_gene_to_unchanging_peaks}

    bedtools closest \
    -a ${up_peaks} \
    -b ${proseq_up_gene_ch} ${proseq_down_gene_ch} ${proseq_unchanging_gene_ch} \
    -D ref \
    -filenames \
    -mdb all \
    > ${closest_proseq_diff_genes_to_up_peaks}

    bedtools closest \
    -a ${down_peaks} \
    -b ${proseq_up_gene_ch} ${proseq_down_gene_ch} ${proseq_unchanging_gene_ch} \
    -D ref \
    -filenames \
    -mdb all \
    > ${closest_proseq_diff_genes_to_down_peaks}

    bedtools closest \
    -a ${unchanging_peaks} \
    -b ${proseq_up_gene_ch} ${proseq_down_gene_ch} ${proseq_unchanging_gene_ch} \
    -D ref \
    -filenames \
    -mdb all \
    > ${closest_proseq_diff_genes_to_nochange_peaks}



    # now dowing the normal up genes that are close to up and down peaks

    # changing this to have all of the genes (up, down, unchanging)
    # we want to see in each peak, what genes are closest to them and are the genes from the up, down or unchanging genes list

    #bedtools closest \
    -a \${up_peaks} \
    -b \${sorted_up_genes} \
    -D ref \
    > \${closest_norm_up_genes_to_up_peaks}


    #bedtools closest \
    -a \${down_peaks} \
    -b \${sorted_up_genes} \
    -D ref \
    > \${closest_norm_up_genes_to_down_peaks}

    #bedtools closest \
    -a \${unchanging_peaks} \
    -b \${sorted_up_genes} \
    -D ref \
    > \${closest_norm_up_genes_to_unchanging_peaks}

    # new version of finding the closest genes to each type of peak

    bedtools closest \
    -a ${up_peaks} \
    -b ${sorted_up_genes} ${sorted_down_genes} ${sorted_nochange_genes} \
    -D ref \
    -filenames \
    -mdb all \
    > ${closest_diff_genes_to_up_peaks}

    bedtools closest \
    -a ${down_peaks} \
    -b ${sorted_up_genes} ${sorted_down_genes} ${sorted_nochange_genes} \
    -D ref \
    -filenames \
    -mdb all \
    > ${closest_diff_genes_to_down_peaks}

    bedtools closest \
    -a ${unchanging_peaks} \
    -b ${sorted_up_genes} ${sorted_down_genes} ${sorted_nochange_genes} \
    -D ref \
    -filenames \
    -mdb all \
    > ${closest_diff_genes_to_nochange_peaks}




    """
}


process get_CpG_islands_in_peaks_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'
    label 'normal_big_resources'

    publishDir "./intersect_CpG_in_peaks", mode: 'copy', pattern: '*'


    input:

    // uncomment when using the outside peaks
    // path(up_peaks_ch)
    
    // path(down_peaks_ch)
    
    // path(unchanging_peaks_ch)
    
    path(cpg_island_unmasked_ch)


    // comment when using outside peaks
    // all true peaks

    tuple val(peak_type), val(exper_type), path(all_true_peaks_ch)
    //tuple val(histone_label), val(condition_label), val(replicate_label), val(bw_names), path(bigwig_filepath), val(peak_type), path(peak_filepath)



    output:

    path("${cpgIslands_in_up_peaks}"), emit: cpg_up_regions
    path("${cpgIslands_in_down_peaks}"), emit: cpg_down_regions
    path("${cpgIslands_in_unchanging_peaks}"), emit: cpg_unchanging_regions
    path("${cpgIslands_in_merged_master_peaks}"), emit: cpg_masterpeak_regions




    script:



    up_peaks_ch = "${all_true_peaks_ch[0]}"
    down_peaks_ch = "${all_true_peaks_ch[1]}"
    unchanging_peaks_ch = "${all_true_peaks_ch[2]}"
    masterpeaks_merged_peaks_ch = "${all_true_peaks_ch[3]}"

    // when i get a peaks channel that has multiple exper types(histones) then put the expr type in the name
    cpgIslands_in_up_peaks = "up_${exper_type}_cpg_islands.bed"
    cpgIslands_in_down_peaks = "down_${exper_type}_cpg_islands.bed"
    cpgIslands_in_unchanging_peaks = "unchanging_${exper_type}_cpg_islands.bed"
    
    // new addition
    cpgIslands_in_merged_master_peaks = "masterpeaks_${exper_type}_cpg_islands.bed"

    

    """
    #!/usr/bin/env bash

    # doing this 3 times 

    # first cpg islands in up peaks
    bedtools intersect \
    -a ${cpg_island_unmasked_ch} \
    -b ${up_peaks_ch} \
    -wa \
    > ${cpgIslands_in_up_peaks}

    # second down peaks
    bedtools intersect \
    -a ${cpg_island_unmasked_ch} \
    -b ${down_peaks_ch} \
    -wa \
    > ${cpgIslands_in_down_peaks}

    # third unchanging peaks
    bedtools intersect \
    -a ${cpg_island_unmasked_ch} \
    -b ${unchanging_peaks_ch} \
    -wa \
    > ${cpgIslands_in_unchanging_peaks}

    # new addition
    # fourth master peaks
    bedtools intersect \
    -a ${cpg_island_unmasked_ch} \
    -b ${masterpeaks_merged_peaks_ch} \
    -wa \
    > ${cpgIslands_in_merged_master_peaks}




    """
}

process plot_over_diff_cpg_regions_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools-3.5.6_rj'

    label 'normal_big_resources'

    publishDir "heatmaps/signal_over_diff_CpG_regions", mode: 'copy', pattern: '*'

    errorStrategy 'ignore'



    input:

    // path(cpg_up_ch)
    
    // path(cpg_down_ch)
    
    // path(cpg_unchanging_ch)

    // // new addition
    // path(cpg_master_ch)
    
    // have to change this to get the correct cpgs with the right bigwigs
    //combined_bigwig_meta_2grouped_ch
    //tuple val(condition_label), val(histone_label), val(replicate_label), val(bw_names), path(bigwig_filepath)

    tuple val(histone_label), val(condition_label), val(replicate_label), val(bw_names), path(bigwig_filepath), path(cpg_peak_filepath)



    output:

    path("*.{png,svg}"), emit: cpg_peak_heatmaps

    script:

    out_matrix_scores_expr_signal_cpg_in_peaks = "matrix_${histone_label}_signal_over_cpg_in_peaks.mat.gz"

    svg_heatmap_expr_signal_over_cpg_in_peaks = "${histone_label}_${replicate_label}_bigwig_signal_features_at_cpg_peaks_heatmap.svg"
    png_heatmap_expr_signal_over_cpg_in_peaks = "${histone_label}_${replicate_label}_bigwig_signal_features_at_cpg_peaks_heatmap.png"


    name_list = []

    num_files = bw_names.size()
    
    for (int i = 0; i < num_files; i++) {
        true_basename = "${bw_names[i]}".replaceFirst(/\..*/, '')
        name_list << true_basename

    }

    cpg_up_ch = "${cpg_peak_filepath[0]}"
    cpg_down_ch = "${cpg_peak_filepath[1]}"
    cpg_unchanging_ch = "${cpg_peak_filepath[2]}"
    cpg_master_ch = "${cpg_peak_filepath[3]}"

    """
    #!/usr/bin/env bash


    computeMatrix reference-point -S ${bw_names.join(' ')} \
    -R "${cpg_up_ch}" "${cpg_down_ch}" "${cpg_unchanging_ch}" "${cpg_master_ch}"  \
    --referencePoint center \
    --beforeRegionStartLength 10000 \
    --afterRegionStartLength 10000 \
    --skipZeros \
    --quiet \
    --binSize 200 \
    --numberOfProcessors "max" \
    -o "${out_matrix_scores_expr_signal_cpg_in_peaks}"

    plotHeatmap -m "${out_matrix_scores_expr_signal_cpg_in_peaks}" \
    -out "${svg_heatmap_expr_signal_over_cpg_in_peaks}" \
    --colorMap 'Reds' \
    --zMin 0 \
    --zMax "auto" \
    --regionsLabel "${cpg_up_ch}" "${cpg_down_ch}" "${cpg_unchanging_ch}" "${cpg_master_ch}" \
    --samplesLabel ${name_list.join(' ')} \
    --labelRotation 30 \
    --dpi 300 \
    --sortUsing sum \
    --perGroup \
    --heatmapWidth 9 \
    --heatmapHeight 15 
   # --plotTitle "${histone_label} Bigwig Signal Over CpG Islands in UP,DOWN,UNCHANGING Peaks"

    plotHeatmap -m "${out_matrix_scores_expr_signal_cpg_in_peaks}" \
    -out "${png_heatmap_expr_signal_over_cpg_in_peaks}" \
    --colorMap 'Reds' \
    --zMin 0 \
    --zMax "auto" \
    --regionsLabel "${cpg_up_ch}" "${cpg_down_ch}" "${cpg_unchanging_ch}" "${cpg_master_ch}" \
    --samplesLabel ${name_list.join(' ')} \
    --labelRotation 30 \
    --dpi 300 \
    --sortUsing sum \
    --perGroup \
    --heatmapWidth 9 \
    --heatmapHeight 15 
    #--plotTitle "${histone_label} Bigwig Signal Over CpG Islands in UP,DOWN,UNCHANGING Peaks"




    """
}

process atac_enrich_counts_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'
    label 'normal_big_resources'

    publishDir "./enrichment_of_experiment/", mode: 'copy', pattern: '*'

    debug true

    input:

    // these are the broad histones
    tuple val(roadmap_histones_filename), val(roadmap_histones_names), path(roadmap_histones_path)

    // these are the narrow histones
    tuple val(roadmap_narrowhistones_filename), val(roadmap_narrowhistones_names), path(roadmap_narrowhistones_path)

    // these are our own peak files to find which peaks are enriched in regions that go up in accessibility
    tuple val(pipeline_peaks_filename), val(pipeline_peaks_names), path(pipeline_peaks_path)

    //path(atac_bigwig)
    path(bam_files)


    output:

    path("*.png"), emit: atac_enrichment_plot
    path("*.tab"), emit: raw_enrichment_counts



    script:

    // getting the true base name of the pipeline peak files
    name_list = []
    file_list = []

    num_files = pipeline_peaks_names.size()
    
    for (int i = 0; i < num_files; i++) {

        

        //println("${pipeline_peaks_filename[i]}".readLines().size())
        //true_basename = "${pipeline_peaks_names[i]}".replaceFirst(/r{1,2,3}.*/, '')
        basename_split = "${pipeline_peaks_filename[i]}".split('_')


        condition = basename_split[2]
        exper_type = basename_split[3]

        true_basename = "${basename_split[0]}_${basename_split[1]}_${basename_split[2]}_${basename_split[3]}_"

        name_list << true_basename

        file_list << pipeline_peaks_filename[i]

        // try this
        //matcher = "${pipeline_peaks_names[i]}" =~ / (\d+)\_(\d+)\./
        
    }

    //pipeline_bed_files = name_list.collect{ file("${it}*new_sorted.bed") }.join(' ')
    //pipeline_regions_name = name_list.join(' ')

    //out_npz_file = "${atac_bigwig.baseName}_counts.npz"
    //out_tab_file = "${atac_bigwig.baseName}_counts.tab"

    // making the plot file name that will be outputted
    out_broad_enrich_plot = "${bam_files[1].baseName}_broad_enrichment_plot.png"
    out_broad_enrich_counts = "${bam_files[1].baseName}_broad_enrichment_counts.tab"

    out_narrow_enrich_plot = "${bam_files[1].baseName}_narrow_enrichment_plot.png"
    out_narrow_enrich_counts = "${bam_files[1].baseName}_narrow_enrichment_counts.tab"

    out_pipeline_enrich_plot = "${bam_files[1].baseName}_pipeline_enrichment_plot.png"
    out_pipeline_enrich_counts = "${bam_files[1].baseName}_pipeline_enrichment_counts.tab"

    """
    #!/usr/bin/env bash

    # i might have to clean each bed file so the start and stop are not the same

    file_list=(${roadmap_histones_path.join(' ')})

    for file in \${file_list[@]}; do

        file_basename=\$(basename \${file} .broadPeak)
        awk 'BEGIN{OFS="\t"} \$3 > \$2 {print \$1, \$2, \$3 }' \${file} > \${file_basename}_new.broadPeak

        sort -k1,1 -k2,2n \${file_basename}_new.broadPeak >\${file_basename}_new_sorted.broadPeak
    done

    # now for the narrow peaks

    narrow_list=(${roadmap_narrowhistones_path.join(' ')})

    for file in \${narrow_list[@]}; do

        file_basename=\$(basename \${file} .narrowPeak)
        awk 'BEGIN{OFS="\t"} \$3 > \$2 {print \$1, \$2, \$3 }' \${file} > \${file_basename}_new.narrowPeak

        sort -k1,1 -k2,2n \${file_basename}_new.narrowPeak >\${file_basename}_new_sorted.narrowPeak
    done

    # now for the pipeline peaks

    pipeline_list=(${file_list.join(' ')})

    for file in \${pipeline_list[@]}; do

        if [ \$(wc -l < \${file}) -gt 0 ]; then
            file_basename=\$(basename \${file} .bed)
            awk 'BEGIN{OFS="\t"} \$3 > \$2 {print \$1, \$2, \$3 }' \${file} > \${file_basename}_new.bed

            sort -k1,1 -k2,2n \${file_basename}_new.bed >\${file_basename}_new_sorted.bed
        fi
    done



    broadPeak_files=\$(ls *new_sorted.broadPeak)

    # here I will find the counts for all of the bedfiles per each atac-seq condition


    plotEnrichment \
    --perSample \
    --bamfiles ${bam_files[1]} \
    --BED \${broadPeak_files[@]} \
    --variableScales \
    --outRawCounts ${out_broad_enrich_counts} \
    --plotFile ${out_broad_enrich_plot}

    # getting the plots and counts for narrow peaks

    narrowPeak_files=\$(ls *new_sorted.narrowPeak)

    # here I will find the counts for all of the bedfiles per each atac-seq condition


    plotEnrichment \
    --perSample \
    --bamfiles ${bam_files[1]} \
    --BED \${narrowPeak_files[@]} \
    --variableScales \
    --outRawCounts ${out_narrow_enrich_counts} \
    --plotFile ${out_narrow_enrich_plot}


    # now plotting for the peaks from our own data

    pipelinePeak_files=\$(ls *new_sorted.bed)

    plotEnrichment \
    --perSample \
    --bamfiles ${bam_files[1]} \
    --BED \${pipelinePeak_files[@]} \
    --variableScales \
    --outRawCounts ${out_pipeline_enrich_counts} \
    --plotFile ${out_pipeline_enrich_plot}



    """
}

process r_atac_enrich_plot_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/r_language'
    label 'normal_big_resources'
    publishDir "./enrichment_of_experiment/", mode: 'copy', pattern:'*'

    input:

    // this is where the enrichment tab meta channel will be 
    // the basename, filename, and peakpath will have two elements(conditions) the list, scrm and h1low counts in histone peaks
    // just plot both h1low/scrm and scrm/h1low
    tuple val(peak_type), val(basename), val(filename), path(peakpath)


    output:

    path("*.png"), emit: atac_enrichment_plots


    script:

    condition_one = filename[0]

    condition_two = filename[1]

    // i need to automate the out file names
    atac_enrich_png_v1 = "atac_enrichment_in_${peak_type}_histones_v1.png"
    atac_enrich_png_v2 = "atac_enrichment_in_${peak_type}_histones_v2.png"
    

    """
    #!/usr/bin/env Rscript

    # now just read the tab files into R

    # i have to load the readr package to use read_*

    library(readr)
    library(ggplot2)

    condition_one = read_tsv("./${condition_one}")
    condition_two = read_tsv("./${condition_two}")

    atac_narrow_enrichment_counts_v1 = data.frame(featureType = condition_one\$featureType, log2FC = log2(condition_one\$percent/condition_two\$percent) )

    atac_narrow_enrichment_counts_v2 = data.frame(featureType = condition_one\$featureType, log2FC = log2(condition_two\$percent/condition_one\$percent) )

    # using ifelse is a conditional that does this 
    # setting another column for both the versions of calculating
    atac_narrow_enrichment_counts_v1\$color = ifelse( atac_narrow_enrichment_counts_v1\$log2FC > 0, 'up', 'down')

    atac_narrow_enrichment_counts_v2\$color = ifelse( atac_narrow_enrichment_counts_v2\$log2FC > 0, 'up', 'down')


    # only need to get the new names one time

    print(atac_narrow_enrichment_counts_v1\$featureType)

    string_test = strsplit(atac_narrow_enrichment_counts_v1\$featureType, "_")

    new_names = sapply(string_test, function(x) {
        paste(na.omit(x[c(1:4,12,13)]), collapse = "_")
    })

    new_names

    #new_name = paste(na.omit(string_test[[1]][1:5-12-13]), collapse = "_")
    #new_name

    # now here is where I change the featureType names to the new names that should be shorter, and do it in both versions of the experiment design
    atac_narrow_enrichment_counts_v1["featureType"] = new_names
    atac_narrow_enrichment_counts_v2["featureType"] = new_names

    print(atac_narrow_enrichment_counts_v1\$featureType)
                                         

    ggplot(data = atac_narrow_enrichment_counts_v1, aes(x = featureType, y = log2FC, fill = color) )+
        geom_bar(stat = "identity")+
        scale_fill_manual(values=c( "up" = "green", "down" = "gray"))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9), plot.caption = element_text(hjust = 0, size = 9),
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10) )+
        labs(title = paste0('ATAC-seq ', "${peak_type} ", 'log2FC ', condition_one\$file[1], '_vs_', condition_two\$file[1]), caption = "This is k562 ChIP-seq narrow peak data from roadmap. \nDeeptools was used to generate counts from ATAC-seq H1low and scrambled bam files,\n to get the number of reads that align in each of the roadmap narrow peaks. The percentage of \nreads in peaks vs total reads is the column I use to get the accessibility enrichment score.\n I take log2 of the percentage of H1low reads divided by the percentage of Scrm reads ratio,\n which gives the log fold enrichment of accessibility")+
        scale_y_continuous(labels = scales::label_number())

    ggsave("${atac_enrich_png_v1}", plot=last_plot(), width = 12, height = 8, units = "in", dpi = 300)


    # now version 2

    ggplot(data = atac_narrow_enrichment_counts_v2, aes(x = featureType, y = log2FC, fill = color) )+
        geom_bar(stat = "identity")+
        scale_fill_manual(values=c( "up" = "green", "down" = "gray"))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9), plot.caption = element_text(hjust = 0, size = 9),
        plot.margin = margin(t = 10, r = 10, b = 30, l = 10) )+
        labs(title = paste0('ATAC-seq ', "${peak_type} ", 'log2FC ', condition_two\$file[1], '_vs_',condition_one\$file[1]), caption = "This is k562 ChIP-seq narrow peak data from roadmap. \nDeeptools was used to generate counts from ATAC-seq H1low and scrambled bam files,\n to get the number of reads that align in each of the roadmap narrow peaks. The percentage of \nreads in peaks vs total reads is the column I use to get the accessibility enrichment score.\n I take log2 of the percentage of H1low reads divided by the percentage of Scrm reads ratio,\n which gives the log fold enrichment of accessibility")+
        scale_y_continuous(labels = scales::label_number())

    ggsave("${atac_enrich_png_v2}", plot=last_plot(), width = 12, height = 8, units = "in", dpi = 300)





    """
}

process get_atacPeaks_in_roadmapPeaks_process {

    // this process should be using bedtools intersect
    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/bedtools_rj'
    label 'normal_big_resources'

    publishDir "./intersect_atacPeaks_in_histone_peaks", mode: 'copy', pattern: '*'
    
    input:

    path(nochange_atac_peaks)
    path(up_atac_peaks)
    path(broad_peak)
    path(narrow_peak)
    //tuple val(broad_name), val(broad_basename), path(roadmap_broad_histone_file)
    //tuple val(narrow_name), val(narrow_basename), path(roadmap_narrowhistone_file)




    output:

    path("*.bed"), emit: all_atac_peak_intersections



    script:

    out_up_broad_intersect_file = "up_atacPeaks_${broad_peak}.bed"
    out_up_narrow_intersect_file = "up_atacPeaks_${narrow_peak}.bed"

    out_nochange_broad_intersect_file = "nochange_atacPeaks_${broad_peak}.bed"
    out_nochange_narrow_intersect_file = "nochange_atacPeaks_${narrow_peak}.bed"

    

    """
    #!/usr/bin/env bash

    # first find the up atac peaks in the broad

    bedtools intersect \
    -a ${up_atac_peaks} \
    -b ${broad_peak} \
    -wa \
    -f 0.20 \
    > ${out_up_broad_intersect_file}

    # second find up atac peaks in narrow

    bedtools intersect \
    -a ${up_atac_peaks} \
    -b ${narrow_peak} \
    -wa \
    -f 0.20 \
    > ${out_up_narrow_intersect_file}

    # third find nochange atac peaks in broad

    bedtools intersect \
    -a ${nochange_atac_peaks} \
    -b ${broad_peak} \
    -wa \
    -f 0.20 \
    > ${out_nochange_broad_intersect_file}

    # fourth find nochange atac peaks in narrow
    bedtools intersect \
    -a ${nochange_atac_peaks} \
    -b ${narrow_peak} \
    -wa \
    -f 0.20 \
    > ${out_nochange_narrow_intersect_file}



    """
}
