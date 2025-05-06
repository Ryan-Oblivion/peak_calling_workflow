

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
        --scaleFactor 1 

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
        --scaleFactor 1 

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
    --sortUsing sum







    """
}



process macs2_call_peaks_process_both {

    debug true

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


process find_idr_in_replicates_process {

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/idr-2.0_rj'
    label 'super_big_resources'
    publishDir "idr_results/${histone[0]}/${condition[0]}", mode: 'copy', pattern:'*'


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
    --heatmapWidth 8 \
    --heatmapHeight 28 \
    --labelRotation 30 \
    --samplesLabel ${true_bw_name} \
    --regionsLabel "${true_bw_name} called peaks"

    plotProfile \
    -m "${out_matrix_name}" \
    --outFileName "${profile_out_name}" \
    --plotWidth 8 \
    --plotHeight 10 \
    --labelRotation 30 \
    --samplesLabel ${true_bw_name}






    """
}