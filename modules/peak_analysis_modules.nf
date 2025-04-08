

process make_alignment_bw_process_control {

    label 'normal_small_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'

    publishDir "./${condition_label[0]}/${histone_label}_bigwigs/", mode: 'copy', pattern: "*${histone_label}*.bigwig"


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

    rep1_file = "${bam_file_name[0]}"
    rep2_file = "${bam_file_name[1]}"
    rep3_file = "${bam_file_name[2]}"

    // rep1_meta_name = "${meta_name[0]}"
    // rep2_meta_name = "${meta_name[1]}"
    // rep3_meta_name = "${meta_name[2]}"


    // rep1_k27_file = "${bam_file_name_k27[2]}"
    // rep2_k27_file = "${bam_file_name_k27[1]}"
    // rep3_k27_file = "${bam_file_name_k27[0]}"

    // rep1_k9_file = "${bam_file_name_k9[2]}"
    // rep2_k9_file = "${bam_file_name_k9[1]}"
    // rep3_k9_file = "${bam_file_name_k9[0]}"






    """
    #!/usr/bin/env bash

    ###### Deeptools parameters ##########




    #######################################


    # making a bigwig file for each bam file
    for bam in "${rep1_file}" "${rep2_file}" "${rep3_file}"; do
        
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

// exact same process as above but for wt not control

process make_alignment_bw_process_wt {

    label 'normal_small_resources'

    conda '/ru-auth/local/home/rjohnson/miniconda3/envs/deeptools_rj'

    publishDir "./${condition_label[0]}/${histone_label}_bigwigs/", mode: 'copy', pattern: "*${histone_label}*.bigwig"


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

    rep1_file = "${bam_file_name[0]}"
    rep2_file = "${bam_file_name[1]}"
    rep3_file = "${bam_file_name[2]}"


    // rep1_k27_file = "${bam_file_name_k27[2]}"
    // rep2_k27_file = "${bam_file_name_k27[1]}"
    // rep3_k27_file = "${bam_file_name_k27[0]}"

    // rep1_k9_file = "${bam_file_name_k9[2]}"
    // rep2_k9_file = "${bam_file_name_k9[1]}"
    // rep3_k9_file = "${bam_file_name_k9[0]}"






    """
    #!/usr/bin/env bash

    ###### Deeptools parameters ##########




    #######################################


    # making a bigwig file for each bam file
    for bam in "${rep1_file}" "${rep2_file}" "${rep3_file}"; do
        
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

    """
    #!/usr/bin/env bash

    ########### deeptools params ##########



    #######################################

    echo ' this is the list of names: "${name_list}"'

    # fix the genebody file
    awk  '\$2!=\$3 {print \$0}' "${wtvslowup_genebody}" > wtvslowup_genebody_noZerolength.bed

    computeMatrix scale-regions -S ${bw_names.join(' ')} \
    -R wtvslowup_genebody_noZerolength.bed \
    --beforeRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 3000 \
    --skipZeros \
    --numberOfProcessors "max" \
    -o "${out_matrix_scores_upgenes}"

    
    plotHeatmap -m "${out_matrix_scores_upgenes}" \
    -out "${png_heatmap_upgenes}" \
    --samplesLabel ${name_list.join(' ')} \
    --labelRotation 30 \
    --heatmapWidth 8 \
    --sortUsing sum



    # now to do this with the other genebody file 'down unchanging'

    awk  '\$2!=\$3 {print \$0}' "${wtvslowdown_nochange_genebody}" > wtvslowdown_nochange_genebody_noZerolength.bed

    computeMatrix scale-regions -S ${bw_names.join(' ')} \
    -R wtvslowdown_nochange_genebody_noZerolength.bed \
    --beforeRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 3000 \
    --skipZeros \
    --numberOfProcessors "max" \
    -o "${out_matrix_scores_down_unchange_genes}"

    
    plotHeatmap -m "${out_matrix_scores_down_unchange_genes}" \
    -out "${png_heatmap_down_unchanging_genes}" \
    --samplesLabel ${name_list.join(' ')} \
    --labelRotation 30 \
    --heatmapWidth 8 \
    --sortUsing sum



    # now plotting both together

    computeMatrix scale-regions -S ${bw_names.join(' ')} \
    -R wtvslowup_genebody_noZerolength.bed wtvslowdown_nochange_genebody_noZerolength.bed \
    --beforeRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 3000 \
    --skipZeros \
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


process macs2_call_peaks_process_control {

    label 'normal_big_resources'

    conda ''

    publishDir " "



    input:



    output:



    script:


    """
    #!/usr/bin/env bash




    """
}