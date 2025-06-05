process SORTED_BAM_STEP {
    publishDir "$project_dir/output/bam", mode: 'copy'

     
    input:
    tuple val(pair_id), path(sam_file)

    output:
    tuple val(pair_id), path("${pair_id}.sorted.bam"), emit: bam
 
    script:
    """
    samtools sort -@ ${task.cpus} -O BAM -o ${pair_id}.sorted.bam $sam_file
    """
}

process PICARD_RMDUPS_BAM {
    publishDir "${project_dir}/output/bam", mode: 'copy'
    input:
    tuple val(pair_id), path(bam_file), path(bai_file)

    output: 
    tuple val(pair_id), path("${pair_id}.unique.sorted.rmdup.bam"), emit: bam
    tuple val(pair_id), path("${pair_id}.rmdup.log"), emit: log
    tuple val(pair_id), path("${pair_id}_marked_dup_metrics.txt"), emit: qc

    script:
    stringency = params.experiment_type == "cutrun" ? "" : "-VALIDATION_STRINGENCY SILENT "
    if ( params.experiment_type == "cutrun" )
    """
    picard MarkDuplicates REMOVE_DUPLICATES=true I=${bam_file} O=${pair_id}.unique.sorted.rmdup.bam M=${pair_id}_marked_dup_metrics.txt VALIDATION_STRINGENCY=SILENT 2> ${pair_id}.rmdup.log
    """
    else if ( params.experiment_type != "cutrun" )
    """
    picard MarkDuplicates REMOVE_DUPLICATES=true I=${bam_file} O=${pair_id}.unique.sorted.rmdup.bam M=${pair_id}_marked_dup_metrics.txt 2> ${pair_id}.rmdup.log
    """

}


process CALC_INSERT_SIZE {
    publishDir "${project_dir}/output/logs", mode: 'copy'
    input:
    tuple val(pair_id), path(bam_file)
    
    output:
    tuple val(pair_id), path("${pair_id}_insert_size_metrics.txt"), path("${pair_id}_insert_size_histogram.pdf"), emit: log
    tuple val(pair_id), path("${pair_id}_insert_size_metrics.txt"), emit: qc

    script:
    // Calculate insert size (only for paired-end libraries) 
    """
    picard CollectInsertSizeMetrics \
        I=${bam_file} \
        O=${pair_id}_insert_size_metrics.txt \
        H=${pair_id}_insert_size_histogram.pdf \
        M=0.05 VALIDATION_STRINGENCY=SILENT
    """
}


process RMDUP_TO_CHRBAM {
    publishDir "${project_dir}/output/bam", mode: 'copy'
   
    input:
    tuple val(pair_id), path(bam_file)

    output:
    tuple val(pair_id), path("${pair_id}.unique.sorted.rmdup.chr.bam"), path("${pair_id}.unique.sorted.rmdup.chr.bam.bai"), emit: bam

    script:
    // Eliminate extraneous chromosomes from bam file 
    """ 
    samtools view -H ${bam_file} | sed -e 's/SN:\\([0-9XY]\\)/SN:chr\\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - ${bam_file} > ${pair_id}.unique.sorted.rmdup.chr.bam
    samtools index ${pair_id}.unique.sorted.rmdup.chr.bam ${pair_id}.unique.sorted.rmdup.chr.bam.bai
    """ 
}



