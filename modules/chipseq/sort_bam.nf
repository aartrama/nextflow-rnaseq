process SORT_BAM_CHIPSEQ {
    publishDir "${project_dir}/output/bam", mode: 'copy'
    cpus 1
     
    input:
    tuple val(pair_id), path(sam_file)

    output:
    tuple val(pair_id), path("${pair_id}.sorted.bam"), path("${pair_id}.sorted.bam.bai"), emit: bam
 
    script:
    paired_specific_command = !pairedEnd ? "" : "-f 3 -F 8"
    """
    samtools view -bh -F 4 -F 256 ${paired_specific_command} ${sam_file} > ${pair_id}_filtered.bam
    samtools sort -@ ${task.cpus} -O BAM -o ${pair_id}.sorted.bam ${pair_id}_filtered.bam
    samtools index ${pair_id}.sorted.bam ${pair_id}.sorted.bam.bai
    """
}
