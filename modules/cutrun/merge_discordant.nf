process MERGE_DISCORDANT_CUTRUN {
    publishDir "${project_dir}/output/bam", mode: 'copy'
    cpus 1
     
    input:
    tuple val(pair_id), path(clip_bam)
    tuple val(pair_id_2), path(mapped_bam), path(dis_bam)

    output:
    tuple val(pair_id), path("${pair_id}.sorted.bam"), path("${pair_id}.sorted.bam.bai"), emit: bam
 
    script:
    """
    ## combine non-discord and filtered discordant pairs
    samtools merge -f -@ ${task.cpus} ${pair_id}_filtered.bam \
                ${mapped_bam} \
                ${clip_bam}
    ## sort combined bam file
    samtools sort -@ ${task.cpus} -O BAM -o ${pair_id}.sorted.bam ${pair_id}_filtered.bam
    samtools index ${pair_id}.sorted.bam ${pair_id}.sorted.bam.bai
    """

}
