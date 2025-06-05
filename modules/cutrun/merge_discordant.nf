process MERGE_DISCORDANT_CUTRUN {
    publishDir "${project_dir}/output/bam", mode: 'copy'
    cpus 1
     
    input:
    tuple val(pair_id), path(mapped_bam) 
    tuple val(pair_id), path(clipped_dis_bam)

    output:
    tuple val(pair_id), path("${pair_id}.sorted.bam"), path("${pair_id}.sorted.bam.bai"), emit: bam
 
    script:
    """
    ## combine non-discord and filtered discordant pairs
    samtools merge -f -@ 1 ${pair_id}_filtered.bam \
                ${mapped_bam} \
                ${clipped_dis_bam}
    ## sort combined bam file
    samtools sort -@ 1 -O BAM -o ${pair_id}.sorted.bam ${pair_id}_filtered.bam
    samtools index ${pair_id}.sorted.bam ${pair_id}.sorted.bam.bai
    """

}
