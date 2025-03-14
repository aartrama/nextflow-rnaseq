process CLIP_DISCORDANT_CUTRUN {
    publishDir "${project_dir}/output/bam", mode: 'copy'     
    input:
    tuple val(pair_id), path(mapped_bam), path(dis_bam)

    output:
    tuple val(pair_id), path("${pair_id}_clipped_dis.bam"), emit: bam
 
    script:
    """
    ## remove overlapping bases from discordant pairs
    samtools sort -@ ${task.cpus} -O BAM -o ${pair_id}_mapped_dis_sorted.bam ${dis_bam}
    bam clipOverlap --in ${pair_id}_mapped_dis_sorted.bam --out ${pair_id}_clipped_dis.bam --stats --overlapsOnly
    """
}
