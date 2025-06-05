process CLIP_DISCORDANT_SAMTOOL_CUTRUN {
    publishDir "${project_dir}/output/bam", mode: 'copy'     
    input:
    tuple val(pair_id), path(dis_bam)

    output:
    tuple val(pair_id), path("${pair_id}_mapped_dis_sorted.bam"), emit: mapped_dis_bam_sorted
 
    script:
    """
    ## remove overlapping bases from discordant pairs
    samtools sort -@ 1 -O BAM -o ${pair_id}_mapped_dis_sorted.bam ${dis_bam}
    """
}

process CLIP_DISCORDANT_BAMUTIL_CUTRUN {
    publishDir "${project_dir}/output/bam", mode: 'copy'         
    input:
    tuple val(pair_id), path(mapped_dis_bam)

    output:
    tuple val(pair_id), path("${pair_id}_clipped_dis.bam"), emit: bam 
 
    script:
    """ 
    bam clipOverlap --in ${mapped_dis_bam} --out ${pair_id}_clipped_dis.bam --stats --overlapsOnly
    """ 
}
