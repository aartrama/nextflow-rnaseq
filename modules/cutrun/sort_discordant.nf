process SORT_DISCORDANT_CUTRUN {
    publishDir "${project_dir}/output/bam", mode: 'copy'
    input:
    tuple val(pair_id), path(bam_file)

    output:
    tuple val(pair_id), path("${pair_id}_mapped_dis.bam"), emit: mapped_dis_bam
    tuple val(pair_id), path("${pair_id}_mapped_paired.bam"), emit: mapped_paired_bam
 
    script:
    """
    ## Split file mapped pairs into discordant/cordant components
    samtools view -bh -f 3 -F 4 -F 8 -F 256 ${bam_file} > ${pair_id}_mapped_paired.bam
    samtools view -bh -f 1 -F 2 -F 4 -F 8 -F 256 ${bam_file} > ${pair_id}_mapped_dis.bam
    """

}
