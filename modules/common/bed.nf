process CREATE_BED {
    publishDir "${project_dir}/output/bed", mode: 'copy'
    input:
    tuple val(pair_id), path(bam_file), path(bai_file)

    output:
    path("${pair_id}.unique.sorted.rmdup.chr.bed"), emit: bed

    script:
    // Create bed file from bam file (main chromosomes only)
    """
    bedtools bamtobed -i ${bam_file} > ${pair_id}.unique.sorted.rmdup.chr.bed 2> ${pair_id}_bed.log
    """
}
