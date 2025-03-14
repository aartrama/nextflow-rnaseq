process CREATE_BW {
    publishDir "${project_dir}/output/bw", mode: 'copy'
    
    input:
    tuple val(pair_id), path(bam_file), path(bai_file)

    output:
    path("${pair_id}.unique.sorted.rmdup.chr.bw")

    script:
    """
    bamCoverage -b ${bam_file} -o ${pair_id}.unique.sorted.rmdup.chr.bw --binSize 10 --normalizeUsing RPKM
    """
}
