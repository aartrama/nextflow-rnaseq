process COUNTS_STEP {
    publishDir "$project_dir/output/counts", mode: 'copy'

    input:
    tuple val(pair_id), path(sorted_bam_file)

    output:
    path("${pair_id}.exon.txt"), emit: counts_exonic
    path("${pair_id}.gene.txt"), emit: counts_genic
 
    script:

    paired_end = !params.pairedEnd ? "" : "-p"
    count_scheme = params.count_unique ? "" : (params.count_fraction ? "-O --fraction" : "-O")

    """
    featureCounts ${paired_end} -t exon ${count_scheme} -a ${params.gtf_file} -o ${pair_id}.exon.txt $sorted_bam_file 
    featureCounts ${paired_end} -t gene ${count_scheme} -a ${params.gtf_file} -o ${pair_id}.gene.txt $sorted_bam_file
    """
}
