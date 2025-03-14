process ALIGNMENT_STEP  {
    publishDir "${project_dir}/output/bam", mode: 'copy', pattern: "*.txt"
     
    input:
    tuple val(pair_id), path(trimmed_reads)
 
    output:
    tuple val(pair_id), path("${pair_id}.sam"), emit: sam
    path("${pair_id}.txt"), emit: log
 
    script:
    sample_names = !pairedEnd ? "-U ${trimmed_reads[0]}" : "-1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]}"

    """
    hisat2 -p ${task.cpus} -x ${index_files}/${index_basename} ${sample_names} -S ${pair_id}.sam --summary-file ${pair_id}.txt --temp-directory \${PWD} 
    """
}
