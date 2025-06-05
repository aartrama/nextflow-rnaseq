process READ_TRIM_UMI{
    publishDir "${project_dir}/output/trim_fastq", mode: 'copy', pattern: "_trimming_report.txt"

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*.fq.gz"), emit: trimmed_reads
 
    script:
    if( !params.pairedEnd )
        """
        umi_tools extract --stdin=${reads} --bc-pattern=${params.umi_1} --stdout ${pair_id}_R1_trimmed.fq.gz
        """
    else
        """
        umi_tools extract -I ${reads[0]} --bc-pattern=${params.umi_1} --bc-pattern2=${params.umi_2} --read2-in=${reads[1]} --stdout=${pair_id}_R1_trimmed.fq.gz --read2-out=${pair_id}_R2_trimmed.fq.gz
        """
}

process UMI_RMDUPS_BAM {
    publishDir "${project_dir}/output/bam", mode: 'copy'

    input:
    tuple val(pair_id), path(bam_file), path(bai_file)

    output: 
    tuple val(pair_id), path("${pair_id}.unique.sorted.rmdup.bam"), emit: bam
    tuple val(pair_id), path("${pair_id}.rmdup.log"), emit: qc

    script:
    pe =  params.pairedEnd ? "--paired --unmapped-reads=use" : ""
    """
    umi_tools dedup --stdin=$bam_file --log=${pair_id}.rmdup.log $pe > ${pair_id}.unique.sorted.rmdup.bam
    """
}
