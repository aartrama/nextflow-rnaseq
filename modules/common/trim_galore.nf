process READ_TRIM_GALORE_BASIC {

    publishDir "${project_dir}/output/trim_fastq", mode: 'copy', pattern: "*_trimming_report.txt"

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*_trimmed.fastq.gz"), emit: trimmed_reads
    tuple path("*_trimming_report.txt"), emit: report

    script:
    def paired_end = pairedEnd ? / --paired / : ''
    def cr_status = experiment_type == "cutrun" ? / --clip_R1 6 --clip_R2 6 / :  ''

    if( pairedEnd )
        """
        trim_galore ${cr_status} ${paired_end} --trim-n ${reads}
        mv ${pair_id}_1_val_1.fq.gz ${pair_id}_1_trimmed.fastq.gz
        mv ${pair_id}_2_val_2.fq.gz ${pair_id}_2_trimmed.fastq.gz
        """
    
    else if ( !pairedEnd )
        """
        trim_galore ${cr_status} ${paired_end} --trim-n ${reads}
        mv ${pair_id}_trimmed.fq.gz ${pair_id}_trimmed.fastq.gz
        """
}


process READ_TRIM_GALORE_POLYA {

    publishDir "${project_dir}/output/trim_fastq", mode: 'copy', pattern: "*_trimming_report.txt"

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("*fq.gz"), emit: trimmed_reads

    script:
    def paired_end = pairedEnd ? / --paired / : ""
    def output_files = pairedEnd ? 
            "${pair_id}_1_trimmed.fastq.gz ${pair_id}_2_trimmed.fastq.gz" :
            "${pair_id}_trimmed.fastq.gz"
    """
    trim_galore ${paired_end} --polyA ${reads} 
    """
}
