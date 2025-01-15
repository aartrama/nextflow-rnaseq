nextflow.enable.dsl=2

// Set parameters
project_dir = params.project_dir
index_basename = params.index_basename
input_dir = params.input_dir
gtf_file = params.gtf_file
pairedEnd = params.pairedEnd

// Define the input channel 
read_pairs_ch = Channel.fromFilePairs("${input_dir}/*_{R1,R2,1,2}{,_001,_S[0-9]+}{,_001}{,.fastq,.fq}{,.gz}", checkIfExists: true)

process FASTQC {
    publishDir "$project_dir/output/fastqc", mode: 'copy'
    cpus 4

    input:
    tuple val(pair_id), path(reads)

    output:
    path("*.html"), emit: html_files
    path("*.zip"), emit: zip_files

    script:
    """
    fastqc ${reads} --outdir .
    """ 

}


process READ_TRIM_GALORE {

    publishDir "$project_dir/output/trim_galore", mode: 'copy'
    cpus 4
    
    input:
    tuple val(pair_id), path(reads)
    
    output:
    tuple val(pair_id), path("*_val_1.fq.gz"), path("*_val_2.fq.gz"), emit: trimmed_reads
    
    script:

    paired_end = !pairedEnd ? "" : " --paired"
    """
    trim_galore ${reads} ${paired_end}
    """
}


process ALIGNMENT_STEP  {
    publishDir "$project_dir/output/bam", mode: 'copy'
    cpus 4

    input:
    tuple val(pair_id), path(trimmed_read1), path(trimmed_read2)
 
    output:
    tuple val(pair_id), path("${pair_id}.sam"), emit: sam
    path("${pair_id}.txt"), emit: log

    script:
    """
    hisat2 -p $task.cpus -x ${index_basename} -1 ${trimmed_read1} -2 ${trimmed_read2} -S ${pair_id}.sam --summary-file ${pair_id}.txt --temp-directory \${PWD}
    """
}

process SORTED_BAM_STEP {
    publishDir "$project_dir/output/bam", mode: 'copy'
    cpus 4

     
    input:
    tuple val(pair_id), path(sam_file)

    output:
    tuple val(pair_id), path("${pair_id}.sorted.bam"), emit: bam
 
    script:
    """
    samtools sort -@ ${task.cpus} -O BAM -o ${pair_id}.sorted.bam $sam_file
    """
}

process COUNTS_STEP {
    publishDir "$project_dir/output/counts", mode: 'copy'
    cpus 4

    input:
    tuple val(pair_id), path(sorted_bam_file)

    output:
    path("${pair_id}.featureCounts.txt"), emit: counts
 
    script:

    paired_end = !pairedEnd ? "" : "-p"

    """
    featureCounts ${paired_end} -t exon -a ${gtf_file} -o ${pair_id}.featureCounts.txt $sorted_bam_file 
    """
}

process COUNTS_MATRIX {
    
    publishDir "$project_dir/output/counts", mode: 'copy'
    cpus 4

    input:
    path(count_files)

    output:
    path("featureCounts.txt")

    script:
    """
    #!/usr/bin/env python

    import sys
    import os 
    import pandas as pd
    dict_of_counts = {}
    
    for files in os.listdir("$project_dir/output/counts"):
        sample = files.split(".")[0]
        print(sample)
        dict_of_counts[sample] = {}
        with open("$project_dir/output/counts/"+files, "r") as infile:
            next(infile)
            next(infile)
            for lines in infile:
                lines = lines.strip().split("\\t")
                dict_of_counts[sample][lines[0]] = int(float(lines[-1]))
    
    dataframe = pd.DataFrame(dict_of_counts)
    dataframe.to_csv("featureCounts.txt", sep='\\t')
    """


}

process MULTIQC_STEP {
    publishDir "$project_dir/output", mode: 'copy'
    cpus 4

    input:
    path("output")

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"), emit: data

    script:
    """
    multiqc ${project_dir}
    """
}

// Define whole workflow
workflow {
    fastqc_files_ch = FASTQC(read_pairs_ch)
    trimmed_files_ch = READ_TRIM_GALORE(read_pairs_ch)
    sam_files_ch = ALIGNMENT_STEP(trimmed_files_ch.trimmed_reads)
    sorted_bam_files_ch = SORTED_BAM_STEP(sam_files_ch.sam)
    count_files_ch = COUNTS_STEP(sorted_bam_files_ch.bam)
    counts_matrix = COUNTS_MATRIX(count_files_ch)
    multiqc_report_file = MULTIQC_STEP(count_files_ch)
}


