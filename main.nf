nextflow.enable.dsl=2

// Set parameters
project_dir = params.project_dir
index_basename = params.index_basename
input_dir = "$project_dir"+'/fastq'
gtf_file = params.gtf_file
pairedEnd = params.pairedEnd
multiqc_config = params.multiqc_config
count_unique = params.count_unique
count_fraction = params.count_fraction

// Define the input channel 
if (pairedEnd) {
    read_pairs_ch = Channel.fromFilePairs("${input_dir}/*_{R1,R2,1,2}{,_001,_S[0-9]+}{,_001}{,.fastq,.fq}{,.gz}", checkIfExists: true)
} else {
    read_pairs_ch = Channel.fromPath("${input_dir}/*.{fastq,fq}{,.gz}", checkIfExists: true)
        .map { file -> tuple(file.simpleName, file) }
}


process FASTQC {
    publishDir "$project_dir/output/fastqc", mode: 'copy'

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
    
    input:
    tuple val(pair_id), path(reads)
    
    output:
    tuple val(pair_id), path("*fq.gz"), emit: trimmed_reads
    
    script:

    paired_end = !pairedEnd ? "" : " --paired"
    """
    trim_galore ${reads} ${paired_end}
    """
}


process ALIGNMENT_STEP  {
    publishDir "$project_dir/output/bam", mode: 'copy'

    input:
    tuple val(pair_id), path(trimmed_reads)
 
    output:
    tuple val(pair_id), path("${pair_id}.sam"), emit: sam
    path("${pair_id}.txt"), emit: log

    script:

    paired_end = !pairedEnd ? "-U ${trimmed_reads[0]}" : "-1 ${trimmed_reads[0]} -2 ${trimmed_reads[1]}"

    """
    hisat2 -p $task.cpus -x ${index_basename} ${paired_end} -S ${pair_id}.sam --summary-file ${pair_id}.txt --temp-directory \${PWD}
    """
}

process SORTED_BAM_STEP {
    publishDir "$project_dir/output/bam", mode: 'copy'

     
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

    input:
    tuple val(pair_id), path(sorted_bam_file)

    output:
    path("${pair_id}.exon.txt"), emit: counts_exonic
    path("${pair_id}.gene.txt"), emit: counts_genic
 
    script:

    paired_end = !pairedEnd ? "" : "-p"
    count_scheme = count_unique ? "" : (count_fraction ? "-O --fraction" : "-O")

    """
    featureCounts ${paired_end} -t exon ${count_scheme} -a ${gtf_file} -o ${pair_id}.exon.txt $sorted_bam_file 
    featureCounts ${paired_end} -t gene ${count_scheme} -a ${gtf_file} -o ${pair_id}.gene.txt $sorted_bam_file
    """
}

process COUNTS_MATRIX {
    
    publishDir "$project_dir/output/counts", mode: 'copy'

    input:
    path("output/counts/*txt")

    output:
    path("featureCounts_genic.txt"), emit: count_matrix_genic
    path("featureCounts_exonic.txt"), emit: count_matrix_exonic

    script:
    """
    #!/usr/bin/env python

    import sys
    import os 
    import pandas as pd

    
    gene_counts = {}
    exon_counts = {}
    
    for files in "output/counts/*txt":
        if "gene.txt" in files:
            counts_dict = gene_counts
        elif "exon.txt" in files:
            counts_dict = exon_counts
        else:
            continue
            
        sample = files.split(".")[0] 
        counts_dict[sample] = {}
        with open(files, "r") as infile:
            next(infile)
            next(infile)
            for lines in infile:
                lines = lines.strip().split("\\t")
                counts_dict[sample][lines[0]] = int(float(lines[-1]))
    
    gene_dataframe = pd.DataFrame(gene_counts)
    exon_dataframe = pd.DataFrame(exon_counts)
    gene_dataframe.to_csv("featureCounts_genic.txt", sep='\\t')
    exon_dataframe.to_csv("featureCounts_exonic.txt", sep='\\t')

    """


}


process MULTIQC_STEP {
    publishDir "$project_dir/output/multiqc", mode: 'copy'

    input:
    path("output")

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"), emit: data
    path("multiqc_summary_text.txt")

    script:
    """
    multiqc ${project_dir} --config ${multiqc_config} --force &> multiqc_summary_text.txt
    """
}

// Define whole workflow
workflow {
    fastqc_files_ch = FASTQC(read_pairs_ch)
    trimmed_files_ch = READ_TRIM_GALORE(read_pairs_ch)
    sam_files_ch = ALIGNMENT_STEP(trimmed_files_ch.trimmed_reads)
    sorted_bam_files_ch = SORTED_BAM_STEP(sam_files_ch.sam)
    count_files_ch = COUNTS_STEP(sorted_bam_files_ch.bam)
    counts_matrix = COUNTS_MATRIX(count_files_ch.counts_exonic)
    multiqc_report_file = MULTIQC_STEP(counts_matrix.count_matrix_genic)
}


