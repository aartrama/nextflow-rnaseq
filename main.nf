nextflow.enable.dsl=2

// Set parameters
index_basename = params.index_basename
input_dir = params.input_dir
gtf_file = params.gtf_file
pairedEnd = params.pairedEnd
temp_dir = params.temp_dir
// Print parameters
log.info"""\
	index_basename: ${index_basename}
    input_dir: ${input_dir}
	gtf_file: ${gtf_file}
	"""

// Define the input channel
read_pairs_ch = Channel.fromFilePairs("${input_dir}/*_R{1,2}_001.{fastq,fq}{,.gz}", 
                                        checkIfExists: true)

process ALIGNMENT_STEP  {
    cpus 4

    input:
    tuple val(pair_id), path(reads)
 
    output:
    tuple val(pair_id), path("${pair_id}.sam"), emit: sam
    path("${pair_id}.txt"), emit: log

    script:
    """
    hisat2 -p $task.cpus -x ${index_basename} -1 ${reads[0]} -2 ${reads[1]} -S ${pair_id}.sam --summary-file ${pair_id}.txt --temp-directory ${temp_dir}
    """
}

process SORTED_BAM_STEP {

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



// Define whole workflow
workflow {
    sam_files_ch = ALIGNMENT_STEP(read_pairs_ch)
    sorted_bam_files_ch = SORTED_BAM_STEP(sam_files_ch.sam)
    count_files = COUNTS_STEP(sorted_bam_files_ch.bam)
}


