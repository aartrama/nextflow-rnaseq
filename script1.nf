nextflow.enable.dsl=2

// Set parameters
index_path = params.index_path
index_basename = params.index_basename
read_pairs = params.read_pairs 

// Print parameters
log.info"""\
	index_path: ${index_path}
	index_basename: ${index_basename}
	read_pairs	: ${read_pairs}
	"""

// Define the input channel
read_pairs_ch = Channel.fromFilePairs(params.read_pairs, checkIfExists: true)

process quantification {

	cpus 4
     
    input:
    tuple val(pair_id), path(reads)
 
    output:
    path("${pair_id}.bam"), emit: bam
    path("${pair_id}.txt"), emit: log
 
    script:
    """
    hisat2 -p $task.cpus -x ${index_path}/${index_basename} -1 ${reads[0]} -2 ${reads[1]} -S ${pair_id}.bam 2> ${pair_id}.txt
    """
}


// Define whole workflow
workflow {
	index_ch = quantification(read_pairs_ch)
}
