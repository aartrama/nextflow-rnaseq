log.info"""\
	index_path: ${params.index_path}
	index_basename: ${params.index_basename}
	reads	: ${params.reads}
	"""
	.stripIndent(true)

process FASTQ {
	cpus 4

	input:
	path reads
	path index_path
	val index_basename
	

	output:
	path 'output.bam'

	script:
	"""
	hisat2 -p $task.cpus -x ${index_path}/${index_basename} -U $reads -S output.bam
	"""
}

workflow {
	index_ch = FASTQ(params.reads, 
					params.index_path, 
					params.index_basename)
}
