#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { CUTRUN } from './workflows/cutrun'
include { CHIPSEQ } from './workflows/chipseq'
include { RNASEQ } from './workflows/rnaseq'

// Set parameters
params.experiment_type = 'rnaseq'
params.index_files = 'path/to/index/files'
params.index_basename = 'index_basename'
params.project_dir = 'path/to/project'
params.pairedEnd = true
params.gtf_file = 'path/to/gtf/file'
params.multiqc_config = 'path/to/multiqc/config'
params.count_unique = true
params.count_fraction = false
params.umi_present = false
params.polyA = false

// Define the input channel 
if (params.pairedEnd) {
    read_pairs_ch = Channel.fromFilePairs("${params.project_dir}/fastq/*_{R1,R2,1,2},{,_001,_S[0-9]+},{,_001},{,.fastq,.fq},{,.gz}", checkIfExists: true)
} else {
    read_pairs_ch = Channel.fromPath("${params.project_dir}/fastq/*.{fastq,fq},{,.gz}", checkIfExists: true)
        .map { file -> tuple(file.simpleName, file) }
}

workflow {
    switch(params.experiment_type) {
        case 'cutrun':
            CUTRUN(read_pairs_ch)
            break
        case 'chipseq':
            CHIPSEQ(read_pairs_ch)
            break
        case 'rnaseq':
            RNASEQ(read_pairs_ch)
            break
        default:
            error "Invalid experiment type: ${params.experiment_type}"
    }
}
