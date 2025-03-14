include { FASTQC } from '../modules/common/fastqc'
include { READ_TRIM_UMI } from '../modules/common/umi_tools'
include { READ_TRIM_GALORE_BASIC } from '../modules/common/trim_galore'
include { READ_TRIM_GALORE_POLYA } from '../modules/common/trim_galore'
include { ALIGNMENT_STEP } from '../modules/common/alignment'
include { SORTED_BAM_STEP } from '../modules/common/bam_processing'
include { COUNTS_STEP } from '../modules/rnaseq/counts'
include { COUNTS_MATRIX } from '../modules/rnaseq/counts_matrix'
include { RUN_MULTIQC } from '../modules/common/multiqc'

workflow RNASEQ {
    take:
    read_pairs_ch

    main:
    fastqc_files_ch = FASTQC(read_pairs_ch)
    read_trimming_ch = params.umi_present ? READ_TRIM_UMI(read_pairs_ch) : READ_TRIM_GALORE_BASIC(read_pairs_ch)
    if (params.polyA) read_trimming_ch = READ_TRIM_GALORE_POLYA(read_trimming_ch.trimmed_reads)
    sam_files_ch = ALIGNMENT_STEP(read_trimming_ch.trimmed_reads)
    sorted_bam_files_ch = SORTED_BAM_STEP(sam_files_ch.sam)
    count_files_ch = COUNTS_STEP(sorted_bam_files_ch.bam)
    counts_matrix = COUNTS_MATRIX(count_files_ch.counts_exonic)
    multiqc_report_file = RUN_MULTIQC(counts_matrix.count_matrix_genic)

    emit:
    multiqc_report = multiqc_report_file.report

}
