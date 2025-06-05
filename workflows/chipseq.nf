include { FASTQC } from '../modules/common/fastqc'
include { READ_TRIM_UMI } from '../modules/common/umi_tools'
include { READ_TRIM_GALORE_BASIC } from '../modules/common/trim_galore'
include { READ_TRIM_GALORE_POLYA } from '../modules/common/trim_galore'
include { ALIGNMENT_STEP } from '../modules/common/alignment'
include { SORT_BAM_CHIPSEQ } from '../modules/chipseq/sort_bam'
include { UMI_RMDUPS_BAM } from '../modules/common/umi_tools'
include { PICARD_RMDUPS_BAM } from '../modules/common/bam_processing'
include { RMDUP_TO_CHRBAM } from '../modules/common/bam_processing'
include { CREATE_BW } from '../modules/common/bigwig'
include { CREATE_BED } from '../modules/common/bed'
include { CALC_INSERT_SIZE } from '../modules/common/bam_processing'
include { RUN_MULTIQC } from '../modules/common/multiqc'

workflow CHIPSEQ {
    take:
    read_pairs_ch

    main:

    fastqc_res = FASTQC(read_pairs_ch)
    read_trimming_ch = params.umi_present ? READ_TRIM_UMI(read_pairs_ch) : READ_TRIM_GALORE_BASIC(read_pairs_ch)
    if (params.polyA) read_trimming_ch = READ_TRIM_GALORE_POLYA(read_trimming_ch.trimmed_reads)
    alignment_ch = ALIGNMENT_STEP(read_trimming_ch.trimmed_reads)
    sam_sorting = SORT_BAM_CHIPSEQ(alignment_ch.sam)
    rmdup_bam = params.umi_present ? UMI_RMDUPS_BAM(sam_sorting.bam) : PICARD_RMDUPS_BAM(sam_sorting.bam)
    chrbam = RMDUP_TO_CHRBAM(rmdup_bam.bam)
    bw_ch = CREATE_BW(chrbam.bam)
    bed_files_ch = CREATE_BED(chrbam.bam)
    if (params.pairedEnd) insert_size_ch = CALC_INSERT_SIZE(rmdup_bam.bam)
    multiqc_ch = RUN_MULTIQC(bed_files_ch)

    emit:
    multiqc_report = multiqc_ch.report

}
