include { FASTQC } from '../modules/common/fastqc'
include { READ_TRIM_UMI } from '../modules/common/umi_tools'
include { READ_TRIM_GALORE_BASIC } from '../modules/common/trim_galore'
include { ALIGNMENT_STEP } from '../modules/common/alignment'
include { SORT_DISCORDANT_CUTRUN } from '../modules/cutrun/sort_discordant'
include { CLIP_DISCORDANT_SAMTOOL_CUTRUN } from '../modules/cutrun/clip_discordant'
include { CLIP_DISCORDANT_BAMUTIL_CUTRUN } from '../modules/cutrun/clip_discordant'
include { MERGE_DISCORDANT_CUTRUN } from '../modules/cutrun/merge_discordant'
include { UMI_RMDUPS_BAM } from '../modules/common/umi_tools'
include { PICARD_RMDUPS_BAM } from '../modules/common/bam_processing'
include { RMDUP_TO_CHRBAM } from '../modules/common/bam_processing'
include { CREATE_BW } from '../modules/common/bigwig'
include { CREATE_BED } from '../modules/common/bed'
include { CALC_INSERT_SIZE } from '../modules/common/bam_processing'
include { RUN_MULTIQC } from '../modules/common/multiqc'

workflow CUTRUN {
    take:
    read_pairs_ch

    main:
    fastqc_res=FASTQC(read_pairs_ch)
    read_trimming_ch = params.umi_present ? READ_TRIM_UMI(read_pairs_ch) : READ_TRIM_GALORE_BASIC(read_pairs_ch)
    sam_files_ch = ALIGNMENT_STEP(read_trimming_ch.trimmed_reads)
    cutrun_sorting=SORT_DISCORDANT_CUTRUN(sam_files_ch.sam)

    cutrun_clip_samtools=CLIP_DISCORDANT_SAMTOOL_CUTRUN(cutrun_sorting.mapped_dis_bam)
    cutrun_clip_bamutil=CLIP_DISCORDANT_BAMUTIL_CUTRUN(cutrun_clip_samtools.mapped_dis_bam_sorted)
    cutrun_merge=MERGE_DISCORDANT_CUTRUN(cutrun_sorting.mapped_paired_bam, cutrun_clip_bamutil.bam)

    rmdup_bam = params.umi_present ? UMI_RMDUPS_BAM(cutrun_merge.bam) : PICARD_RMDUPS_BAM(cutrun_merge.bam)
    chrbam=RMDUP_TO_CHRBAM(rmdup_bam.bam)
    bigwig_files_ch = CREATE_BW(chrbam.bam)
    bed_files_ch = CREATE_BED(chrbam.bam)
    if (params.pairedEnd) insert_size_ch = CALC_INSERT_SIZE(rmdup_bam.bam)
    multiqc_ch = RUN_MULTIQC(bed_files_ch.bed)

    emit:
    multiqc_report = multiqc_ch.report
}
