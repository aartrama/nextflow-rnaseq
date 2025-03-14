process RUN_MULTIQC {
    publishDir "${project_dir}/output/multiqc", mode: 'copy'

    input:
    path("output")

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_data"), emit: data
    path("multiqc_summary_text.txt")

    script:
    qc_config = multiqc_config
    """
    multiqc ${project_dir} --config ${qc_config} --force &> multiqc_summary_text.txt
    """
}
