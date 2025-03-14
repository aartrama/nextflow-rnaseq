process FASTQC {
    publishDir "${project_dir}/output/fastqc", mode: 'copy'

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
