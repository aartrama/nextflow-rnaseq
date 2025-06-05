SINGULARITY_CACHE_DIR="/sc/arion/projects/Nestlerlab/shenl03_rna/aarthi/2025/li/test_nextflow_jun4/singularity_cache"
SINGULARITY_TMP_DIR="/sc/arion/projects/Nestlerlab/shenl03_rna/aarthi/2025/li/test_nextflow_jun4/tmp"
mkdir -p $SINGULARITY_CACHE_DIR
mkdir -p $SINGULARITY_TMP_DIR  

module load singularity/3.11.0
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --name rnaseq.sif docker://nfcore/rnaseq:1.4.2
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --name fastqc.sif docker://staphb/fastqc
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --name trim-galore.sif docker://dceoy/trim_galore:latest
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --name hisat2.sif docker://aarthir239/hisat2
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --name samtools.sif docker://biocontainers/samtools:v1.9-4-deb_cv1
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --name featurecounts.sif docker://thatdnaguy/featurecounts:v2.0.6_02
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --name umitools.sif docker://getwilds/umitools:latest
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --name bamutil.sif docker://marrip/bamutil:1.0.15
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --name bedtools.sif docker://staphb/bedtools:latest
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --name deeptools.sif docker://mgibio/deeptools:3.5.3
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --name picard.sif docker://mgibio/picard:v3.4.0-noble
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --tmpdir $SINGULARITY_TMP_DIR --name pandas.sif docker://biocontainers/pandas:1.5.1_cv1
singularity pull --dir $SINGULARITY_CACHE_DIR --disable-cache --tmpdir $SINGULARITY_TMP_DIR --name multiqc.sif docker://multiqc/multiqc:latest
