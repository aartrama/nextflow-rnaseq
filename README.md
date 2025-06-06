# Nextflow Pipeline 

This repository contains a Nextflow pipeline for processing bulk RNA-seq, ChIP-seq and Cut&Run datasets. 

RNA-seq workflow performs quality control, read trimming, alignment, feature counting, and summarization using MultiQC.

ChIP-seq and Cut&Run workflows perform quality control, read trimming, alignment, generation of sorted bam, bigwig, tdf, indexed bam files and summarization using MultiQC.


## Dependency

[Nextflow](https://www.nextflow.io/)
[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)

## Features

Support for single-end and paired-end RNA-Seq data
Singularity containerized execution.
Flexible customization for local or LSF cluster-based execution

## One Time Setup

Start an interactive session using the following command -

```bash
bsub -P acc_Nestlerlab -q interactive -R span[hosts=1] -n 8 -W 00:30 -Ip /bin/bash
```

Copy the setup.sh file from this repository to your project folder and set the singularity temp and cache folder in it. You do not need to modify anything else. Run the script using -

```bash
sh setup.sh
```

Setup should take ~ 10-15 minutes.


## Configuration Parameters

Copy the nextflow.config file from this repository and place it in your project folder. Your project folder should contain a directory called 'fastq' in which all fastq files are placed. Following should be the project folder structure :
```
.
├── fastq
│   ├── gut2_R1_001.fastq.gz
│   ├── gut2_R2_001.fastq.gz
│   ├── gut3_1.fastq.gz
│   ├── gut3_2.fastq.gz
│   ├── gut_R1.fastq.gz
│   └── gut_R2.fastq.gz
└── nextflow.config

```

Update the following parameters in your Nextflow config file:

| Parameter           | Description                                                                                  | Default Value |
|---------------------|----------------------------------------------------------------------------------------------|---------------|
| `project_dir`       | Base directory for the pipeline output.                                                     | '/path/to/project' |
| `index_files`    | Path to the HISAT2 index files.                                                          | '/path/to/index' |
| `index_basename`    | The index basename                                                          | 'mm39' |
| `experiment_type`    | Choose one among chipseq, rnaseq or cutrun                                                         | 'chipseq' |
| `multiqc_config`    | Path to multiqc config file. You can find one in this repo                                                         | "/path/to/multiqc_config.yaml" |

| `umi_present`          | Whether or not the library is UMI based                                                            | false |
| `umi_1`          | UMI sequence extracted from Read 1                                                             | "X" |
| `umi_2`          | UMI sequence extracted from Read 2                                                           | "NNNNNNNN" |
| `gtf_file`          | Path to the GTF annotation file.                                                            | "/path/to/gtf/mouse_genome_mm39.gtf" |
| `pairedEnd`         | Whether the data is paired-end (`true` or `false`).                                         | true |
| `count_unique`      | Count only uniquely mapped reads.                                                           | true |
| `count_fraction`    | Count fractional values for multi-mapped reads.                                             | false |
| `polyA`    | Whether the library is polyA based                                             | false |
| `singularity cacheDir`    | Set it to the same one you used in setup.sh                                            | "/path/to/singularity/cachedir" |



## Usage

**Run the pipeline**  
Use the following command to execute the pipeline:

For interactive node execution:

```bash
module load singularity/3.11.0
module load nextflow
module load git
nextflow run https://github.com/aartrama/nextflow-rnaseq -r main -profile local
```

For cluster execution:

```bash
module load singularity/3.11.0
module load nextflow
module load git
nextflow run https://github.com/aartrama/nextflow-rnaseq -r main -profile minerva
```


## References

[Nextflow Documentation](https://www.nextflow.io/docs/latest/index.html)

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

[HISAT2](http://daehwankimlab.github.io/hisat2/manual/)

[FeatureCounts](https://subread.sourceforge.net/featureCounts.html)

[Samtools](http://www.htslib.org/)

[UMI-tools](https://umi-tools.readthedocs.io/en/latest/)

[BamUtil](https://genome.sph.umich.edu/wiki/BamUtil)

[Bedtools](https://bedtools.readthedocs.io/en/latest/)

[deepTools](https://deeptools.readthedocs.io/en/develop/)

[Picard](https://broadinstitute.github.io/picard/)

[Pandas](https://pandas.pydata.org/docs/)

[MultiQC](https://multiqc.info/)




