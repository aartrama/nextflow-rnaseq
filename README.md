# Nextflow RNA-Seq Pipeline 

This repository contains a Nextflow pipeline for processing bulk RNA-seq data. The workflow performs quality control, read trimming, alignment, feature counting, and summarization using MultiQC.

## Dependency

- [Nextflow](https://www.nextflow.io/)
- [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)

## Features

- Support for single-end and paired-end RNA-Seq data
- Singularity containerized execution.
- Flexible customization for local or LSF cluster-based execution

## One Time Setup

Start an interactive session using the following command -

```bash
bsub -P acc_Nestlerlab -q interactive -R span[hosts=1] -n 8 -W 00:30 -Ip /bin/bash
```

Copy the setup.sh file from this repository to your project folder and set the temp and singularity cache folder. Run the script using -

```bash
sh setup.sh
```

Setup should take about 10 minutes.


## Configuration Parameters

Copy the nextflow.config file from this repository and place it in your project folder. Your project folder should contain a directory called 'fastq' in which all fastq files are placed. Update the following parameters in your Nextflow config file:

| Parameter           | Description                                                                                  | Default Value |
|---------------------|----------------------------------------------------------------------------------------------|---------------|
| `project_dir`       | Base directory for the pipeline output.                                                     | `/path/to/project` |
| `index_files`    | Path to the HISAT2 index files.                                                          | `/path/to/index` |
| `index_basename`    | The index basename                                                          | `mm39` |
| `gtf_file`          | Path to the GTF annotation file.                                                            | `/path/to/gtf/mouse_genome_mm39.gtf` |
| `pairedEnd`         | Whether the data is paired-end (`true` or `false`).                                         | `true` |
| `multiqc_config`    | Path to MultiQC configuration file.                                                         | `/path/to/multiqc_config.yaml` |
| `count_unique`      | Count only uniquely mapped reads.                                                           | `true` |
| `count_fraction`    | Count fractional values for multi-mapped reads.                                             | `false` |

Place nextflow.config file in your project directory. Your project directory should also contain a folder called 'fastq' containing your fastq files. Following should be the project directory structure :
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

## Pipeline Steps

1. **FASTQC**  
   Performs quality control on raw sequencing reads using FastQC.

2. **Trim Galore**  
   Trims low-quality regions and adapter sequences from reads.

3. **HISAT2 Alignment**  
   Aligns trimmed reads to a reference genome using HISAT2.

4. **SAM to BAM Conversion**  
   Sorts and converts SAM files to BAM format using Samtools.

5. **Feature Counts**  
   Generates gene- and exon-level count matrices using featureCounts.

6. **Counts Matrix Generation**  
   Compiles counts from all samples into final matrix files.

7. **MultiQC Report**  
   Summarizes results from previous steps into an HTML report.

## References

[Nextflow Documentation](https://www.nextflow.io/docs/latest/index.html)

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[Trim Galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

[HISAT2](http://daehwankimlab.github.io/hisat2/manual/)

[FeatureCounts](https://subread.sourceforge.net/featureCounts.html)


