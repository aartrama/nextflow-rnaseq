# Nextflow RNA-Seq Pipeline 

This repository contains a Nextflow pipeline for processing bulk RNA-seq data. The workflow performs quality control, read trimming, alignment, feature counting, and summarization using MultiQC.

## Dependency

- [Nextflow](https://www.nextflow.io/)
- [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)

## Features

- Support for single-end and paired-end RNA-Seq data
- Singularity containerized execution.
- Flexible customization for local or LSF cluster-based execution


### Configuration Parameters

Update the following parameters in your Nextflow config file:

| Parameter           | Description                                                                                  | Default Value |
|---------------------|----------------------------------------------------------------------------------------------|---------------|
| `project_dir`       | Base directory for the pipeline output.                                                     | `/path/to/project` |
| `input_dir`         | Input directory containing FASTQ files.                                                     | `/path/to/input` |
| `index_basename`    | Path to the HISAT2 index basename.                                                          | `/path/to/index` |
| `gtf_file`          | Path to the GTF annotation file.                                                            | `/path/to/gtf` |
| `pairedEnd`         | Whether the data is paired-end (`true` or `false`).                                         | `false` |
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

For local execution, ensure nextflow and singularity is installed on your computer.

```bash
nextflow run https://github.com/shenlab-sinai/NGS-Data-Charmer -r main -profile local
```

For cluster execution:

```bash
module load singularity/3.11.0
module load nextflow
nextflow run https://github.com/shenlab-sinai/NGS-Data-Charmer -r main -profile minerva
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


