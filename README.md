RNA-Seq Transcriptomics Pipeline

This repository provides a Linux-based RNA-Seq transcriptomics analysis pipeline that performs end-to-end processing of RNA sequencing data, starting from raw FASTQ files to a merged gene expression count matrix. The pipeline is designed to be reproducible, interactive, and easy to understand for students and researchers.

1️⃣ What is RNA-Seq?

RNA sequencing (RNA-Seq) is a high-throughput sequencing technique used to study the transcriptome, i.e., all RNA molecules expressed in a biological sample at a given time.

RNA-Seq allows us to:

Measure gene expression levels

Identify differentially expressed genes

Discover novel transcripts and splice variants

Study functional genomics and regulatory mechanisms

Typical RNA-Seq workflow:

RNA isolation

Library preparation

Sequencing (Illumina, etc.)

Quality control

Read trimming

Alignment to a reference genome

Gene-level read counting

Downstream statistical analysis (DESeq2, edgeR, etc.)

This repository covers steps 4–7.

2️⃣ Getting Raw RNA-Seq Data from NCBI

Public RNA-Seq datasets can be downloaded from NCBI Sequence Read Archive (SRA).

Example workflow:

Identify a project (e.g., PRJNA394193)

Install SRA Toolkit

Download and convert data to FASTQ format

fasterq-dump SRRXXXXXXX --split-files --threads 8


This pipeline assumes that each sample has its own folder containing paired-end FASTQ files:

PROJECT_DIR/
├── 001/
│   ├── sample_1.fastq.gz
│   └── sample_2.fastq.gz
├── 002/
│   ├── sample_1.fastq.gz
│   └── sample_2.fastq.gz

3️⃣ install_tools.sh – Tool Installation & Setup

This script prepares the system and pipeline environment.

What it does:

Updates system package list

Installs core dependencies (Java, Python, wget, unzip, etc.)

Installs bioinformatics tools:

HISAT2 (read alignment)

SAMtools (BAM/SAM processing)

FastQC (quality control)

MultiQC (QC summarization)

Downloads and installs Trimmomatic (v0.39) locally inside the pipeline folder

Creates required directories:

genome_data/

hisat2_index/

Why this script is important:

Ensures reproducibility

Avoids manual dependency errors

Safe to run multiple times

Run once:

./install_tools.sh

4️⃣ fastqc_trim.sh – Quality Control & Trimming

This script performs quality control and optional trimming on raw sequencing reads.

Step-by-step actions:

Automatically detects paired-end FASTQ files

Runs FastQC on raw reads

Pauses and asks the user whether trimming should be performed

If approved:

Runs Trimmomatic for adapter and quality trimming

Runs FastQC again on trimmed reads

Key features:

Interactive (user-controlled trimming)

Automatic sample detection

Per-sample output organization

Prevents reprocessing of completed samples

Outputs are stored inside each sample folder:

sample/
├── fastqc_raw/
├── trimmed/
│   ├── *_paired.fastq.gz
│   └── *_unpaired.fastq.gz

5️⃣ run_alignment.sh – Alignment & Read Counting

This script aligns trimmed reads to the reference genome and generates gene-level counts.

Major steps:

Automatically detects:

Genome FASTA

GTF annotation

Builds HISAT2 index (if not already present)

Aligns trimmed reads using HISAT2

Converts SAM → sorted BAM

Indexes BAM files

Performs HTSeq-count to generate gene expression counts

Output structure:
sample/
├── hisat2/
│   └── sample.sorted.bam
├── htseq_counts/
│   └── sample_counts.txt


This script ensures:

No redundant index rebuilding

No re-alignment if BAM already exists

Clean and organized outputs

6️⃣ merge.py – Merge Count Files

This Python script merges individual sample count files into a single expression matrix.

What it does:

Automatically detects sample folders

Reads HTSeq count files from each sample

Merges all counts by gene ID

Exports a single Excel file

Final output:
All_samples_counts.xlsx


This file can be directly used for:

DESeq2

edgeR

Limma

Other downstream transcriptomic analyses

7️⃣ Folder Descriptions
📁 genome_data/

Contains reference genome files:

Genome FASTA (.fa, .fasta, .fa.gz)

Annotation file (.gtf, .gtf.gz)

Example:

genome_data/
├── genome.fa
└── annotation.gtf

📁 hisat2_index/

Stores HISAT2 genome index files (*.ht2) generated automatically.

Purpose:

Speeds up alignment

Avoids rebuilding index every run

8️⃣ Recommended Pipeline Execution Order
./install_tools.sh
./fastqc_trim.sh
./run_alignment.sh
python3 merge.py
