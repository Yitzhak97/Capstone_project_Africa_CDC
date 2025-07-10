# Capstone Project – eboVar: A Bash Pipeline for EBOV Variant Profiling and Analysis (Group 3)

## Description

This capstone project challenges participants to design and implement a comprehensive bioinformatics pipeline for analyzing raw Illumina sequencing data of Ebola virus isolates. The project emphasizes command-line scripting, data manipulation, variant calling, data wrangling and analysis, and containerization for reproducibility.

### Tools used in the pipeline
1. FastQC         
2. fastp          
3. MultiQC        
4. BWA            
5. Samtools
6. BCFtools


## Pipeline Steps

The pipeline (eboVar.sh) performs:

1. FastQC on raw reads

2. Trimming with fastp

3. FastQC on trimmed reads

4. Alignment with BWA

5. Index BAM

6. Variant calling: bcftools

7. Sort + index VCF


## Requirements

Linux environment

Conda

Apptainer (for container use)

Tools listed in containers/ebovar.yml


## How to Run

📥 Clone repository & prepare

git clone https://github.com/<your_username>/Capstone_EboVar.git
cd Capstone_EboVar

✅ Run locally with conda

conda env create -f containers/ebovar.yml
conda activate ebovar

bash scripts/eboVar.sh -i data/raw -o results -r data/reference/ebov_ref.fa -t 4


## Containerized workflow

Build Apptainer container

apptainer build containers/ebovar.sif containers/ebovar.def

Run pipeline inside container

apptainer run containers/ebovar.sif -i data/raw -o results -r data/reference/ebov_ref.fa -t 4


## Downstream analysis 📊

1️⃣ Filter high-quality variants:

bash scripts/test.sh

2️⃣ Run R script for summary and plots:

Rscript scripts/variant_analysis.R


## Outputs 📑

MultiQC HTML report (results/qc/)

Sorted BAM files & indices (results/bam/)

VCF files (results/vcf/)

Filtered TSV files per sample

Merged high-quality variant table (hq_allvariants.tsv)

R Markdown/PDF report with plots


## Documentation 📖

All scripts are well-commented.
README.md explains:

Purpose & structure

Usage instructions

Containerization details

Expected outputs

---

✏ Authors & Training

Created as part of the Africa CDC Bioinformatics Fellowship Program (ACDC) – Kigali, July 2025.
Participants: 
1. Makonk Najah: NMIMR, Ghana
2. Stella E. Nabirye: UVRI, Uganda
3. Isaac Adison: NPHL, South Sudan

---

✅ License

Open-source, educational use.

