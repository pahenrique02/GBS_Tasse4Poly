# GBS_Tasse4Poly

Automated pipeline for **Genotyping-by-Sequencing (GBS)** analysis in polyploids using **TASSEL-4-Poly**.

---

## Description

This pipeline, written in **Shell** and **Perl**, processes raw sequencing data (FASTQ) all the way to variant calling (VCF).  
It includes steps for demultiplexing, tag counting, building intermediate files (TBT, tagCounts, topm), and running TASSEL-4-Poly, producing reproducible results with detailed logs.

---

## Requirements

- Linux/Unix system  
- Java (>= 1.8)  
- Perl  
- TASSEL-4-Poly  
- FASTQ files of genotypes  
- Barcode file  
- Reference genome (FASTA)

---

## Installation / Setup

```bash
# Clone the repository
git clone https://github.com/pahenrique02/GBS_Tasse4Poly.git
cd GBS_Tasse4Poly

# Make the pipeline executable
chmod +x GBS_Pipeline.sh



