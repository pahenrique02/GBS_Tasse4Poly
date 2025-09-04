#!/bin/bash

# Author: PHSS
# Date: 04/09/2025

echo "I am at:" $PWD
echo "First, cloning Tassel4Poly repository"

# git clone https://github.com/guilherme-pereira/tassel4-poly.git

# Add to the PATH to avoid need to use absolute path
export PATH=$PATH:$HOME/Documentos/DataScience/BackupGBS_ServerUNESP/tassel4-poly

# Creating the required directories for the analysis
mkdir -p topm tbt tagCounts reference mergedTagCounts mergedTBT logs hapmap fastq barcodes

# Place your reference genome at reference folder, your barcode at barcodes folder 
# and in the fastq folder with your fastq files using the flowcell name and lanes 
# with the extension _fastq.txt.gz.

echo 'Starting the pipeline'

run_pipeline.pl -Xmx12g -fork1 -FastqToTagCountPlugin \
  -i ./fastq \
  -k ./barcodes/Key161.txt \
  -e PstI-MseI -c 1 \
  -o tagCounts \
  -endPlugin -runfork1 | tee ./logs/VCFCallingstep01.log.txt

echo 'Merging TagCounts'
run_pipeline.pl  -Xmx12g -fork1 -TagCountToFastqPlugin \
-i ./tagCounts/AHVLYCDRX3_1.cnt  \
-o ./mergedTagCounts/myMasterGBSTagsAHVLYCDRX3.fq \
-c 1  -endPlugin -runfork1 | tee ./logs/VCFCallingStep02.log.txt