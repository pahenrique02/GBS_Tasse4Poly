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

#Prepping the reference genome
echo 'Indexing the reference genome'
# If you have a fasta file with multiple sequences, you need to create a .fai

bwa index -a bwtsw SofficinarumxspontaneumR570_771_v2.0.hardmasked.fa

bwa aln -t 6 reference/R570Header_Simplified.fa mergedTagCounts/myMasterGBSTags.fq > mergedTagCounts/myMasterGBSTagsR570.sai

bwa samse reference/CollapsedMethyl.fa mergedTagCounts/myMasterGBSTagsAHVLYCDRX3.sai mergedTagCounts/myMasterGBSTagsAHVLYCDRX3.fq > mergedTagCounts/myMasterGBSTagsAHVLYCDRX3.sam 

samtools flagstat mergedTagCounts/myMasterGBSTags.sam

run_pipeline.pl -Xmx36g -fork1 -SAMConverterPlugin -i mergedTagCounts/alignment.sam -o topm/alignment.topm -endPlugin -runfork1 | tee VCF_NEW_CallingStep05.log.txt

run_pipeline.pl  -Xmx36g -fork1 -FastqToTBTPlugin -i fastq -k ./barcodes/Key161.txt -e PstI-MseI -o tbt/NEW -sh -m topm/alignment.topm -endPlugin -runfork1 | tee VCF_NEW_CallingSte06.log.txt

run_pipeline.pl  -Xmx20g -fork1 -MergeTagsByTaxaFilesPlugin  -i tbt/R570-o mergedTBT//alignment.tbt.shrt -x -endPlugin -runfork1 | tee VCF_NEW_CallingSte07.log.txt


run_pipeline.pl -Xmx48g -fork1 -DiscoverySNPCallerPlugin -i mergedTBT/AHVLYCDRX3+HFLKYBGX9myStudy.tbt.shrt -sh -m topm/AHVLYCDRX3+HFLKYBGX9myMasterGBSTags.topm -mUpd topm/AHVLYCDRX3+HFLKYBGX9myMasterTagsWithVariants.topm -o hapmap/raw/AHVLYCDRX3+HFLKYBGX9myGBSGenos_chr+.hmp.txt -vcf    -inclGaps  -ref reference/CollapsedMethyl.fa  -sC 1 -eC 1 -endPlugin -runfork1 | tee VCFCallingSte08.log.txt

run_pipeline.pl -Xmx32g -fork1 -MergeDuplicateSNPsPlugin -hmp hapmap/raw/AHVLYCDRX31HFLKYBGX9myGBSGenos_chr1.hmp.txt -o hapmap/mergedSNPs/myGBSGenos_mergedSNPs_chr+.hmp.txt   -misMat 0.3 -callHets  -sC 1 -eC 1 -endPlugin -runfork1  | tee VCFCallingSte09.log.txt

run_pipeline.pl -Xmx32g -fork1 -MergeIdenticalTaxaPlugin -hmp hapmap/mergedSNPs/AHVLYCDRX3myGBSGenos_mergedSNPs_chr+.hmp.txt -o hapmap/mergedTaxa/AHVLYCDRX3myGenos_taxaMerged_chr+.hmp.txt.gz -hetFreq 0.8 -sC 1 -eC 1 -endPlugin -runfork1  | tee VCFCallingSte10.log.txt

run_pipeline.pl -Xmx5g -fork1 -h AHVLYCDRX3+HFLKYBGX9myGBSGenos_chr+.hmp.txt -export -exportType VCF -runfork1

run_pipeline.pl -Xmx6g -fork1 -GBSHapMapFiltersPlugin -hmp ./AHVLYCDRX3myGBSGenos_chr1.hmp.txt -hLD -mnTCov 0.05 -mnSCov 0.05 -mnMAF 0.05 -o  DepthFiltered2.hmp.txt -sC 1 -eC 1 -endPlugin -runfork1