#!/bin/bash

# Author: PHSS
# Date: 04/09/2025

echo "I am at:" $PWD
echo "First, cloning Tassel4Poly repository"

# git clone https://github.com/guilherme-pereira/tassel4-poly.git

# Add to the PATH to avoid need to use absolute path
export PATH=$PATH:$HOME/Documentos/DataScience/BackupGBS_ServerUNESP/tassel4-poly

# Creating the required directories for the analysis
mkdir -p topm tbt tagCounts reference mergedTagCounts mergedTBT logs hapmap/raw fastq barcodes hapmap/mergedSNPs

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

#### IN THIS STEP USE YOUR REFERENCE GENOME  DIRECTLY IN THE "reference" directory, THE SYMBOLIC LINK IS ONLY FOR SCRIPT TESTING DUE REFERENCE GENOME SIZE NOT BE ABLE TO BE PUSHED TO GITHUB ####

## If use the methyl filtered reference genome from Grativol et al 2014, fix the header sed -i 's/>Chr/>chr/g' CollapsedMethyl.fa
## In the further stepes Samtools will recognize the chromosomes only if the header start with "chr" folloed by numbers.
#bwa index -a bwtsw ./reference/Reference

#echo 'Aligning the tags to the reference genome'

#bwa aln -t 6 reference/Reference mergedTagCounts/myMasterGBSTagsAHVLYCDRX3.fq > mergedTagCounts/myMasterGBSTagsR570.sai

#bwa samse reference/Reference mergedTagCounts/myMasterGBSTagsR570.sai mergedTagCounts/myMasterGBSTagsAHVLYCDRX3.fq > mergedTagCounts/myMasterGBSTagsAHVLYCDRX3.sam 

#echo "Alignment statistics"

#samtools flagstat mergedTagCounts/myMasterGBSTagsAHVLYCDRX3.sam

#run_pipeline.pl -Xmx12g -fork1 -SAMConverterPlugin -i mergedTagCounts/myMasterGBSTagsAHVLYCDRX3.sam -o topm/alignment.topm -endPlugin -runfork1 | tee ./logs/VCF_NEW_CallingStep05.log.txt

#run_pipeline.pl  -Xmx16g -fork1 -FastqToTBTPlugin -i fastq -k ./barcodes/Key161.txt -e PstI-MseI -o tbt -sh -m topm/alignment.topm -endPlugin -runfork1 | tee ./logs/VCF_NEW_CallingSte06.log.txt

#run_pipeline.pl  -Xmx20g -fork1 -MergeTagsByTaxaFilesPlugin  -i tbt/ -o mergedTBT/alignment.tbt.shrt -x -endPlugin -runfork1 | tee ./logs/VCF_NEW_CallingSte07.log.txt

echo 'SNP Calling'

#run_pipeline.pl -Xmx20g -fork1 -DiscoverySNPCallerPlugin -i mergedTBT/alignment.tbt.shrt -sh -m topm/alignment.topm -mUpd topm/alignment.topm -o hapmap/raw/myGBSGenos_chr+.hmp.txt -vcf    -inclGaps  -ref reference/Reference  -sC 1 -eC 1 -endPlugin -runfork1 | tee ./logs/VCFCallingSte08.log.txt

#run_pipeline.pl -Xmx20g -fork1 -MergeDuplicateSNPsPlugin -hmp hapmap/raw/myGBSGenos_chr1.hmp.txt -o hapmap/mergedSNPs/myGBSGenos_mergedSNPs_chr+.hmp.txt   -misMat 0.3 -callHets  -sC 1 -eC 1 -endPlugin -runfork1  | tee ./logs/VCFCallingSte09.log.txt
#run_pipeline.pl -Xmx20g		-fork1 -MergeDuplicateSNPsPlugin -vcf hapmap/raw/myGBSGenos_chr1.vcf -o hapmap/mergedSNPs/myGBSGenos_mergedSNPs_chr+.vcf -misMat 0.3  -callHets -sC 1 -eC 1 -endPlugin -runfork1 
#run_pipeline.pl -Xmx20g -fork1 -MergeIdenticalTaxaPlugin -hmp hapmap/mergedSNPs/myGBSGenos_mergedSNPs_chr+.hmp.txt -o hapmap/mergedTaxa/myGenos_taxaMerged_chr+.hmp.txt -hetFreq 0.8 -sC 1 -eC 1 -endPlugin -runfork1  | tee ./logs/VCFCallingSte10.log.txt

echo 'Converting to VCF'

run_pipeline.pl -Xmx15g -fork1 -h hapmap/mergedTaxa/myGenos_taxaMerged_chr+.hmp.txt -export -exportType VCF -runfork1

echo 'Pipeline finished'

echo 'The final VCF files will be filtered with VCFTOOLS'

vcftools \
  --vcf YOUR_RAW_VCF.vcf \
  --maf 0.05 \
  --max-missing 0.75 \
  --min-alleles 2 \
  --max-alleles 2 \
  --minDP 50 \
  --out YOUR_Filtered \
  --recode --recode-INFO-all
