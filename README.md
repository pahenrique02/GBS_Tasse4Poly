# GBS Pipeline

GBS_Tassel4Poly is an automated pipeline written in Shell that uses Perl scripts from the original Tassel Pipeline. It also specialises in plugins for polyploids, which were developed by Pereira et al. (2018) for analysing genotyping-by-sequencing (GBS) data using TASSEL-4-Poly. It processes small subset of raw FASTQ data from sugarcane, performs barcode demultiplexing and generates TASSEL-compatible intermediate files, such as TBT, tagCounts and topm. It also performs variant calling (VCF) with traceability through logs. This can be modified for polyploids according to your need and availablity (optional) of reference genomes.

 REFERENCES
Glaubitz JC, Casstevens TM, Lu F, Harriman J, Elshire RJ, Sun Q, Buckler ES. (2014) TASSEL-GBS: A High Capacity Genotyping by Sequencing Analysis Pipeline. PLoS ONE 9(2): e90346. https://doi.org/10.1371/journal.pone.0090346
Pereira GS, Garcia AAF, Margarido GRA. (2018) A fully automated pipeline for quantitative genotype calling from next generation sequencing data in autopolyploids. BMC Bioinformatics 19:398. https://doi.org/10.1186/s12859-018-2433-6.







