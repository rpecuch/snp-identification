# Single Nucleotide Polymorphism Identification
# Author: Rita Pecuch 

The goal of this pipeline is to identify single nucleotide polymorphisms (SNPs) present in the input DNA samples.

Input files: 
- hawkins_pooled_sequences.fastq: all samples to analyze in FASTQ format
- harrington_clinical_data.txt: tab-delimited text file containing information about each DNA sample
- dgorgon_reference.fa: FASTA file containing reference genome
- utils.py: Python file containing several utility functions needed by the pipeline

Output files: 
- fastqs/: contains FASTQ file for each DNA sample with barcode and regions of low quality scores trimmed
- bamfiles/: contains sorted BAM file for each DNA sample
- report.txt: contains the name, color, and number of reads for each DNA sample as well as the percentage of reads in the sample containing the identified SNP

Execution:

1) Copy over the following files into the directory that you would like to run the script:
- hawkins_pooled_sequences.fastq
- harrington_clinical_data.txt
- dgorgon_reference.fa
- utils.py
- pipeline.py

2) Execute the following command: python3 pipeline.py --fastq hawkins_pooled_sequences.fastq --txt harrington_clinical_data.txt --fa dgorgon_reference.fa --utils utils.py