## Run script: python3 pipeline.py --fastq hawkins_pooled_sequences.fastq --txt harrington_clinical_data.txt --fa dgorgon_reference.fa --utils utils.py

## Import needed modules
import argparse
import csv
import os
import pandas as pd
import pysam

## Run pipeline
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # Define required command line args
    parser.add_argument("-f", "--fastq", required=True, help="Relative path to fastq file containing all sequences to analyze")
    parser.add_argument("--txt", required=True, help="Relative path to clinical data txt file")
    parser.add_argument("--fa", required=True, help="Relative path to fa file containing reference genome")
    parser.add_argument("--utils", required=True, help="Relative path to utility function file")
    # Parse command line args
    args = parser.parse_args()

    ## Import needed classes and functions
    exec(open(args.utils).read())

    ## Read clinical data file into dictionary
    clinical_dict = read_clinical_data(args.txt)

    ## Read fastq file and perform trimming

    # Create ParseFastQ object
    fastqfile = ParseFastQ(args.fastq)

    # Loop through entries in fastq file
    for fastq_obj in fastqfile:
        # Get entry information
        barcode, entry = get_entry_info(fastq_obj)

        # Use barcode to match sequence to patient
        clinical_dict[barcode]['sequences'].append(entry)

    # Create output directory
    make_dir('fastqs')

    # Write new fastq files
    for key in clinical_dict:
        entry = clinical_dict[key]
        file_path = f'fastqs/{entry["name"]}_trimmed.fastq'
        write_trimmed_fastq(entry, file_path)

    ## Perform alignment

    # Index reference file
    os.system(f'bwa index {args.fa}')
    # Create output dir
    make_dir('samfiles')
    # Retrieve input files
    input_fastqs = os.listdir('fastqs')
    # Perform alignment
    for fastq in input_fastqs:
        bwa_align(fastq, args.fa, 'fastqs', 'samfiles')

    ## Convert samfiles to bamfiles

    # Create output dir
    make_dir('bamfiles')
    # Retrieve input files
    input_sams = os.listdir('samfiles')
    # Convert
    for samfile in input_sams:
        sam_to_bam(samfile, 'samfiles', 'bamfiles')
    # Remove old sam files
    os.system('rm -r samfiles/')

    ## Variant discovery and generate report

    # Retrieve input files
    input_bams = os.listdir('bamfiles')

    # Convert dictionary to dataframe to easily retrieve needed info for report
    clinical_df = pd.DataFrame.from_dict(clinical_dict, orient='index')
    
    for sorted_bam in input_bams:
        # Do not run on indexed files, as indexes will be opened automatically when sorted BAM files are opened
        if sorted_bam.endswith('bam'):
            # Run variant discovery function on each sorted BAM
            sample_dict = pileup(sorted_bam, 'bamfiles')

            # Generate report content for sample
            file_content = generate_file_content(sorted_bam, clinical_df, sample_dict)

            # Add to report
            append_file(file_content, 'report.txt')
