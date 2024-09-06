# Class to read data from fastq file
class ParseFastQ(object):
    # Initialize instance of class
    def __init__(self,filePath,headerSymbols=['@','+']):
        # Assign properties
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            # Open file in read mode
            # U: universal newlines
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
         
    # Initalize interator object
    def __iter__(self):
        return self
     
    # Define what is returned on next iteration
    def __next__(self):
        # Get next 4 lines and place in list
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
         
        # Count number of true lines and number of None
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # Check for acceptable end of file
        if nones == 4:
            raise StopIteration
        # Make sure we got 4 full lines of data
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # Make sure we are in the correct "register"
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # Make sure the seq line and qual line have equal lengths
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
         
        # Return fastQ data as tuple
        # (header, sequence, separator, quality score)
        return tuple(elemList)

# Function to read clinical data file into dictionary
def read_clinical_data(txt_file):
    with open(txt_file) as f:
        reader = csv.reader(f, delimiter="\t")
        # Skip header row
        next(reader, None)
        clinical_list = list(reader)
    clinical_dict = {}
    for entry in clinical_list:
        # Use barcodes as keys in dictionary
        clinical_dict[entry[2]] = {'name': entry[0], 'color': entry[1], 'sequences': []}
    return clinical_dict

# Function to trim barcode from sequence and quality score
def trim_barcode(sequence, qual_score):
    # Read barcode
    barcode = sequence[0:5]
    # Remove from sequence and quality score
    trimmed_seq = sequence[5:len(sequence)]
    trimmed_score = qual_score[5:len(sequence)]
    
    return barcode, trimmed_seq, trimmed_score

# Function to trim sequence and quality score when low quality score
def trim_score(qual_score, sequence):
    options = ['D','F']
    # Initialize index and previous letter
    i=0
    prev_letter=None
    # Loop through quality score
    for letter in qual_score:
        if prev_letter:
            # If consecutive D or F
            if prev_letter in options and letter in options:
                # Identify where to trim string
                split = i-1
                # Trim quality score and sequence
                qual_score = qual_score[0:split]
                sequence = sequence[0:split]
                # Stop after trim
                break
        # Increment
        prev_letter = letter
        i+=1

    return qual_score, sequence

# Function to get information about an entry in fastq file
def get_entry_info(fastq_obj):
    header = fastq_obj[0]
    sequence = fastq_obj[1]
    sep = fastq_obj[2]
    qual_score=fastq_obj[3]

    # Trim off barcode
    barcode, trimmed_seq, trimmed_score = trim_barcode(sequence, qual_score)
    # Trim quality score
    trimmed_score, trimmed_seq = trim_score(trimmed_score, trimmed_seq)

    return barcode, [header, trimmed_seq, sep,trimmed_score]

# Function to create new directory if needed
def make_dir(dir_name):
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

# Function to write trimmed fastq files
def write_trimmed_fastq(entry, file_path):
    content_to_write = []
    for sequence in entry['sequences']:
        content_to_write.append('\n'.join(sequence))
        # Each entry should be on new line
        content_to_write.append('\n')

    with open(file_path, 'w') as f:
        f.writelines(content_to_write)

# Function to perform aligment with bwa
def bwa_align(fastq, fa, fastqs_dir, samfiles_dir):
    name = fastq.split('_')[0]
    os.system(f'bwa mem {fa} {fastqs_dir}/{fastq} > {samfiles_dir}/{name}.sam')

# Function to convert SAM files to BAM files, then sort and index
def sam_to_bam(samfile, samfiles_dir, bamfiles_dir):
    name = samfile.split('.')[0]
    os.system(f'samtools view -bS {samfiles_dir}/{samfile} > {bamfiles_dir}/{name}.bam')
    # Sort
    os.system(f'samtools sort -m 100M -o {bamfiles_dir}/{name}.sorted.bam {bamfiles_dir}/{name}.bam')
    # Index
    os.system(f'samtools index {bamfiles_dir}/{name}.sorted.bam')
    # Remove non-sorted file
    os.system(f'rm {bamfiles_dir}/{name}.bam')

# Discover variants
def pileup(sorted_bam, bamfiles_dir):
    # Read BAM file
    samfile = pysam.AlignmentFile(f'{bamfiles_dir}/{sorted_bam}', "rb")
    
    # Use a dictionary to keep mutations for entire sample
    sample_dict = {}

    # Iterate over all positions
    for pileupcolumn in samfile.pileup():
        # You can uncomment the below line to see the coverage for each base
        # print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))

        # Use a dictionary to count up the bases for specific position
        ntdict = {}
        # Iterate over reads that are part of sequence coverage
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # You can uncomment the below line to see what is happening in the pileup
                # print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))

                # Retrieve base present at specific position
                base = pileupread.alignment.query_sequence[pileupread.query_position]

                # Retrieve base present at same position in reference sequence
                if 'ref_base' not in ntdict.keys():
                    ref_seq = pileupread.alignment.get_reference_sequence()
                    ref_base = ref_seq[pileupread.query_position]
                    ntdict['ref_base'] = ref_base.upper()

                # Populate the ntdict with the counts of each base
                # This dictionary will hold all of the base read counts per nucletoide per position
                if base in ntdict.keys():
                    ntdict[base] += 1
                else:
                    ntdict[base] = 1

        # Add to dictionary of mutations for sample if frequency is not 100% or 0%
        if len(ntdict) > 2:
            sample_dict[pileupcolumn.pos] = ntdict

    samfile.close()
    return sample_dict

# Generate file content
def generate_file_content(sorted_bam, clinical_df, sample_dict):
    # Write sentence one for report
    sample_name = sorted_bam.split('.')[0]
    sample_df = clinical_df.loc[clinical_df['name'] == sample_name].reset_index()
    color = sample_df.loc[0, 'color']
    num_reads = len(sample_df.loc[0, 'sequences'])

    sentences = [f'Sample {sample_name} had a {color} mold with {num_reads} reads. ']

    # Write subsequent sentence for report
    if sample_dict:
        read_pos = list(sample_dict.keys())[0]
        ref_base = sample_dict[read_pos]['ref_base']
        for key in sample_dict[read_pos]:
            if key != 'ref_base' and key != ref_base:
                num_mutations = sample_dict[read_pos][key]
                mutation = key
                # Calculate percent of reads with mutation
                pct_mutation = (num_mutations / num_reads) * 100
                rounded_pct = round(pct_mutation, 2)

                text_string = f'{rounded_pct} % of the reads at position {read_pos} had the mutation {mutation}. '
                sentences.append(text_string)

    # Return all sentences with spacing at end for readability
    return f'{"".join(sentences)}\n\n'

# Write file content
def append_file(file_content, file_name):
    with open(file_name, 'a') as f:
        f.write(file_content)