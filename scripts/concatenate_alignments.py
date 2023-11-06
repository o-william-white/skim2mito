import argparse
import os
import pandas as pd
from functools import reduce
import sys
import shutil

usage = """

python concatenate_alignments.py --input example_data/alignments/ --type DNA --output example_output --prefix example 

"""

description = """
Python script to create a concatenated alignment of aligned gene sequences and partition file suitable for raxaml or iqtree. This script requires the Pandas library. Genes must be (1) aligned, (2) in fasta format, (3) in a single directory, (4) named <gene>.fasta and (5) with sequences named based on <sample> only. Pay special attention to the naming of the fasta files as <gene>.fasta and the sequences as <sample>, as this information is used to build the concatentated alignments. 
"""

# argparse
parser = argparse.ArgumentParser(usage=usage, description=description)
parser.add_argument("--input",     help = "Input directory containing fasta files.", required=True)
parser.add_argument("--type",      help = "Alignment type.", choices = ['DNA', 'PROTEIN'], required=True)
parser.add_argument("--prefix",    help = "Output file prefix.", required=True)
parser.add_argument("--output",    help = "Output directory.", required=True)
parser.add_argument("--overwrite", help = "Overwrite output directory.", required=False, action = "store_true")
args = parser.parse_args()

# read multi-line fasta 
def read_fasta(filename):
    name, seq = None,""
    fasta = open(filename)
    for line in fasta:
        if line.startswith('>') and name == None:
            name = line.rstrip('\n').replace('>','')
        else:
            if line.startswith('>') and name != None:
                yield [name, seq]
                name = line.rstrip('\n').replace('>','')
                seq = ''
            else:
                seq = seq + line.rstrip('\n')
    yield [name, seq]
    fasta.close()

# list fasta files
list_filenames = []
for file in os.listdir(args.input): 
    if file.endswith(".fasta"):
        list_filenames.append(file)

# sort filenames
list_filenames = sorted(list_filenames)

# create empty list to populate with pandas dataframes
list_pandas_df = []

# iterate through fasta files
for filename in list_filenames:
    # get name of gene from the file name
    gene = filename.replace('.fasta','')
    # read fasta as list of sequences
    fasta = read_fasta(f'{args.input}/{filename}')
    # convert list to pandas dataframe
    fasta_df = pd.DataFrame(data=fasta, columns = ['sample', 'sequence'])
    # replace gene names in the sample column
    # fasta_df = fasta_df.replace(f';{gene}', '' , regex = True)
    # set sample names as index
    fasta_df = fasta_df.set_index('sample')
    # rename sequence column as the gene name
    fasta_df = fasta_df.rename(columns = {'sequence':gene})
    # append to list of pandas dataframes
    list_pandas_df.append(fasta_df)

# merge dataframe by indexes
merge_df = reduce(lambda  left, right : pd.merge(left, right, how='outer', left_index=True, right_index=True), list_pandas_df)

# perform alignment checks and replace NaNs with gaps

# iterate through columns
for i in merge_df.columns:
    # get sequences as Series
    seq = merge_df[i]
    # sequence lengths excluding NaNs
    seq_len = seq.dropna().str.len()
    # max length
    max_len = max(seq_len)
    # check all sequences are the same length (aligned)
    if not all(seq_len == max_len):
        sys.exit(f'Error: Sequence lengths for {i} are variable')
    # replace NaN with gaps '-' the same length as the sequence
    seq[seq.isnull()] = '-'*max_len

# get gene lengths as series with name 'gene_lengths'
gene_lengths = merge_df.iloc[1].str.len().rename('gene_lengths')

# create dataframe of partitions
pt = pd.DataFrame(gene_lengths, index = merge_df.columns)
pt['end'] = pt['gene_lengths'].cumsum()
pt['start'] = (pt['end'] - pt['gene_lengths']) + 1
pt = pt[['start','end']]


# create output dir if not already presnet
if not os.path.exists(args.output):
    os.mkdir(args.output)
else:
    if not args.overwrite:
        sys.exit("Error: Output directory already exists. Remove or use --overwrite")
    else:
        shutil.rmtree(args.output)
        os.mkdir(args.output)

# write concatenated fasta
concat_fasta = open(f'{args.output}/{args.prefix}.fasta', 'w')
for index, row in merge_df.iterrows():
    concat_fasta.write(f'>{index}\n')
    concat_fasta.write(f'{row.str.cat(sep="")}\n')
concat_fasta.close()

# write partition file
partition = open(f'{args.output}/{args.prefix}.txt', 'w')
for index, row in pt.iterrows():
    partition.write(f'{args.type}, {index} = {row["start"]}-{row["end"]}\n')
partition.close()


