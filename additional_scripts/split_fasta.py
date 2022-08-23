import argparse
import sys
import gzip
import os

# argparse

argparse_description = '''
Python script to split a fasta file into indivdual_sequences
'''

argparse_usage = '''
python split_fasta.py --inpuut <input_fasta> --output <output_directory>
'''

parser = argparse.ArgumentParser(description=argparse_description, usage=argparse_usage)
parser.add_argument("--input",  help="Input fasta to split",            required=True)
parser.add_argument('--output', help="Output directory to write files to", required=True)
args = parser.parse_args()

# function to edit sequence id
def edit_name(name):
    name = name.replace('|', '_').replace(':', '').replace(',','').replace(' ', '_').replace('(', '').replace(')', '').replace('+', '_').replace('-', '_')
    if name[-1] == '.':
        name = name[:-1]
    return(name)

# function to split fasta and write to files to output folder
def split_fasta(fasta, output_dir):
    fasta = open(fasta, 'r')
    out_fasta = ''
    for line in fasta:
        if line.startswith('>'):
            if out_fasta != '':
                out_fasta.close()
            seq_name = line[1:].rstrip('\n')
            file_name = edit_name(seq_name)
            out_fasta = open(output_dir + '/' + file_name + '.fasta', 'w')
            out_fasta.write('>'+ seq_name +'\n')
        else:
            out_fasta.write(line)
    out_fasta.close()

# split fasta
split_fasta(args.input, args.output)

