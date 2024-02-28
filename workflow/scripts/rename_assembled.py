import os
import argparse

# simple python script to rename assembled sequence
# a two column table is written with new and old names

# argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input",     help = "Input fasta",      required=True)
parser.add_argument("--sample",    help = "Sample name",      required=True)
parser.add_argument("--output",    help = "Output directory", required=True)
args = parser.parse_args()

# function to parse a fasta file
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

# function to write a remaned fasta and summary of name changes 
def rename_fasta(fasta, sample_name, output_prefix):
    # open output files
    new_fasta = open(f"{output_prefix}.fasta", "w")
    new_names = open(f"{output_prefix}.txt",   "w")
    # set contig number to iterate over  
    contig_number = 0
    # iterate through fasta
    for i in fasta:
        name, sequence = i[0], i[1]
        if "circular" in name:
            # define short name if sequence circular
            short_name = f"{sample_name}_circular"
        else: 
            # define short name if contig
            short_name = f"{sample_name}_contig{contig_number}"
        # write to output files
        new_fasta.write(f">{short_name}\n{sequence}\n")
        new_names.write(f"{short_name}\t{name}\n")
        # increase contig iterator
        contig_number += 1
    # close output files
    new_fasta.close()
    new_names.close()

# make output dir if not already present
if not os.path.exists(args.output):
    os.mkdir(args.output)

# read fasta
fas = read_fasta(args.input)

# rename and write fasta and summary file
rename_fasta(fas, args.sample, f"{args.output}/{args.sample}")

