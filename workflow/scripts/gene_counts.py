#!/bin/python 

import argparse
import os
import re

usage = """

python gene_counts.py --input example_data/genes/ --output example_output.txt 

"""

description = """

Python script to count genes across samples. The input should be a directory containing '.fasta' files. The output should be a path to the text file containing the gene counts across samples. The script assumes that the fasta files are named as '<gene>.fasta' and sequences are named as using format of skim2phylo i.e., 'Zet_ZKP_1315332_circular;5628-7163;-;cox1' or 'Zet_ZKP_1315332'. 

"""

# argparse
parser = argparse.ArgumentParser(usage=usage, description=description)
parser.add_argument("--input",     help = "Input directory (one or more) containing '.fasta' files.", nargs = "*", required=True)
parser.add_argument("--output",    help = "Output text file with counts.",           required=True)
args = parser.parse_args()

# functions
def read_fasta(filename):
    name, seq = None,''
    fasta = open(filename, 'r')
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

def format_name(name):
    name = name.split(';')[0]
    name = re.sub("_contig\d*$|","", name)
    name = re.sub("_circular|","", name)
    return name

# create empty dictionary to hold gene (keys) and samples (values)
dict_genes = {}

# list files in input dir
for path in args.input:
    for file in os.listdir(path):
        # files ending with ".fasta"
        if file.endswith("fasta"):
            # get gene from file name
            gene = file.replace(".fasta","")
            # read fasta
            fasta = read_fasta(f"{path}/{file}")
            #iterate through sequnces in fasta
            for name, seq in fasta:
                # if long assembly name found, reformat to sample only
                if re.search("contig", name) or re.search("circular", name):
                    name = format_name(name)
                # add gene and name to dictionary
                if dict_genes.get(gene) == None:
                    dict_genes[gene] = [name]
                else:
                    dict_genes[gene].append(name)

# get dict_genes sorted by key
dict_genes_keys = list(dict_genes.keys())
dict_genes_keys.sort()
dict_genes_sorted = {i: dict_genes[i] for i in dict_genes_keys}

# create empty list to append unique sample names
sample_names = []

# append unique sample names to list
for key, value in dict_genes_sorted.items():
    for name in value:
        if name not in sample_names:
            sample_names.append(name)

# sort sample list
sample_names.sort()

# open output file
out = open(args.output, "w")

# define column names
column_names = ["sample"]
column_names.extend(list(dict_genes_sorted.keys()))

# write column names
out.write("\t".join(column_names) + "\n")

# write gene counts
for sample in sample_names:
    line = [sample]
    for key, value in dict_genes_sorted.items():
        count = value.count(sample)
        line.append(str(count))    
    out.write("\t".join(line) + "\n")
            
# close output file
out.close()

