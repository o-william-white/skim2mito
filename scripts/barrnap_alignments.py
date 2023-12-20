import os
import sys
import subprocess
import re

barrnap_dir = sys.argv[1]
output_dir = sys.argv[2]

# get paths for files in directory with specific endings
def file_paths(directory, ending):
    output_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(ending):
                output_list.append(root + "/" + file)
    return output_list

# e.g.
# bed_paths = file_paths("get_organelle_mt_mitos", ".bed")
# fas_paths = file_paths("get_organelle_mt_mitos", ".fas")

# function to fomat name to match mitos format
def format_name(name):
    gene, sample_name = name.split("::")
    sample_name = sample_name.replace(":", ";")
    return(sample_name + ";" + gene)

# read multi-line fasta
def read_fasta(input_path):
    name, seq = None,""
    with open(input_path, "r") as fasta:
        for line in fasta:
            if line.startswith('>') and name == None:
                name = line.rstrip('\n').replace('>','').replace(' ', '')
                name = format_name(name)
            else:
                if line.startswith('>') and name != None:
                    yield (name, seq)
                    name = line.rstrip('\n').replace('>','').replace(' ', '')
                    name = format_name(name)
                    seq = ''
                else:
                    seq = seq + line.rstrip('\n')
        yield (name, seq)

# read all sequence names and sequences across fastas into a list
seq_list = []
for fas in file_paths(barrnap_dir, ".fas"):
    for i in read_fasta(fas):
        if i[0] != None: # ignore empty fasta files
            seq_list.append(i)

# get a list of all possible genes in annotation 
genes = []
for seq in seq_list:
    gene = seq[0].split(";")[2]
    #gene = gene.split("_")[0].split("-")[0]
    if gene not in genes:
        genes.append(gene)
genes = sorted(genes)

# print(genes)

# create files with sequences for each gene
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
for gene in genes:
    with open(output_dir + "/" + gene + ".fasta", "w") as output_fasta:
        for seq in seq_list:
            if re.search(gene, seq[0]):
                output_fasta.write(">" + seq[0] + "\n" + seq[1] + "\n")

### create a table of gene presence/absence for each contig
table_list = []                        # list to populate
table_list.append(["sample"] + genes)  # column names
# loop through bed files and genes
for bed_path in file_paths(barrnap_dir, ".gff"):
    sample = bed_path.replace(barrnap_dir, "").replace("/result.gff", "").replace("\\", "-")
    sample_list = [sample] # list to populate / row
    for gene in genes: 
        with open(bed_path, "r") as bed:
            count = 0
            for line in bed:
                if not line.startswith("#"):
                    if gene in line.rstrip("\n").split("\t")[8]:
                        count += 1
            sample_list.append(count)
    table_list.append(sample_list)

# write to table
with open(f"{output_dir}/summary.txt", "w") as gene_summary:
    for i in table_list:
        gene_summary.write("\t".join(map(str, i)) + "\n")

