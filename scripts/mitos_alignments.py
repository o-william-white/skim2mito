import os
import sys
import subprocess
import re

mitos_dir = sys.argv[1]
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

# read multi-line fasta
def read_fasta(input_path):
    name, seq = None,""
    with open(input_path, "r") as fasta:
        for line in fasta:
            if line.startswith('>') and name == None:
                name = line.rstrip('\n').replace('>','').replace(' ', '')
            else:
                if line.startswith('>') and name != None:
                    yield (name, seq)
                    name = line.rstrip('\n').replace('>','').replace(' ', '')
                    seq = ''
                else:
                    seq = seq + line.rstrip('\n')
        yield (name, seq)

# read all sequence names and sequences across fastas into a list
seq_list = []
for fas in file_paths(mitos_dir, ".fas"):
    for i in read_fasta(fas):
        seq_list.append(i)

# get a list of all possible genes in annotation 
# note I ignored _1 _2 and -a -b to denote multiple copies or parts
# I also removed annotations for OH, OL and trn
genes = []
for seq in seq_list:
    gene = seq[0].split(";")[3]
    gene = gene.split("_")[0].split("-")[0]
    if gene not in genes and not re.search("trn|OH|OL", gene):
        genes.append(gene)
genes = sorted(genes)

# create files with sequences for each gene
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
for gene in genes:
    with open(output_dir + "/" + gene + ".fasta", "w") as output_fasta:
        for seq in seq_list:
            if re.search(gene+"$|"+gene+"_[0-9]$|"+gene+"-[a-z]$", seq[0]):
                output_fasta.write(">" + seq[0] + "\n" + seq[1] + "\n")

### create a table of gene presence/absence for each contig
table_list = []                        # list to populate
table_list.append(["sample"] + genes)  # column names
# loop through bed files and genes
for bed_path in file_paths(mitos_dir, ".bed"):
    sample = bed_path.replace(mitos_dir, "").replace("/result.bed", "").replace("\\", "-")
    sample_list = [sample] # list to populate / row
    for gene in genes: 
        with open(bed_path, "r") as bed:
            count = 0
            for line in bed:
                if gene in line.rstrip("\n").split("\t")[3]:
                    count += 1
            sample_list.append(count)
    table_list.append(sample_list)

# write to table
with open(f"{output_dir}/summary.txt", "w") as gene_summary:
    for i in table_list:
        gene_summary.write("\t".join(map(str, i)) + "\n")

