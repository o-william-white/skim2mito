import os
import sys
import shutil
import argparse
from Bio import Entrez, SeqIO

# argparse
argparse_description = '''
Python script to takes a tab seaparate table with sample names in the first column and lineage in remaining columns.
For each sample (row), the script will determine the rank at which to download target sequences based on the minimum threshold specified.
Sequences are downloaded from NCBI using using Biopython.
'''
argparse_usage = '''
TBC
'''
parser = argparse.ArgumentParser(description=argparse_description, usage=argparse_usage)
parser.add_argument("--sample_table", help="Tab separate file samples table", required=True)
parser.add_argument("--output_table", help="Path to output table", required=True)
parser.add_argument('--target',       help="Target sequence type.", choices=["chloroplast", "mitochondrion", "ribosomal"], required=True)
parser.add_argument('--email',        help="Email for Entrez.", required=True)
parser.add_argument('--output',       help="Output directory.", required=True)
parser.add_argument('--threshold_n',  help="Threshold for number of sequences to download.", type=int, required=True)
parser.add_argument('--max_n',        help="Maximum number of sequences to download.",       type=int, required=True)
args = parser.parse_args()

# set email
Entrez.email = args.email

# functions

# create dir and overwrite if specified
def create_dir(path, target, overwrite):
    if overwrite == True:
        shutil.rmtree(path)
    os.mkdir(path)
    os.mkdir(f"{path}/fasta")
    if target != "ribosomal":
        os.mkdir(f"{path}/genbank")

def taxonomy_exists(taxonomy):
    handle = Entrez.esearch(db="Taxonomy", term=taxonomy + "[Scientific Name]")
    record = Entrez.read(handle)
    count = int(record["Count"])
    if count != 0:
        return True
    else:
        return False

def taxonomy_search(rank, target, max_n):
    # define search term
    if target == "chloroplast" or target == "mitochondrion":
        search_term = rank + "[Organism] AND " + target + "[Title] AND complete genome[Title]"
    else:
        if target == "ribosomal":
            search_term = rank + "[Organism] AND 28S[Title]"
    # esearch
    handle = Entrez.esearch(db="Nucleotide", term=search_term, retmax=max_n)
    record = Entrez.read(handle)
    return record

def efetch_seq(id, format, output_dir):
    handle = Entrez.efetch(db="Nucleotide", id=id, rettype=format, retmode="text")
    seq_record = SeqIO.read(handle, format)
    output_path = f"{output_dir}/{seq_record.id}.{format}"
    print(f"    Downloading {seq_record.id}.{format}")
    SeqIO.write(seq_record, output_path, format)

# open input table
sample_table = open(args.sample_table, "r")

# make output dir if not already present
if not os.path.exists(args.output):
    os.mkdir(args.output)

# open output table
output_table = open(args.output_table, "w")

# iterate through sample list
for line in sample_table:
    line = line.rstrip("\n").split("\t")
    sample, lineage = line[0], line[1:]
    lineage.reverse()
    threshold_met = False

    print(f"\n\n***Starting new search***")
    print(f"\nCounting the number of {args.target} sequences for {sample} across the lineage:")
    print(f"   {';'.join(lineage)}\n")
    # set target ids list
    # setting target ids list here means that more closely related references that do not meet threshold number are retained
    target_ids = []
    # iterate through lineage
    for taxonomy in lineage:
        if threshold_met == False:
            # check if taxonomy exists in ncbi
            if not taxonomy_exists(taxonomy):
                print(taxonomy + "\tNA")
            else:
                # search for target sequence
                record = taxonomy_search(taxonomy, args.target, args.max_n)
                # get count of target sequence
                count = int(record["Count"])
                print(taxonomy + "\t" + str(count))
                # append target ids
                for record_id in record["IdList"]:
                    target_ids.append(record_id)
                # check if count meets threshold for number of sequences 
                if len(target_ids) >= args.threshold_n:
                    print(f"\nThreshold of {args.threshold_n} met for {taxonomy}")
                    output_table.write(f"{sample}\t{taxonomy}\n")
                    threshold_met = True
                    # print the number of sequences that will be downloaded
                    if len(target_ids) > args.max_n:
                        target_ids = target_ids[0:args.max_n]
                        print(f"\nDownloading the first {args.max_n} sequences")
                    else:
                        print(f"\nDownloading count sequences")
                    # check if directory has already been created and all files downloaded
                    if os.path.exists(f"{args.output}/{taxonomy}"):
                        files = [f for f in os.listdir(f"{args.output}/{taxonomy}/fasta/") if f.endswith('.fasta')] # only checks that all fasta files present
                        if len(files) == count or count > args.max_n and len(files) == args.max_n:
                            print(f"\nDirectory with downloaded files already present for {taxonomy}")
                            break
                        else:
                            print(f"\nOverwriting directory with incomplete number of downloaded files for {taxonomy}")
                            create_dir(path=f"{args.output}/{taxonomy}", target=args.target, overwrite=True)
                    else:
                        print(f"\nCreating new directory for {taxonomy}")
                        create_dir(path=f"{args.output}/{taxonomy}", target=args.target, overwrite=False)
                    # download sequences
                    for target_id in target_ids:
                        efetch_seq(id=target_id, format="fasta", output_dir=f"{args.output}/{taxonomy}/fasta")
                        if args.target != "ribosomal":
                            efetch_seq(id=target_id, format="gb", output_dir=f"{args.output}/{taxonomy}/genbank")      

# close sample list and output table  
sample_table.close()
output_table.close()


