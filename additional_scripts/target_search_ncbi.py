
import argparse
from Bio import Entrez, SeqIO
import os
import shutil
import sys

# set email
Entrez.email = "o.william.white@gmail.com"


# argparse

argparse_description = '''
Python script to count or download organelle assemblies or ribrosomal sequences from NCBI using using Biopython.
'''

argparse_usage = '''

# print counts of chloroplast assemblies for NCBI taxonomic id 3702
python target_search_ncbi.py --taxonomy 3702 --target chloroplast --email o.william.white@gmail.com --action count

# print counts of chloroplast assemblies for scientific name id "Arabidopsis thaliana"
python target_search_ncbi.py --taxonomy "Arabidopsis thaliana" --target chloroplast --email o.william.white@gmail.com --action count

# download chloroplast assemblies for "Camelineae" 
python target_search_ncbi.py --taxonomy "Camelineae" --target chloroplast --email o.william.white@gmail.com --action download --output test

# download the first 5 chloroplast assemblies for "Camelineae" and overwrite the original output 
python target_search_ncbi.py --taxonomy "Camelineae" --target chloroplast --email o.william.white@gmail.com --action download --max_n 5 --output test --overwrite
   
'''
parser = argparse.ArgumentParser(description=argparse_description, usage=argparse_usage)
parser.add_argument("--taxonomy",  help="Taxonomic name or NCBI ID.", required=True)
parser.add_argument('--target',    help="Target sequence type.", choices=["chloroplast", "mitochondrion", "ribosomal"], required=True)
parser.add_argument('--email',     help="Email for Entrez.", required=True)
parser.add_argument('--action',    help="Count or download sequences.", choices=["count", "download"], required=True)
parser.add_argument('--output',    help="Output directory. Only required for '--action download'.")
parser.add_argument('--max_n',     help="Maximum number of sequences to download. Only required for '--action download'.", default=None, type=int)
parser.add_argument('--overwrite', help="Overwrite output directory. Only required for '--action download'.", action="store_true")
parser.add_argument('--verbose',   help="Print additional information for error checking", action="store_true")
args = parser.parse_args()

# argparse checks
if args.action == "download" and args.output == None:
    sys.exit("Must specify output directory '--output' with '--action download'")

# functions

# create dir and overwrite if specified
def create_dir(path, overwrite, verbose=False):
    # path exists and overwrite is false
    if os.path.exists(path) and overwrite == False:
        sys.exit("Output directory already exists, remove directory or use --overwrite")
    # path exists and overwrite is true
    elif os.path.exists(path) and overwrite == True:
        if verbose:
            print("Overwriting output directory '" + path + "'\n")
        shutil.rmtree(path)
        os.mkdir(path)
        os.mkdir(path + "/fasta")
        if args.target != "ribosomal":
            os.mkdir(path + "/genbank")
    # path does not exist
    elif not os.path.exists(path):
        if verbose:
            print("Creating output directory '" + path + "'\n")
        os.mkdir(path)
        os.mkdir(path + "/fasta")
        if args.target != "ribosomal":
            os.mkdir(path + "/genbank")


# get taxonomic id from scientific name
def get_taxonomic_id(taxonomy):
    handle = Entrez.esearch(db="Taxonomy", term=taxonomy + "[Scientific Name]")
    record = Entrez.read(handle)
    return str(record["IdList"][0])


assert get_taxonomic_id("Arabidopsis thaliana") == "3702"


# get scientific name from taxonomic id
def get_scientific_name(taxonomy):
    handle = Entrez.efetch(db="Taxonomy", id=taxonomy)
    record = Entrez.read(handle)
    return record[0]["ScientificName"]


assert get_scientific_name("3702") == "Arabidopsis thaliana"


# get lineage from taxonomy id
def get_lineage(taxonomy_id, taxonomy_name):
    # efetch
    handle = Entrez.efetch(db="Taxonomy", id=taxonomy_id, retmode="xml")
    records = Entrez.read(handle)
    # get lineage
    lineage = records[0]["Lineage"].split("; ")[::-1]
    # append taxonomy_name to lineage
    lineage = [taxonomy_name] + lineage
    return lineage


# esearch for target sequences using taxonomy name
def target_search(rank, target):
    # define search term
    if target == "chloroplast" or target == "mitochondrion":
        search_term = rank + "[Organism] AND " + target + "[Title] AND complete genome[Title]"
    else:
        if target == "ribosomal":
            search_term = rank + "[Organism] AND 28S[Title]"
    if args.verbose:
        print("Search term = '" + search_term + "'", end = "\t")
    # esearch
    handle = Entrez.esearch(db="Nucleotide", term=search_term)
    record = Entrez.read(handle)
    # max retained records is 20 - repeat with retmax set to counts
    counts = record["Count"]
    handle = Entrez.esearch(db="Nucleotide", term=search_term, retmax=counts)
    record = Entrez.read(handle)
    return record


# write target assembly counts across lineage
def lineage_count(lineage, target):
    for rank in lineage:
        record = target_search(rank, target)
        counts = record["Count"]
        # write counts per rank
        print(rank + "\t" + str(counts))


# efetch genbank file using accession id
def efetch_seq(id, format, output_dir, verbose=False):
    # efetch
    handle = Entrez.efetch(db="Nucleotide", id=id, rettype=format, retmode="text")
    seq_record = SeqIO.read(handle, format)
    # define output path
    output_path = output_dir + "/" + seq_record.id + "." + format
    # write file
    if verbose:
        print("   " + seq_record.id + "." + format)
    SeqIO.write(seq_record, output_path, format)


# main code

# summarise input
if args.verbose:
    print("Searching NCBI for:")
    print("   Taxonomy: " + args.taxonomy)
    print("   Organelle: " + args.target)
    print("   Email: " + args.email + "\n")

# create output directory if not present
if args.action == "download":
    create_dir(args.output, args.overwrite, verbose=args.verbose)

# get ncbi taxonomy id and scientific name
try:
    int(args.taxonomy)
    if args.verbose:
        print("Taxonomy input is an NCBI ID, searching for for scientific name")
    taxonomy_id = args.taxonomy
    taxonomy_name = get_scientific_name(taxonomy_id)
except:
    if args.verbose:
        print("Taxonomy input is a name/string, searching for for NCBI ID")
    taxonomy_name = args.taxonomy
    taxonomy_id = get_taxonomic_id(taxonomy_name)

if args.verbose:
    print("   NCBI Id: " + taxonomy_id)
    print("   Scientific name: " + taxonomy_name + "\n")


# get full lineage
lineage = get_lineage(taxonomy_id=taxonomy_id, taxonomy_name=taxonomy_name)
if args.verbose:
    print("Full lineage: " + ", ".join(lineage) + "\n")

# action = count
if args.action == "count":
    lineage_count(lineage=lineage, target=args.target)

# action = download
if args.action == "download":
    target_ids = []
    record = target_search(taxonomy_name, args.target)
    # iterate through record ids and append record ids to target_ids
    for record_id in record["IdList"]:
        target_ids.append(record_id)
    # check if max_n was set
    if args.max_n != None and len(target_ids) > args.max_n:
        if args.verbose:
            print("Downloading the first " + str(args.max_n) + " assemblies for " + taxonomy_name)
        target_ids = target_ids[0:args.max_n]
    else:
        if args.verbose:
            print("Downloading " + str(len(target_ids)) + " assemblies for " + taxonomy_name)
    for target_id in target_ids:
        efetch_seq(id=target_id, format="fasta", output_dir=args.output + "/fasta",  verbose=args.verbose)
        if args.target != "ribosomal":
            efetch_seq(id=target_id, format="gb", output_dir=args.output + "/genbank", verbose=args.verbose)

if args.verbose:
    print("Job complete!")
