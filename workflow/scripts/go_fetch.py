#!/usr/bin/env python3

import os
import sys
import shutil
import argparse
from Bio import Entrez, SeqIO
import time
import urllib
import subprocess
from random import Random

# argparse
argparse_description = """
Simple python script to fetch organelle or ribosomal reference sequences from NCBI for a given taxonomy.
"""

argparse_usage = """
./go_fetch.py --taxonomy 3702--target chloroplast --db genbank --min 5 --max 10 --output arabidopsis_chloroplast --overwrite --getorganelle --email user_email@example.com
"""

# argparse
parser = argparse.ArgumentParser(prog = "go_batch.py", description=argparse_description, usage=argparse_usage)
parser.add_argument("--taxonomy",     help="Taxonomy of rank to search for e.g. \"Arabidopsis\"", type=str, required=True)
parser.add_argument("--target",       help="Target sequence type.", choices=["chloroplast", "mitochondrion", "ribosomal", "ribosomal_complete"], required=True)
parser.add_argument("--db",           help="Database to search. Either refseq (NCBI) or genbank (INSDC). Default=refseq.", choices=["refseq", "genbank"], required=False, default="refseq")
parser.add_argument("--min",          help="Minimum number of target sequences to download.", type=int, required=False)
parser.add_argument("--max",          help="Maximum number of target sequences to download. Must be larger than --min.", type=int, required=False)
parser.add_argument("--seed",         help="Seed used for subsampling.", type=int, required=False)
parser.add_argument("--output",       help="Output directory.", required=False)
parser.add_argument("--overwrite",    help="Overwrite output directory.", action="store_true", required=False)
parser.add_argument("--getorganelle", help="Format seed and gene database for get organelle.", action="store_true", required=False)
parser.add_argument("--email",        help="Email for Entrez.", required=True)
parser.add_argument("--api",          help="API for NCBI.", type=str, required=False)
parser.add_argument("--version",      action="version", version='1.0.0')
args = parser.parse_args()

### additional checks

# additional dependency checks

try:
    cmd = 'python -c "from Bio import Entrez, SeqIO"'
    subprocess.call(cmd.split(" "), stderr=subprocess.DEVNULL)
except FileNotFoundError:
    sys.exit("Error: Biopython not in path")

try:
    subprocess.call(["trf"], stderr=subprocess.DEVNULL)
except FileNotFoundError:
    sys.exit("Error: trf not in path")

try:
    subprocess.call(["get_annotated_regions_from_gb.py"], stdout=subprocess.DEVNULL)
except FileNotFoundError:
    sys.exit("Error: get_annotated_regions_from_gb.py not in path")

# additional parameter checks

if args.max <= args.min:
    sys.exit("Error: --max must be larger than --min")


### set email
print(f"Using email {args.email}")
Entrez.email = args.email


### set api if given
if args.api != None:
    print(f"Using API key: {args.api}") 
    Entrez.api_key = args.api


### increase time between and number of tries used by entrez

# increase sleep time between tries
Entrez.sleep_between_tries = 20
# max tries 
Entrez.max_tries = 20


### functions

# create dir and overwrite if specified
def create_dir(dirpath, overwrite):
    if os.path.exists(dirpath):
        if overwrite == True:
            shutil.rmtree(dirpath)
        else:
            sys.exit(f"{dirpath} already exists. Remove or use --overwrite")
    os.mkdir(dirpath)
    os.mkdir(f"{dirpath}/fasta")
    os.mkdir(f"{dirpath}/genbank")

# get taxonomic id from scientific name
def get_taxonomic_id(taxonomy):
    try:
        handle = Entrez.esearch(db="Taxonomy", term=f"{taxonomy}[Scientific Name]")
        record = Entrez.read(handle)
    except urllib.error.HTTPError as e:
        if e.code == 400:
            print("HTTP Error 400: get_taxonomic_id bad request. Retrying in 10 seconds...")
            time.sleep(10)  # Wait for 10 seconds
            handle = Entrez.esearch(db="Taxonomy", term=f"{taxonomy}[Scientific Name]")
            record = Entrez.read(handle)
        else:
            sys.exit(f"HTTP Error {e}: get_taxonomic_id bad request. Exiting.")
    return str(record["IdList"][0])
assert get_taxonomic_id("Arabidopsis thaliana") == "3702"

# check if taxonomy id exists
def taxid_exists(taxid):
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid)
        record = Entrez.read(handle)
    except urllib.error.HTTPError as e:
        if e.code == 400:
            print("HTTP Error 400: taxid_exists bad request. Retrying in 10 seconds...")
            time.sleep(10)  # Wait for 10 seconds
            handle = Entrez.efetch(db="taxonomy", id=taxid)
            record = Entrez.read(handle)
        else:
            sys.exit(f"HTTP Error {e}: taxid_exists bad request. Exiting.")
    if record:
        return True
    else:
        return False
assert taxid_exists(3701) == True

# get scientific name from taxonomic id
def get_scientific_name(taxonomy):
    try:
        handle = Entrez.efetch(db="Taxonomy", id=taxonomy)
        record = Entrez.read(handle)
    except urllib.error.HTTPError as e:
        if e.code == 400:
            print("HTTP Error 400: get_scientific_name bad request. Retrying in 10 seconds...")
            time.sleep(10)  # Wait for 10 seconds
            handle = Entrez.efetch(db="Taxonomy", id=taxonomy)
            record = Entrez.read(handle)
        else: 
            sys.exit(f"HTTP Error {e}: get_scientific_name. Exiting.")
    return record[0]["ScientificName"]
assert get_scientific_name("3702") == "Arabidopsis thaliana"

# check if taxonomy name exists
def scientific_name_exists(taxonomy):
    try: 
        handle = Entrez.esearch(db="Taxonomy", term=f"{taxonomy}[Scientific Name]")
        record = Entrez.read(handle)
    except urllib.error.HTTPError as e: 
        if e.code == 400:
            print("HTTP Error 400: scientific_name_exists bad request. Retrying in 10 seconds...")
            time.sleep(10)  # Wait for 10 seconds
            handle = Entrez.esearch(db="Taxonomy", term=f"{taxonomy}[Scientific Name]")
            record = Entrez.read(handle)
        else:
            sys.exit(f"HTTP Error {e}: scientific_name_exists bad request. Exiting.")
    if int(record["Count"]) >= 1:
        return True
    else:
        return False
assert scientific_name_exists("Arabidopsis") == True

# get rank from taxid
def get_rank(taxid):
    try:
        # efetch
        handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
        record = Entrez.read(handle)
    except urllib.error.HTTPError as e:
        if e.code == 400:
            print("HTTP Error 400: get_rank bad request. Retrying in 10 seconds...")
            time.sleep(10)  # Wait for 10 seconds
            handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
            record = Entrez.read(handle)
        else:
            sys.exit(f"HTTP Error {e}: get_rank bad request. Exiting.")
    rank = record[0]["Rank"]
    return rank

# get lineage from taxid
def get_lineage(taxid):
    try:
        # efetch
        handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
        record = Entrez.read(handle)
    except urllib.error.HTTPError as e:
        if e.code == 400:
            print("HTTP Error 400: get_lineage bad request. Retrying in 10 seconds...")
            time.sleep(10)  # Wait for 10 seconds
            handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
            record = Entrez.read(handle)
        else:
            sys.exit(f"HTTP Error {e}: get_lineage bad request. Exiting.")
    # get lineage
    lineage = record[0]["Lineage"].split("; ")[::-1]
    # return lineage
    return lineage

#assert get_lineage(3701) == ["Camelineae", "Brassicaceae", "Brassicales", "malvids", "rosids", "Pentapetalae", "Gunneridae", "eudicotyledons", "Mesangiospermae", "Magnoliopsida", "Spermatophyta", "Euphyllophyta", "Tracheophyta", "Embryophyta", "Streptophytina", "Streptophyta", "Viridiplantae", "Eukaryota", "cellular organisms"]

def get_children(taxonomy):
    handle = Entrez.esearch(db="taxonomy", term=f"{taxonomy}[next level]", retmax=9999)
    record = Entrez.read(handle)
    children = []
    print(f"Found {len(record['IdList'])} children")
    if len(record['IdList']) > 10: 
        print("This might take some time :)")
    for i in record["IdList"]:
        children.append(get_scientific_name(i))
    return(children)
#assert get_children("Arabidopsis lyrata") == ["Arabidopsis petraea subsp. umbrosa", "Arabidopsis petraea subsp. septentrionalis", "Arabidopsis lyrata subsp. lyrata", "Arabidopsis lyrata subsp. petraea"]

def print_phylogeny(main_lineage, optional_children = []):
    spacer = ""
    for l in main_lineage:
        print(f"{spacer}{l}")
        if spacer == "":
            spacer = "|_" + spacer
        else:
            spacer = "  " + spacer
    if optional_children is not []:
        for c in optional_children:
            print(f"{spacer}{c}")

# generate search term
# target = "chloroplast", "mitochondrion", "ribosomal", "ribosomal_complete"
def search_term(taxid, target, db):
    # add taxid to term
    term = f"{taxid}[Organism]"
    if target == "chloroplast" or target == "mitochondrion":
        term += f" AND {target}[Title] AND complete genome[Title]"
    if target == "ribosomal":
        term = f"({taxid}[Organism] AND (28S[Title] OR 25S[Title])) OR ({taxid}[Organism] AND 18S[Title]) OR ({taxid}[Organism] AND 5.8S[Title])"
    if target == "ribosomal_complete":
        term += f" AND (28S[Title] OR 25S[Title]) AND 18S[Title] AND 5.8S[Title]"
    # refseq are derived from genbank but not part of
    if db == "refseq":
        term += f" AND refseq[filter]"
    # genbank is part of the International Nucleotide Sequence Database Collaboration (INSDC) along with the European Nucleotide Archive and the DNA Data Bank of Japan (DDBJ)
    if db == "genbank":
        term += f" AND ddbj_embl_genbank[filter]"
    return term 

# count the number of sequences on ncbi using the term generated
def entrez_esearch(input_term):
    # esearch
    handle = Entrez.esearch(db="Nucleotide", term=input_term, retmax=999)
    record = Entrez.read(handle)
    return record["IdList"]

# efetch 
def entrez_efetch(id, format, output_directory):
    handle = Entrez.efetch(db="Nucleotide", id=id, rettype=format, retmode="text")
    seq_record = SeqIO.read(handle, format)
    output_path = f"{output_directory}/{seq_record.id}.{format}"
    print(f"   Downloading {seq_record.id}.{format}")
    SeqIO.write(seq_record, output_path, format)

# subsample a dictionary of taxa (keys), and lists of accessions (values)
def subsample(input_dictionary, sample_limit, seed):
    # create copy of input_dictionary
    input_dictionary_copy = {key: value.copy() for key, value in input_dictionary.items()}
    # create ouput to populate with counts
    output_dictionary = {k: [] for k in input_dictionary_copy.keys()}
    # use a single random instance with proper seeding
    random_instance = Random(seed)
    while sample_limit > 0:
        # shuffle keys to ensure randomness
        keys = list(input_dictionary_copy.keys())
        random_instance.shuffle(keys)
        # iterate through keys
        for key in keys:
            if input_dictionary_copy[key]:
                value = input_dictionary_copy[key].pop(0)
                output_dictionary[key].append(value)
                sample_limit -= 1
                if sample_limit == 0:
                    break
    return output_dictionary

# main recursive search function
def recursive_search(taxonomy, lineage, target, db, min_th, max_th, idlist):

    # define search term
    term = search_term(taxonomy, target, db)
    print(f"Using search term = {term}")

    # count the number of sequences in the input idlist
    count_idlist_input = len(idlist)

    # esearch using term and return idlist of matching accessions
    idlist_esearch = entrez_esearch(term)

    # count the number of sequences in the idlist
    count_idlist_esearch = len(idlist_esearch)

    # get idlist for new sequences
    idlist_new = []
    for i in idlist_esearch:
        if i not in idlist:
            idlist_new.append(i)

    # count new sequences
    count_idlist_new = len(idlist_new)

    # combine input and new idlist
    idlist_combined = idlist.copy()
    idlist_combined.extend(idlist_new)

    # count running total
    count_idlist_total = len(idlist_combined)

    # print counts
    print(f"Sequences input: {count_idlist_input}")
    print(f"Sequences found: {count_idlist_esearch}")
    print(f"Sequences new:   {count_idlist_new}")
    print(f"Sequences total: {count_idlist_total}")

    # does the number of sequences meet the minimum threshold
    if count_idlist_total < min_th:
        print("Minimum threhold not reached, moving up within lineage")
        try: 
            taxonomy = lineage[lineage.index(taxonomy)+1]
            print(f"Next level is {taxonomy}\n")
            return recursive_search(taxonomy, lineage, target, db, min_th, max_th, idlist_combined)
        except IndexError:
            sys.exit("Error: cannot search any higher within lineage")
    else:

        print("Minimum threhold reached")

        # if maximum not exceeded, download all
        if count_idlist_total <= max_th:
            print(f"Maximum threshold not exceeded. Downloading {count_idlist_total} sequences\n")
            
            print("\nCreating output directory")
            create_dir(args.output, args.overwrite)

            # efetch
            for i in idlist_combined:
                entrez_efetch(i, "fasta", f"{args.output}/fasta")
                entrez_efetch(i, "gb",    f"{args.output}/genbank")

        # if maximum exceeded, download subsample
        else:

            # get subsample number required
            count_idlist_subsample = max_th - count_idlist_input

            # get taxonomic id
            taxid = get_taxonomic_id(taxonomy)

            # get taxonomic rank
            rank = get_rank(taxid)

            print(f"Maximum threshold exceeded. Subsampling {count_idlist_subsample} sequences from children\n")

            # get children
            children = get_children(taxonomy)

            if len(children) == 0 or rank == "species": 
                
                if len(children) == 0:

                    print("No children lineages. Must be a terminal rank, i.e species\n")

                if rank == "species":
                    # note that a entrez search term with a subspecies taxonomy e.g. "Lutra lutra chinensis" does not return any results
                    # avoid searching below species rank
                    print("No children linages at species rank or above.\n")

                    if len(children) >= 1: 
                        print_phylogeny([taxonomy], children)

                print(f"\nDownloading the first {count_idlist_subsample} sequences\n")

                print("Creating output directory")
                create_dir(args.output, args.overwrite)

                # efetch
                for i in idlist_combined[:count_idlist_subsample]:
                    entrez_efetch(i, "fasta", f"{args.output}/fasta")
                    entrez_efetch(i, "gb",    f"{args.output}/genbank")
                
            else:

                print("All children identified:\n")
 
                # print phylogeny
                print_phylogeny([taxonomy], children)

                # create dictionary of sequence ids from children
                dictionary_children = {}

                # iterate through children
                for c in children:

                    if c != "environmental samples":

                        term = search_term(c, target, db)

                        # esearch using term and return idlist of matching accessions
                        esearch_idlist = entrez_esearch(term)

                        if len(esearch_idlist) != 0:
                            for i in esearch_idlist:
                                if i not in idlist:
                                    # print(i)
                                    if dictionary_children.get(c) is None:
                                        dictionary_children[c] = [i]
                                    else:
                                        dictionary_children[c].append(i)

                # subsample dictionary
                dictionary_children_subsample = subsample(dictionary_children, count_idlist_subsample, args.seed)

                # print how many sequences found across children
                print("\nIdentified sequences:")
                for key, value, in dictionary_children.items():
                    print(f"   {len(value)}  {key}")

                # print how many sequences subsampled across children
                print("\nSampled sequences:")
                list_subsample = []
                for key, value, in dictionary_children_subsample.items():
                    print(f"   {len(value)}  {key}")
                    for v in value:
                        list_subsample.append(v)

                list_subsample

                print("\nCreating output directory")
                create_dir(args.output, args.overwrite)
            
                list_download = idlist.copy()
                list_download.extend(list_subsample)

                print(f"\nDownloading {len(list_download)} sequences")
                for i in list_download:
                    entrez_efetch(i, "fasta", f"{args.output}/fasta")
                    entrez_efetch(i, "gb",    f"{args.output}/genbank")


# cat files
def cat_files(input_dir, file_ending, output_file):
    # cmd to concatenate files
    cmd_cat = ["cat"]
    for f in os.listdir(input_dir):
        if f.endswith(file_ending):
             cmd_cat.append(f"{input_dir}/{f}")
    # subprocess run
    result_cat = subprocess.run(cmd_cat, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # capture stdout to text file
    with open(output_file, "w") as log:
        log.write(f"{result_cat.stdout}\n")

# format seed
def format_seed(path):
    print("\nFormating seed database")
    print("   Running trf")
    # get pwd
    pwd = os.getcwd()
    # cd to dir containing fasta files
    os.chdir(f"{path}/fasta")
    # iterate through fasta files
    for fasta in os.listdir():
        if fasta.endswith(".fasta"):
            # define cmd
            cmd_trf = f"trf {fasta} 2 7 7 80 10 50 500 -f -d -m -h"
            print(f"   {cmd_trf}")
            # subprocess run
            result_trf = subprocess.run(cmd_trf.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            # capture stdout and stderr to text file
            with open(f"{fasta}.log", "w") as log:
                log.write(f"stdout:\n{result_trf.stdout}\n")
                log.write(f"stderr:\n{result_trf.stderr}\n")
    # change back to orginal dir
    os.chdir(pwd)
    # cat masked fastas
    cat_files(f"{path}/fasta", ".mask",  f"{path}/seed.fasta")
   
# format gene
def format_gene(path, target):
    print("\nFormating gene database")
    print("    Running get_annotated_regions_from_gb.py")
    cmd_gar = ["get_annotated_regions_from_gb.py"]
    for gb in os.listdir(f"{path}/genbank/"):
        if gb.endswith(".gb"):
            cmd_gar.append(f"{path}/genbank/{gb}")
    if target == "mitochondrion" or target == "chloroplast":
        cmd_gar.extend(["-o", f"{path}/annotated_regions", "-t", "CDS", "--mix"])
    else:
        if target == "ribosomal" or target == "ribosomal_complete":
            cmd_gar.extend(["-o", f"{path}/annotated_regions", "-t", "rRNA", "--mix"])
    # subprocess run
    result_gar = subprocess.run(cmd_gar, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # capture stdout and stderr to text file
    with open(f"{path}/annotated_regions/get_annotated_regions_from_gb.log", "w") as log:
        log.write(f"stdout:\n{result_gar.stdout}\n")
        log.write(f"stderr:\n{result_gar.stderr}\n")
    # cp gene file to sample dir
    subprocess.run(["cp", f"{path}/annotated_regions/gene/gene.fasta", f"{path}/"])


### main

# get and check ncbi taxonomy_id and scientific_name
print("Running go_fetch")
print("Checking user input\n")
try:
    int(args.taxonomy)
    print("Taxonomy input is a numeric NCBI ID.")
    # check if taxonomy id exists
    if not taxid_exists(args.taxonomy):
        sys.exit(f"Taxonomy ID {args.taxonomy} does not exist. Please check for the correct taxonomy id on https://www.ncbi.nlm.nih.gov/taxonomy")
    else:
        print("Taxonomy ID exists on NCBI")
    taxonomy_id = args.taxonomy
    taxonomy_name = get_scientific_name(taxonomy_id)
except:
    print("Taxonomy input is a scientific name.")
    if not scientific_name_exists(args.taxonomy):
        sys.exit(f"Scientific name {args.taxonomy} does not exist. Please check for the correct taxonomy id on https://www.ncbi.nlm.nih.gov/taxonomy")
    taxonomy_name = args.taxonomy
    taxonomy_id = get_taxonomic_id(taxonomy_name)

print(f"NCBI Id: {taxonomy_id}")
print(f"Scientific name: {taxonomy_name}")

# get lineage
print("\nChecking lineage\n")

lineage = get_lineage(taxonomy_id)
lineage.insert(0, taxonomy_name)

# get children
# children = get_children(taxonomy_name)
# only request child lineages if needed to reduce entrez searches

# print phylogeny
# print_phylogeny(lineage[::-1], optional_children = children)
print_phylogeny(lineage[::-1])

print("\nStarting search\n")

# start recursive search function                
recursive_search(taxonomy_name, lineage, args.target, args.db, args.min, args.max, [])

# format gene and seed database if requested
if args.getorganelle:
    format_seed(args.output)
    format_gene(args.output, args.target)

print("\ngo_fetch complete!")


