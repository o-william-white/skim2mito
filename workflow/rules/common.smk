import pandas as pd
import sys
from snakemake.utils import min_version


# define min version
min_version("8.4.12")


# set configfile
configfile: "config/config.yaml"


# configfile parameters
user_email = config["user_email"]
go_reference = config["go_reference"]
forward_adapter = config["forward_adapter"]
reverse_adapter = config["reverse_adapter"]
fastp_dedup = config["fastp_dedup"]
mitos_refseq = config["mitos_refseq"]
mitos_code = config["mitos_code"]
alignment_trim = config["alignment_trim"]
missing_threshold = config["missing_threshold"]
outgroup = config["outgroup"]
plot_height = config["plot_height"]
plot_width = config["plot_width"]

# read sample data
if os.path.exists(config["samples"]):
    sample_data = pd.read_csv(config["samples"]).set_index("ID", drop=False)
else:
    sys.exit(f"Error: samples.csv file '{config['samples']}' does not exist")


# functions to get metadata sample list
def get_forward(wildcards):
    return sample_data.loc[wildcards.sample, "forward"]


def get_reverse(wildcards):
    return sample_data.loc[wildcards.sample, "reverse"]


def get_fastq(wildcards):
    fwd = sample_data.loc[wildcards.sample, "forward"]
    rev = sample_data.loc[wildcards.sample, "reverse"]
    return [fwd, rev]


def get_taxid(wildcards):
    return sample_data.loc[wildcards.sample, "taxid"]


def get_seed(wildcards):
    return sample_data.loc[wildcards.sample, "seed"]


def get_gene(wildcards):
    return sample_data.loc[wildcards.sample, "gene"]


# functions for checkpoints
def get_seqkit_output(wildcards):
    ck_output = checkpoints.assembled_sequence.get(**wildcards).output[0]
    return expand(
        rules.seqkit.output,
        sample=glob_wildcards(os.path.join(ck_output, "{sample}.fasta")).sample,
    )

def get_blobtools_output(wildcards):
    ck_output = checkpoints.assembled_sequence.get(**wildcards).output[0]
    return expand(
        rules.blobtools.output,
        sample=glob_wildcards(os.path.join(ck_output, "{sample}.fasta")).sample,
    )

def get_assess_assembly_output(wildcards):
    ck_output = checkpoints.assembled_sequence.get(**wildcards).output[0]
    return expand(
        rules.assess_assembly.output,
        sample=glob_wildcards(os.path.join(ck_output, "{sample}.fasta")).sample,
    )

def get_annotated_samples(wildcards):
    ck_output = checkpoints.assembled_sequence.get(**wildcards).output[0]
    return expand(
        rules.annotations.output,
        sample=glob_wildcards(os.path.join(ck_output, "{sample}.fasta")).sample,
    )

def get_sucessfully_annotated_samples(wildcards):
    ck_output = checkpoints.extract_annotated_genes.get(**wildcards).output[1]
    # get the unique values (before the "/" character) in the first column, exclude the header
    annotated_samples = list(set(l.strip().split('\t')[0].split('/')[0]
                                 for i, l in enumerate(open(ck_output).readlines()) if i != 0))
    return expand(
        rules.annotations.output,
        sample=annotated_samples,
    )

def get_mafft_output(wildcards):
    checkpoint_output = checkpoints.extract_annotated_genes.get(**wildcards).output[0]
    return expand(
        "results/mafft/{i}.fasta",
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i,
    )

def get_plot_tree_output(wildcards):
    checkpoint_output = checkpoints.extract_annotated_genes.get(**wildcards).output[0]
    return expand(
        "results/plot_tree/{i}.png",
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i,
    )

def get_mafft_filtered_output(wildcards):
    checkpoint_output = checkpoints.extract_annotated_genes.get(**wildcards).output[0]
    return expand(
        "results/mafft_filtered/{i}.fasta",
        i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i,
    )


# config paramter checks
if go_reference == "go_fetch" and user_email == "user@example_email.com":
    sys.exit(
        f"Error: if using go_fetch to download references, please change the example email provided in the config file'"
    )
if go_reference != "go_fetch" and go_reference != "custom":
    sys.exit(f"Error: go_reference must be 'go_fetch' or 'custom'")
if not isinstance(fastp_dedup, bool):
    sys.exit(f"Error: fastp_dedup must be 'True' or 'False'")
if mitos_refseq not in [
    "refseq39",
    "refseq63f",
    "refseq63m",
    "refseq63o",
    "refseq89f",
    "refseq89m",
    "refseq89o",
]:
    sys.exit(
        "Error: mitos_refseq must be one of 'refseq39', 'refseq63f', 'refseq63m', 'refseq63o', 'refseq89f', 'refseq89m', 'refseq89o'"
    )
if mitos_code not in [2, 4, 5, 9, 13, 14]:
    sys.exit("Error: mitos_code must be one of 2, 4, 5, 9, 13, 14")
if (
    not isinstance(missing_threshold, float)
    or missing_threshold < 0.0
    or missing_threshold > 1.0
):
    sys.exit("Error: missing_threshold must be a float between 0.0 and 1.0")
if alignment_trim not in ["gblocks", "clipkit"]:
    sys.exit("Error: alignment_trim must be 'gblocks' or 'clipkit'")

# samples.csv check
if any(sample_data["ID"].duplicated()):
    sys.exit(
        f"Error: duplicated sample names present: {list(sample_data['ID'][sample_data['ID'].duplicated()])}"
    )
for i in sample_data["forward"]:
    if not os.path.exists(i):
        sys.exit(f"Error: forward reads path '{i}' does not exist")
for i in sample_data["reverse"]:
    if not os.path.exists(i):
        sys.exit(f"Error: reverse reads path '{i}' does not exist")
if go_reference == "custom":
    for i in sample_data["seed"].unique():
        if not os.path.exists(i):
            sys.exit(f"Error: seed database path '{i}' does not exist")
    for i in sample_data["gene"].unique():
        if not os.path.exists(i):
            sys.exit(f"Error: gene database path '{i}' does not exist")

wildcard_constraints:
    sample=r"[^*/~]+",