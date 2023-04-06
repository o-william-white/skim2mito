
# Genome skimming pipeline

This repo is a Snakemake pipeline to assemble organelle or ribsomal sequences from genome skimming data using GetOrganelle, check assembly quality, annotate genes and perfrom basic phylogenetic analyses. This is under active development so issues are likey. Currently, it can assemble and annotate mitochondrial and ribosomal sequences. Chloroplast sequences can be assembled but annotation is not possible yet.  

## Setup 

The pipeline uses conda environments to install the necessary tools. The pipeline also requires reference data for GetOrganelle, a blastn database, an NCBI taxdump and reference data for MITOS2 gene annotation. Scripts are provided to install the necessary reference databases. 

```
# get github repo
git clone https://github.com/o-william-white/genome_skimming_pipeline

# change dir
cd genome_skimming_pipeline

# setup conda env
conda env create -n genome_skimming_pipeline -f genome_skimming_pipeline/envs/conda_env.yaml

# get blast database
# TBC - already available on NHM HPC :)

# get new_taxdump
bash additional_scripts/wget_new_taxdump.sh

# get mitos2 annotation
bash additional_scripts/fetch_mitos2_reference_data.sh
```

## Download example data

We will use a subset of eight Apis mellifera samples from Parejo et al 2020 https://academic.oup.com/gbe/article/12/12/2535/5900668

```
# set up sra-tools environment
conda env create -n sra-tools -c bioconda

# download reads
bash additional_scripts/fastq_dump_example_data.sh
```

## Download GetOrganelle reference data

GetOrganelle requires reference data in the format of seed and gene reference fasta files. You can use the default reference data for GetOrganelle, but I would recomend using custom reference databases where possible. See here for details of how to set up your own databases https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#how-to-assemble-a-target-organelle-genome-using-my-own-reference 

I have shared a basic python script called go_fetch in another repo https://github.com/o-william-white/go_fetch to download and format reference data formated for GetOrganelle. Go fetch downloads the reference from NCBI using biopython, removes repetive sequences using trf, and formats the data for GetOrganelle.

```
# set up go_fetch environment
conda create -n go_fetch -c bioconda getorganelle biopython trf

# activate environment
conda activate go_fetch

# download refseq mitochondrion sequences for Apis mellifera
# NOTE: You will need to change the path go_fetch.py and the email address used by Biopython. 
python /PATH/TO/go_fetch/go_fetch.py \
   --taxonomy "Apis mellifera" \
   --target mitochondrion \
   --download \
   --min 2 \
   --max 10 \
   --name mitochondrion \
   --output apis_mellifera_db \
   --email o.william.white@gmail.com

# deactivate environment
conda deactivate
```

You should now have the following refernece data: `apis_mellifera_db/mitochondrion/seed.fasta` and `apis_mellifera_db/mitochondrion/gene.fasta`

## Setup Snakemake config files

Sankemake requires a config.yaml and samples.csv to define input paramters and sequence data for each sample.

```
# setup example config.yaml and samples.csv
bash additional_scripts/setup_example_config.sh
```

This will create a file called config/config_example.yaml that should look like this: 
```
# path to sample sheet csv with columns for ID,forward,reverse
samples: config/samples_example.csv

# name of output directory
output_dir: results_example

# GetOrganelle target sequence type (animal_mt,embplant_cp,anonym)
target_type: animal_mt

# path to blast db
blast_db: /workspaces/groups/database/nt-2022-06-22_blastdb/nt

# path to new_taxdump
taxdump: taxdump

# path to mitos refseq
mitos_refseq: mitos2_reference_data

# mito code
mitos_code: 5

# barrnap kindgom (Bacteria:bac, Archaea:arc, Eukaryota:euk, Metazoan Mitochondria:mito)
barrnap_kingdom: NA

# number of threads to use
threads: 6
```

It will also create csv table of samples that looks like this: 

|ID|forward|reverse|seed|gene|
|--|------ | ----- | -- | -- |
|A.mellifera_ERR5665185|example_data/ERR5665185_1.fastq.gz|example_data/ERR5665185_2.fastq.gz|apis_mellifera_db/mitochondrion/seed.fasta|apis_mellifera_db/mitochondrion/gene.fasta|
|A.mellifera_ERR5665186|example_data/ERR5665186_1.fastq.gz|example_data/ERR5665186_2.fastq.gz|apis_mellifera_db/mitochondrion/seed.fasta|apis_mellifera_db/mitochondrion/gene.fasta|
|A.mellifera_ERR5665187|example_data/ERR5665187_1.fastq.gz|example_data/ERR5665187_2.fastq.gz|apis_mellifera_db/mitochondrion/seed.fasta|apis_mellifera_db/mitochondrion/gene.fasta|
|A.mellifera_ERR5665184|example_data/ERR5665184_1.fastq.gz|example_data/ERR5665184_2.fastq.gz|apis_mellifera_db/mitochondrion/seed.fasta|apis_mellifera_db/mitochondrion/gene.fasta|
|A.mellifera_ERR5665183|example_data/ERR5665183_1.fastq.gz|example_data/ERR5665183_2.fastq.gz|apis_mellifera_db/mitochondrion/seed.fasta|apis_mellifera_db/mitochondrion/gene.fasta|
|A.mellifera_ERR5665182|example_data/ERR5665182_1.fastq.gz|example_data/ERR5665182_2.fastq.gz|apis_mellifera_db/mitochondrion/seed.fasta|apis_mellifera_db/mitochondrion/gene.fasta|
|A.mellifera_ERR5665181|example_data/ERR5665181_1.fastq.gz|example_data/ERR5665181_2.fastq.gz|apis_mellifera_db/mitochondrion/seed.fasta|apis_mellifera_db/mitochondrion/gene.fasta|
|A.mellifera_ERR5665180|example_data/ERR5665180_1.fastq.gz|example_data/ERR5665180_2.fastq.gz|apis_mellifera_db/mitochondrion/seed.fasta|apis_mellifera_db/mitochondrion/gene.fasta|


## Run Snakemake

Now we have set up our conda environemt, downloaded the reference data, downloaded the input data and formated the snakemake config, we are ready to run the pipeline. 

The first time you run the pipeline, it will start by setting up the singularity and conda environemts. If you run this Snakemake file again, it will not need to setup these environments again. 

You can run the pipeline by entering:
```
source activate genome_skimming_pipeline

snakemake \
   --cores 6 \
   --use-conda \
   --use-singularity \
   --configfile config/config_example.yaml
```

If you have access to an HPC you can submit this as a job. For example, see the slurm script below. 
```
#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_snakemake_example_%j.out
#SBATCH --error=job_snakemake_example_%j.err
#SBATCH --mem=18G
#SBATCH --cpus-per-task=18

source activate genome_skimming_pipeline

snakemake \
   --cores 18 \
   --use-conda \
   --use-singularity \
   --configfile config/config_example.yaml

echo Complete!
```

