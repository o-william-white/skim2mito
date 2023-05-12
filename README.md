
# Genome skimming pipeline

This repo is a Snakemake pipeline to assemble organelle or ribsomal sequences from genome skimming data using GetOrganelle, check assembly quality, annotate genes and perform basic phylogenetic analyses. This is under active development so issues are likey. Currently, it can assemble and annotate mitochondrial and ribosomal sequences. Chloroplast sequences can be assembled but annotation is not possible yet.  

## Setup 

The pipeline uses conda environments to install the necessary tools. The pipeline also requires reference data for GetOrganelle, a blastn database, an NCBI taxdump and reference data for MITOS2 gene annotation. Scripts are provided to install the necessary reference databases. 

```
# get github repo
git clone https://github.com/o-william-white/genome_skimming_pipeline

# change dir
cd genome_skimming_pipeline

# setup conda env
conda env create -n genome_skimming_pipeline -f envs/conda_env.yaml

# get blast database
# TBC - already available on NHM HPC :)

# get new_taxdump
bash additional_scripts/fetch_new_taxdump.sh

# get mitos2 annotation
bash additional_scripts/fetch_mitos2_reference_data.sh
```

## Download example data

We will use a subset of nine bank vole (Myodes glareolus) samples from Baker et al (2017) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5680428/

```
# set up sra-tools environment
conda create -n sra-tools -c bioconda sra-tools

# download reads
bash additional_scripts/fastq_dump_example_data.sh
```

## Download GetOrganelle reference data

GetOrganelle requires reference data in the format of seed and gene reference fasta files. You can use the default reference data for GetOrganelle, but I would recomend using custom reference databases where possible. See here for details of how to set up your own databases https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#how-to-assemble-a-target-organelle-genome-using-my-own-reference 

I have shared a basic python script called go_fetch.py in another repo https://github.com/o-william-white/go_fetch to download and format reference data formatted for GetOrganelle. Go fetch downloads the reference from NCBI using biopython, removes repetitive sequences using trf, and formats the data for GetOrganelle.

```
# set up go_fetch environment
conda create -n go_fetch -c bioconda getorganelle biopython trf

# activate environment
conda activate go_fetch

# clone go_fetch.py script, NOTE: in this case in $HOME
cd ~ && git clone https://github.com/o-william-white/go_fetch.git && cd -

# download refseq mitochondrion sequences for Myodes
# NOTE: You will need to change the email address used by Biopython. 
python $HOME/go_fetch/go_fetch.py \
   --taxonomy "Myodes" \
   --target mitochondrion \
   --download \
   --min 2 \
   --max 10 \
   --name mitochondrion \
   --output myodes_db \
   --email o.william.white@gmail.com

# deactivate environment
conda deactivate
```

You should now have the following reference data: `myodes_db/mitochondrion/seed.fasta` and `myodes_db/mitochondrion/gene.fasta`

## Setup Snakemake config files

Snakemake requires a config.yaml and samples.csv to define input paramters and sequence data for each sample.

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
|--|-------|-------|----|----|
|M.glareolus_SRR5201684|example_data/SRR5201684_1.fastq.gz|example_data/SRR5201684_2.fastq.gz|myodes_db/mitochondrion/seed.fasta|myodes_db/mitochondrion/gene.fasta|
|M.glareolus_SRR5201683|example_data/SRR5201683_1.fastq.gz|example_data/SRR5201683_2.fastq.gz|myodes_db/mitochondrion/seed.fasta|myodes_db/mitochondrion/gene.fasta|
|M.glareolus_SRR5201682|example_data/SRR5201682_1.fastq.gz|example_data/SRR5201682_2.fastq.gz|myodes_db/mitochondrion/seed.fasta|myodes_db/mitochondrion/gene.fasta|
|M.glareolus_SRR5201681|example_data/SRR5201681_1.fastq.gz|example_data/SRR5201681_2.fastq.gz|myodes_db/mitochondrion/seed.fasta|myodes_db/mitochondrion/gene.fasta|
|M.glareolus_SRR5201680|example_data/SRR5201680_1.fastq.gz|example_data/SRR5201680_2.fastq.gz|myodes_db/mitochondrion/seed.fasta|myodes_db/mitochondrion/gene.fasta|
|M.glareolus_SRR5201679|example_data/SRR5201679_1.fastq.gz|example_data/SRR5201679_2.fastq.gz|myodes_db/mitochondrion/seed.fasta|myodes_db/mitochondrion/gene.fasta|
|M.glareolus_SRR5201678|example_data/SRR5201678_1.fastq.gz|example_data/SRR5201678_2.fastq.gz|myodes_db/mitochondrion/seed.fasta|myodes_db/mitochondrion/gene.fasta|
|M.glareolus_SRR5201677|example_data/SRR5201677_1.fastq.gz|example_data/SRR5201677_2.fastq.gz|myodes_db/mitochondrion/seed.fasta|myodes_db/mitochondrion/gene.fasta|
|M.glareolus_SRR5201676|example_data/SRR5201676_1.fastq.gz|example_data/SRR5201676_2.fastq.gz|myodes_db/mitochondrion/seed.fasta|myodes_db/mitochondrion/gene.fasta|


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
#SBATCH --mem=50G
#SBATCH --cpus-per-task=18

source activate genome_skimming_pipeline

snakemake \
   --cores 18 \
   --use-conda \
   --use-singularity \
   --configfile config/config_example.yaml

echo Complete!
```

