# skim2phylo  

**skim2phylo** is a snakemake pipeline for the batch assembly, annotation, and phylogenetic analysis of mitochondrial genomes and ribosomal genes from low coverage "genome skims". The pipeline was designed specifically to work with sequence data from museum collections. 

An additional pipeline **gene2phylo** is provided which can be used to implement phylogenetic analysis for a given set of aligned genes.  

## Setup 

The pipeline is written in Snakemake and uses conda environments and singularity images to install the necessary tools. 

It is reccomended to install conda using Mambaforge. See details here https://snakemake.readthedocs.io/en/stable/getting_started/installation.html 

Once conda is installed, you can pull the github repo and set up the base conda environment.

```
# get github repo
git clone https://github.com/o-william-white/skim2phylo

# change dir
cd skim2phylo

# setup conda env
conda env create -n skim2phylo -f envs/conda_env.yaml
```

## Example usage

Before you run skim2phylo on your own data, it is recommended to run at least one of the example datasets provided. This will confirm there are no user-specific issues with the setup and it also installs all the dependencies. 

The example data includes simulated mitochondrial and ribosomal reads from 25 different butterfly speices. 

### Mitochondrial example

To run the mitochondrial example data, run the code below. *Note that the first time you run this, it will take some time to install each of the conda environments*, so it is a good time to take a tea break :).
```
source activate skim2phylo

snakemake \
   --snakefile skim2phylo.smk \
   --cores 6 \
   --use-conda \
   --use-singularity \
   --configfile example_data/config_mitochondrion.yaml \
   --rerun-incomplete
``` 

If you have access to a HPC, you can submit snakemake workflows as jobs. For example, using sbatch:
```
sbatch --cpus-per-task=24 --mem=16G run_mitochondrion_example.sh
```

Snakemake requires a config.yaml and samples.csv to define input paramters and sequence data for each sample. For the mitochondrion example data provided, the config file is located here `example_data/config_mitochondrion.yaml` and it looks like this: 
```
# path to sample sheet csv with columns for ID,forward,reverse
samples: example_data/samples_mitochondrion.csv

# name of output directory
output_dir: results_mitochondrion_example

# fastp depulication (True/False)
fastp_dedup: True

# GetOrganelle target sequence type (mitochondrion, ribosomal)
target_type: mitochondrion

# mitos refseq database (refseq39, refseq63f, refseq63m, refseq63o, refseq89f, refseq89m, refseq89o)
mitos_refseq: refseq39

# mito code (2 = Vertebrate, 4 = Mold, 5 = Invertebrate, 9 = Echinoderm, 13 = Ascidian, 14 = Alternative flatworm
mitos_code: 5

# barrnap kindgom (Bacteria:bac, Archaea:arc, Eukaryota:euk, Metazoan Mitochondria:mito)
barrnap_kingdom: NA

# alignment trimming method to use (gblocks or clipkit)
alignment_trim: gblocks

# alignment missing data threshold for alignment (0.0 - 1.0)
missing_threshold: 0.5

# name of outgroup sample (optional)
# use "NA" if there is no obvious outgroup
# if more than one outgroup use a comma separated list i.e. "sampleA,sampleB"
outgroup: "Eurema_blanda"

# plot dimensions (cm)
plot_height: 20
plot_width: 20

# number of threads to use
threads: 6
```

The samples.csv file is located here `example_data/samples_mitochondrion.csv` and the first 10 lines look like  this: 

| ID | forward | reverse | seed | gene |
|----|---------|---------|------|------|
| Adelpha_iphiclus           | example_data/mitochondrion/Adelpha_iphiclus_1.fq.gz           | example_data/mitochondrion/Adelpha_iphiclus_2.fq.gz           | example_data/seed_mitochondrion.fasta|example_data/gene_mitochondrion.fasta |
| Anartia_jatrophae_saturata | example_data/mitochondrion/Anartia_jatrophae_saturata_1.fq.gz | example_data/mitochondrion/Anartia_jatrophae_saturata_2.fq.gz | example_data/seed_mitochondrion.fasta|example_data/gene_mitochondrion.fasta |
| Araschnia_levana           | example_data/mitochondrion/Araschnia_levana_1.fq.gz           | example_data/mitochondrion/Araschnia_levana_2.fq.gz           | example_data/seed_mitochondrion.fasta|example_data/gene_mitochondrion.fasta |
| Auzakia_danava             | example_data/mitochondrion/Auzakia_danava_1.fq.gz             | example_data/mitochondrion/Auzakia_danava_2.fq.gz             | example_data/seed_mitochondrion.fasta|example_data/gene_mitochondrion.fasta |
| Baeotus_beotus             | example_data/mitochondrion/Baeotus_beotus_1.fq.gz             | example_data/mitochondrion/Baeotus_beotus_2.fq.gz             | example_data/seed_mitochondrion.fasta|example_data/gene_mitochondrion.fasta |
| Catacroptera_cloanthe      | example_data/mitochondrion/Catacroptera_cloanthe_1.fq.gz      | example_data/mitochondrion/Catacroptera_cloanthe_2.fq.gz      | example_data/seed_mitochondrion.fasta|example_data/gene_mitochondrion.fasta |
| Chalinga_pratti            | example_data/mitochondrion/Chalinga_pratti_1.fq.gz            | example_data/mitochondrion/Chalinga_pratti_2.fq.gz            | example_data/seed_mitochondrion.fasta|example_data/gene_mitochondrion.fasta |
| Diaethria_gabaza_eupepla   | example_data/mitochondrion/Diaethria_gabaza_eupepla_1.fq.gz   | example_data/mitochondrion/Diaethria_gabaza_eupepla_2.fq.gz   | example_data/seed_mitochondrion.fasta|example_data/gene_mitochondrion.fasta |
| Doleschallia_melana        | example_data/mitochondrion/Doleschallia_melana_1.fq.gz        | example_data/mitochondrion/Doleschallia_melana_2.fq.gz        | example_data/seed_mitochondrion.fasta|example_data/gene_mitochondrion.fasta |
| Eurema_blanda              | example_data/mitochondrion/Eurema_blanda_1.fq.gz              | example_data/mitochondrion/Eurema_blanda_2.fq.gz              | example_data/seed_mitochondrion.fasta|example_data/gene_mitochondrion.fasta |

### Ribosomal example

To run the ribosomal example data. run the code below. 
```
snakemake \
   --snakefile skim2phylo.smk \
   --cores 6 \
   --use-conda \
   --use-singularity \
   --configfile example_data/config_ribosomal.yaml \
   --rerun-incomplete
```
As above, the ribosomal example can be submitted using sbatch
```
sbatch --cpus-per-task=24 --mem=16G run_ribosomal_example.sh
```

### Main output files

TBC


### Filtering putative contaminants 

If you are working with museum collections, it is possible that you may assemble and annotate sequences from contaminant/non-target species. To remove putative contaminants, and generate the files required for a final phylogenetic analysis, a python script is provided: 
```
python additional_scripts/remove_contaminants.py \
   --input results_mitochondrion_example/mafft_filtered/ results_ribosomal_example/mafft_filtered/ \
   --cont Catacroptera_cloanthe 5.8S \
   --output results_genes_example \
   --overwrite
```

### Running you own data

You can generate your own config.yaml and samples.csv files. 

GetOrganelle requires reference data in the format of seed and gene reference fasta files. You can use the default reference data for GetOrganelle, but I would recomend using custom reference databases where possible. See here for details of how to set up your own databases https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#how-to-assemble-a-target-organelle-genome-using-my-own-reference 

I have shared a basic python script called go_fetch.py in another repo https://github.com/o-william-white/go_fetch to download and format reference data formatted for GetOrganelle. Go fetch downloads the reference from NCBI using biopython, removes repetitive sequences using trf, and formats the data for GetOrganelle.

