# skim2phylo  

**skim2phylo** is a snakemake pipeline for the batch assembly, annotation, and phylogenetic analysis of mitochondrial genomes and ribosomal genes from low coverage "genome skims". The pipeline was designed specifically to work with sequence data from museum collections. 

An additional pipeline **gene2phylo** is provided which can be used to implement further phylogenetic analysis for set of aligned genes, using individual genes, a partitioned gene alignment and astral.

## Contents
[Setup](#setup)

[Example usage](#example-usage)
   - [Mitochondrial example](#mitochondrial-example)
   - [Main output files](#main-output-files)
   - [Ribosomal example](#ribosomal-example)
   - [Counting annotated mitochondrial and ribosomal genes](#counting-annotated-mitochondrial-and-ribosomal-genes)
   - [Filtering putative contaminants and formatting alignments](#filtering-putative-contaminants-and-formatting-alignments)
   - [Re-running phylogenetic analysis](#re-running-phylogenetic-analysis)
[Barcode only](#barcode-only)

[Running your own data](#running-your-own-data)

<br/>
<div align="right">
    <b><a href="#skim2phylo">↥ back to top</a></b>
</div>
<br/>

## Setup 

The pipeline is written in Snakemake and uses conda environments and singularity images to install the necessary tools. 

It is recommended to install conda using Mambaforge. See details here https://snakemake.readthedocs.io/en/stable/getting_started/installation.html 

Once conda is installed, you can pull the github repo and set up the base conda environment.

```
# get github repo
git clone https://github.com/o-william-white/skim2phylo

# change dir
cd skim2phylo

# setup conda env
conda env create -n skim2phylo -f envs/conda_env.yaml
```

<br/>
<div align="right">
    <b><a href="#skim2phylo">↥ back to top</a></b>
</div>
<br/>

## Example usage

Before you run skim2phylo on your own data, it is recommended to run the example datasets provided. This will confirm there are no user-specific issues with the setup and it also installs all the dependencies. 

The example data includes simulated mitochondrial and ribosomal reads from 25 different butterfly species.

<br/>
<div align="right">
    <b><a href="#skim2phylo">↥ back to top</a></b>
</div>
<br/>

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
sbatch --cpus-per-task=24 --mem=30G run_mitochondrion_example.sh
```

Snakemake requires a config.yaml and samples.csv to define input parameters and sequence data for each sample. For the mitochondrion example data provided, the config file is located here `example_data/config_mitochondrion.yaml` and it looks like this: 
```
# path to sample sheet csv with columns for ID,forward,reverse,seed,gene
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

The example samples.csv file is located here `example_data/samples_mitochondrion.csv` and the first 10 lines look like  this: 

| ID | forward | reverse | seed | gene |
|----|---------|---------|------|------|
| Adelpha_iphiclus | example_data/mitochondrion/Adelpha_iphiclus_1.fq.gz | example_data/mitochondrion/Adelpha_iphiclus_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Anartia_jatrophae_saturata | example_data/mitochondrion/Anartia_jatrophae_saturata_1.fq.gz | example_data/mitochondrion/Anartia_jatrophae_saturata_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Araschnia_levana | example_data/mitochondrion/Araschnia_levana_1.fq.gz | example_data/mitochondrion/Araschnia_levana_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Auzakia_danava | example_data/mitochondrion/Auzakia_danava_1.fq.gz | example_data/mitochondrion/Auzakia_danava_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Baeotus_beotus | example_data/mitochondrion/Baeotus_beotus_1.fq.gz | example_data/mitochondrion/Baeotus_beotus_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Catacroptera_cloanthe | example_data/mitochondrion/Catacroptera_cloanthe_1.fq.gz | example_data/mitochondrion/Catacroptera_cloanthe_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Chalinga_pratti | example_data/mitochondrion/Chalinga_pratti_1.fq.gz | example_data/mitochondrion/Chalinga_pratti_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Diaethria_gabaza_eupepla | example_data/mitochondrion/Diaethria_gabaza_eupepla_1.fq.gz | example_data/mitochondrion/Diaethria_gabaza_eupepla_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Doleschallia_melana | example_data/mitochondrion/Doleschallia_melana_1.fq.gz | example_data/mitochondrion/Doleschallia_melana_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Eurema_blanda | example_data/mitochondrion/Eurema_blanda_1.fq.gz | example_data/mitochondrion/Eurema_blanda_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Hypolimnas_usambara | example_data/mitochondrion/Hypolimnas_usambara_1.fq.gz | example_data/mitochondrion/Hypolimnas_usambara_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Junonia_villida | example_data/mitochondrion/Junonia_villida_1.fq.gz | example_data/mitochondrion/Junonia_villida_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Kallima_paralekta | example_data/mitochondrion/Kallima_paralekta_1.fq.gz | example_data/mitochondrion/Kallima_paralekta_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Kallimoides_rumia | example_data/mitochondrion/Kallimoides_rumia_1.fq.gz | example_data/mitochondrion/Kallimoides_rumia_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Litinga_cottini | example_data/mitochondrion/Litinga_cottini_1.fq.gz | example_data/mitochondrion/Litinga_cottini_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Mallika_jacksoni | example_data/mitochondrion/Mallika_jacksoni_1.fq.gz | example_data/mitochondrion/Mallika_jacksoni_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Moduza_procris | example_data/mitochondrion/Moduza_procris_1.fq.gz | example_data/mitochondrion/Moduza_procris_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Parasarpa_zayla | example_data/mitochondrion/Parasarpa_zayla_1.fq.gz | example_data/mitochondrion/Parasarpa_zayla_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Phaedyma_columella | example_data/mitochondrion/Phaedyma_columella_1.fq.gz | example_data/mitochondrion/Phaedyma_columella_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Precis_pelarga | example_data/mitochondrion/Precis_pelarga_1.fq.gz | example_data/mitochondrion/Precis_pelarga_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Protogoniomorpha_temora | example_data/mitochondrion/Protogoniomorpha_temora_1.fq.gz | example_data/mitochondrion/Protogoniomorpha_temora_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Salamis_cacta | example_data/mitochondrion/Salamis_cacta_1.fq.gz | example_data/mitochondrion/Salamis_cacta_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Smyrna_blomfildia | example_data/mitochondrion/Smyrna_blomfildia_1.fq.gz | example_data/mitochondrion/Smyrna_blomfildia_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Tacola_larymna | example_data/mitochondrion/Tacola_larymna_1.fq.gz | example_data/mitochondrion/Tacola_larymna_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |
| Yoma_algina | example_data/mitochondrion/Yoma_algina_1.fq.gz | example_data/mitochondrion/Yoma_algina_2.fq.gz | example_data/seed_mitochondrion.fasta | example_data/gene_mitochondrion.fasta |

<br/>
<div align="right">
    <b><a href="#skim2phylo">↥ back to top</a></b>
</div>
<br/>

### Main output files

Below is a table summarising all of the output files generated by the pipeline.

| Directory             | Description                                                                                            |
|-----------------------|--------------------------------------------------------------------------------------------------------|
| fastqc                | Fastqc reports for input reads                                                                         |
| fastp                 | Fastp reports for input reads                                                                          |
| getorganelle          | GetOrganelle output with a directory for each sample                                                   |
| assembled_sequence    | Assembled sequences selected from GetOrganelle output and renamed                                      |
| seqkit                | Seqkit summary of the assembly                                                                         |
| blastn                | Blastn output of each assembly                                                                         |
| minimap               | Mapping output of quality filtered reads against each assembly                                         |
| blobtools             | Blobtools assembly summary collating blastn and mapping output                                         |
| assess_assembly       | Plots of annotations, mean depth, GC content and proportion mismatches                                 |
| annotations           | Annotation outputs of mitos or barrnap depending on target type                                        |
| summary               | Summary per sample (seqkit stats) and contig (GC content, length, coverage, taxonomy and annotations)  |
| protein_coding_genes  | Unaligned fasta files of annotated genes identified across all samples                                 |
| mafft                 | Mafft aligned fasta files of annotated genes identified across all samples                             |
| mafft_filtered        | Mafft aligned fasta files after the removal of sequences based on a missing data threshold             |
| alignment_trim        | Ambiguous parts of alignment removed using either gblocks or clipkit                                   |
| iqtree                | Iqtree phylogenetic analysis of annotated genes                                                        |
| plot_tree             | Plots of phylogenetic trees                                                                            |
| logs                  | Log files for each step in pipeline                                                                    |

<br/>
<div align="right">
    <b><a href="#skim2phylo">↥ back to top</a></b>
</div>
<br/>

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
sbatch --cpus-per-task=24 --mem=30G run_ribosomal_example.sh
```

### Counting annotated mitochondrial and ribosomal genes

Now we have assembled and annotated mitochondrial and ribosomal genes, it is useful to count the number of annotated genes across samples. This can be useful for identifying samples or genes with large amounts of missing data for downstream analyses.  

```
python scripts/gene_counts.py \
   --input results_mitochondrion_example/mafft_filtered/ results_ribosomal_example/mafft_filtered/ \
   --output gene_counts.txt
```

The output `gene_counts.txt` will look like this: 

| sample | 18S_rRNA | 28S_rRNA | 5_8S_rRNA | atp6 | atp8 | cob | cox1 | cox2 | cox3 | nad1 | nad2 | nad3 | nad4 | nad4l | nad5 | nad6 | rrnL | rrnS |
|--------|----------|----------|-----------|------|------|-----|------|------|------|------|------|------|------|-------|------|------|------|------|
| Adelpha_iphiclus | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Anartia_jatrophae_saturata | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Araschnia_levana | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Auzakia_danava | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 0 |
| Baeotus_beotus | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Catacroptera_cloanthe | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 |
| Chalinga_pratti | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Diaethria_gabaza_eupepla | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Doleschallia_melana | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Eurema_blanda | 0 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Hypolimnas_usambara | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Junonia_villida | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0 | 1 | 0 |
| Kallima_paralekta | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Kallimoides_rumia | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Litinga_cottini | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Mallika_jacksoni | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Moduza_procris | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Parasarpa_zayla | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Phaedyma_columella | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Precis_pelarga | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Protogoniomorpha_temora | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Salamis_cacta | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Smyrna_blomfildia | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Tacola_larymna | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |
| Yoma_algina | 1 | 0 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 |

<br/>
<div align="right">
    <b><a href="#skim2phylo">↥ back to top</a></b>
</div>
<br/>

### Filtering putative contaminants and formatting alignments 

If you are working with museum collections, it is possible that you may assemble and annotate sequences from contaminant/non-target species. *Contaminant sequences can be identified based on the blast search output or unusual placement in the phylogenetic trees* (see blobtools and plot_tree outputs). A python script `format_alignments.py `is provided to remove putative contaminants from alignments, and format the alignments required for a final phylogenetic analysis. 

For example, let's say we wanted to remove all sequences from the sample "Kallima_paralekta" and 5.8S ribosomal sequence, you could run the script as shown below. The script works by identifying and removing sequences that have names with  `Kallima_paralekta` or `5_8S` in the sequence names. The filtered alignments are written to a new output directory `filter_alignments_output`. 

```
python scripts/format_alignments.py  \
   --input results_mitochondrion_example/mafft_filtered/ results_ribosomal_example/mafft_filtered/ \
   --cont Kallima_paralekta 5_8S \
   --output filter_alignments_output \
   --overwrite
```

If there were no obvious contaminants but you would like to run more detailed phylogenetic analysis, you can run the same script without the `--cont` parameter.
```
python scripts/format_alignments.py \
   --input results_mitochondrion_example/mafft_filtered/ results_ribosomal_example/mafft_filtered/ \
   --output filter_alignments_output \
   --overwrite
```

*Note that the output fasta files have been reformatted so each alignment file is named after the gene and each sequence is named after the sample.* This is necessary to ensure that sequences are concatenated correctly downstream.  

<br/>
<div align="right">
    <b><a href="#skim2phylo">↥ back to top</a></b>
</div>
<br/>

### Re-running phylogenetic analysis

You should now have a single directory containing a combination of mitochondrial and ribosomal genes across multiple samples. The next optional step is to repeat the phylogenetic analysis using the supplementary pipeline **gene2phylo**, which will generate individual gene trees and a partitioned gene tree using iqtree and a species tree using astral. To run the example data provided, run the code below.
```
snakemake \
   --snakefile gene2phylo.smk \
   --cores 4 \
   --use-conda \
   --configfile example_data/config_genes.yaml \
   --rerun-incomplete
```
As above, the example can be submitted using sbatch
```
sbatch --cpus-per-task=24 --mem=10G run_genes_example.sh
```

<br/>
<div align="right">
    <b><a href="#skim2phylo">↥ back to top</a></b>
</div>
<br/>

## Barcode only 

If you are only interested in the assembly of genes and potentially barcode sequences without the phylogenetic analysis, you can stop the pipeline from running the gene alignment and other downstream steps using the `--omit-from` parameter. 
```
snakemake \
   --snakefile skim2phylo.smk \
   --cores 6 \
   --use-conda \
   --use-singularity \
   --configfile example_data/config_mitochondrion.yaml \
   --rerun-incomplete \
   --omit-from mafft
```

<br/>
<div align="right">
    <b><a href="#skim2phylo">↥ back to top</a></b>
</div>
<br/>
  
## Running your own data

The first thing you need to do is generate your own config.yaml and samples.csv files. 

GetOrganelle requires reference data in the format of seed and gene reference fasta files. You can use the default reference data for GetOrganelle, but I would recommend using custom reference databases where possible. See here for details of how to set up your own databases https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#how-to-assemble-a-target-organelle-genome-using-my-own-reference 

I have shared a basic python script called go_fetch.py in another repo https://github.com/o-william-white/go_fetch to download and format reference data formatted for GetOrganelle. Go fetch downloads the reference from NCBI using biopython, removes repetitive sequences using trf, and formats the data for GetOrganelle.

If you have any questions, please do get in touch in the issues or by email o.william.white@gmail.com

