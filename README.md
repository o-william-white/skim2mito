# skim2mito (with modifications for museomics)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**skim2mito** is a snakemake pipeline for the batch assembly, annotation, and phylogenetic analysis of mitochondrial genomes from low-coverage genome skims. The pipeline was designed to work with sequence data from museum collections. This fork includes numerous QoL, configurability, and museomics enhancements. It also contains a supplementary `filter_mitogenes.py` script that can be run after the pipeline to extract and validate gene annotations to ensure high quality.

## Contents
 - [Setup](#setup)
 - [Example data](#example-data)
 - [Input](#input)
 - [Output](#output)
 - [Running your own data](#running-your-own-data)
 - [Citations](#citations)

## Setup

The pipeline is written in Snakemake and uses conda to install the necessary tools.

It is *strongly recommended* to install conda using Mambaforge. See details here https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

Once conda is installed, you can pull the github repo and set up the base conda environment.

```bash
# get github repo
git clone https://github.com/SchistoDan/skim2mito

# change dir
cd skim2mito

# setup conda env
conda env create -n skim2mito -f skim2mito.yaml
```

<br/>
<div align="right">
    <b><a href="#skim2mito">↥ back to top</a></b>
</div>
<br/>

## Example data

Before you run your own data, it is recommended to run the example datasets provided. This will confirm there are no user-specific issues with the setup and it also installs all the dependencies. The example data includes simulated mitochondrial data from 25 different butterfly species.

To run the example data, use the code below. The first time you run the pipeline, it will take some time to install each of the conda environments, so it is a good time to take a tea break :).
```bash
conda activate skim2mito

snakemake --profile workflow/profiles/test
```

The resources requested for the test data can be found in the snakemake profile `workflow/profiles/test/config.yaml`. This requires 16G of memory and 4 cores, with the analysis taking approximately three hours. I have purposely tried to keep the memory requirements low so it can be run on most systems. However, the runtime can be improved by adjusting the profile to allow more memory and cores.

## Input

Snakemake requires a `config.yaml` and `samples.csv` to define input parameters and sequence data for each sample. 

The config file is located here `config/config.yaml`

The example samples.csv file is located here `config/samples.csv` and it looks like this (note that the seed and gene columns are only required if the custom getorganelle database option is specified in the config file):


 ID | forward | reverse | taxid | seed | gene 
----|---------|---------|-------|------|------
Adelpha_iphiclus | .test/reads/Adelpha_iphiclus_1.fq.gz | .test/reads/Adelpha_iphiclus_2.fq.gz | 100750 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Anartia_jatrophae_saturata | .test/reads/Anartia_jatrophae_saturata_1.fq.gz | .test/reads/Anartia_jatrophae_saturata_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Araschnia_levana | .test/reads/Araschnia_levana_1.fq.gz | .test/reads/Araschnia_levana_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Auzakia_danava | .test/reads/Auzakia_danava_1.fq.gz | .test/reads/Auzakia_danava_2.fq.gz | 100750 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Baeotus_beotus | .test/reads/Baeotus_beotus_1.fq.gz | .test/reads/Baeotus_beotus_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Catacroptera_cloanthe | .test/reads/Catacroptera_cloanthe_1.fq.gz | .test/reads/Catacroptera_cloanthe_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Chalinga_pratti | .test/reads/Chalinga_pratti_1.fq.gz | .test/reads/Chalinga_pratti_2.fq.gz | 100750 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Diaethria_gabaza_eupepla | .test/reads/Diaethria_gabaza_eupepla_1.fq.gz | .test/reads/Diaethria_gabaza_eupepla_2.fq.gz | 127268 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Doleschallia_melana | .test/reads/Doleschallia_melana_1.fq.gz | .test/reads/Doleschallia_melana_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Eurema_blanda | .test/reads/Eurema_blanda_1.fq.gz | .test/reads/Eurema_blanda_2.fq.gz | 42450 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Hypolimnas_usambara | .test/reads/Hypolimnas_usambara_1.fq.gz | .test/reads/Hypolimnas_usambara_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Junonia_villida | .test/reads/Junonia_villida_1.fq.gz | .test/reads/Junonia_villida_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Kallima_paralekta | .test/reads/Kallima_paralekta_1.fq.gz | .test/reads/Kallima_paralekta_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Kallimoides_rumia | .test/reads/Kallimoides_rumia_1.fq.gz | .test/reads/Kallimoides_rumia_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Litinga_cottini | .test/reads/Litinga_cottini_1.fq.gz | .test/reads/Litinga_cottini_2.fq.gz | 100750 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Mallika_jacksoni | .test/reads/Mallika_jacksoni_1.fq.gz | .test/reads/Mallika_jacksoni_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Moduza_procris | .test/reads/Moduza_procris_1.fq.gz | .test/reads/Moduza_procris_2.fq.gz | 100750 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Parasarpa_zayla | .test/reads/Parasarpa_zayla_1.fq.gz | .test/reads/Parasarpa_zayla_2.fq.gz | 100750 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Phaedyma_columella | .test/reads/Phaedyma_columella_1.fq.gz | .test/reads/Phaedyma_columella_2.fq.gz | 100750 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Precis_pelarga | .test/reads/Precis_pelarga_1.fq.gz | .test/reads/Precis_pelarga_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Protogoniomorpha_temora | .test/reads/Protogoniomorpha_temora_1.fq.gz | .test/reads/Protogoniomorpha_temora_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Salamis_cacta | .test/reads/Salamis_cacta_1.fq.gz | .test/reads/Salamis_cacta_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Smyrna_blomfildia | .test/reads/Smyrna_blomfildia_1.fq.gz | .test/reads/Smyrna_blomfildia_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Tacola_larymna | .test/reads/Tacola_larymna_1.fq.gz | .test/reads/Tacola_larymna_2.fq.gz | 100750 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta
Yoma_algina | .test/reads/Yoma_algina_1.fq.gz | .test/reads/Yoma_algina_2.fq.gz | 40040 | .test/seed_mitochondrion.fasta | .test/gene_mitochondrion.fasta


<br/>
<div align="right">
    <b><a href="#skim2mito">↥ back to top</a></b>
</div>
<br/>

## Output

All output files are saved to the `results` directory. Below is a table summarising all of the output files generated by the pipeline.

| Directory             | Description               |
|-----------------------|---------------------------|
| fastqc_raw            | Fastqc reports for raw input reads |
| fastp                 | Fastp reports from quality control of raw reads |
| fastqc_qc             | Fastqc reports for quality controlled reads |
| go_fetch              | Optional output containing reference databasesused by GetOrganelle |
| getorganelle          | GetOrganelle output with a directory for each sample |
| assembled_sequence    | Assembled sequences selected from GetOrganelle output and renamed |
| seqkit                | Seqkit summary of each assembly |
| blastn                | Blastn output of each assembly |
| minimap               | Mapping output of quality filtered reads against each assembly |
| blobtools             | Blobtools assembly summary collating blastn and mapping output |
| assess_assembly       | Plots of annotations, mean depth, GC content and proportion mismatches |
| annotations           | Annotation outputs of mitos |
| summary               | Summary per sample (seqkit stats), contig (GC content, length, coverage, taxonomy and annotations) and annotated gene counts |
| annotated_genes  | Unaligned fasta files of annotated genes identified across all samples |
| mafft                 | Mafft aligned fasta files of annotated genes identified across all samples |
| mafft_filtered        | Mafft aligned fasta files after the removal of sequences based on a missing data threshold |
| alignment_trim        | Ambiguous parts of alignment removed using either gblocks or clipkit |
| iqtree                | Iqtree phylogenetic analysis of annotated genes |
| plot_tree             | Plots of phylogenetic trees |
| multiqc               | Multiqc summary report |

<br/>
<div align="right">
    <b><a href="#skim2mito">↥ back to top</a></b>
</div>
<br/>


## Running your own data

The first thing you need to do is generate your own config.yaml and samples.csv files, using the files provided as a template.

GetOrganelle requires reference data in the format of seed and gene reference fasta files. You can generate these files yourself, or use go_fetch.py to download and format reference data formatted for GetOrganelle. go_fetch.py works by searching NCBI based on the NCBI taxonomy specified by the taxid column in the samples.csv file. Note that the seed and gene columns in the samples.csv file are only required if you want to provide your own custom GetOrganelle seed and gene reference databases. 

You can use the default reference data for GetOrganelle, but I would recommend using custom reference databases where possible. See here for details of how to set up your own databases https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#how-to-assemble-a-target-organelle-genome-using-my-own-reference


## Citations

If you use the pipeline, please cite the paper: 10.1111/1755-0998.14036

Since the pipeline is a wrapper for several other bioinformatic tools we also ask that you cite the tools used by the pipeline:
 - Fastqc https://github.com/s-andrews/FastQC
 - Fastp https://doi.org/10.1093/bioinformatics/bty560
 - GetOrganelle https://doi.org/10.1186/s13059-020-02154-5
 - Blastn https://doi.org/10.1186/1471-2105-10-421
 - Minimap2 https://doi.org/10.1093/bioinformatics/bty191
 - Blobtools https://doi.org/10.12688/f1000research.12232.1
 - Seqkit https://doi.org/10.1371/journal.pone.0163962
 - MITOS2 https://doi.org/10.1016/j.ympev.2012.08.023
 - Gblocks (default) https://doi.org/10.1093/oxfordjournals.molbev.a026334
 - Clipkit (optional) https://doi.org/10.1371/journal.pbio.3001007
 - Mafft https://doi.org/10.1093/molbev/mst010
 - Iqtree https://doi.org/10.1093/molbev/msu300
 - ete3 https://doi.org/10.1093/molbev/msw046
 - ggtree https://doi.org/10.1111/2041-210X.12628
 - Multiqc https://doi.org/10.1093/bioinformatics/btw354

<br/>
<div align="right">
    <b><a href="#skim2mito">↥ back to top</a></b>
</div>
<br/>
