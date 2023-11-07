#!/bin/bash

source activate skim2phylo

snakemake \
   --snakefile gene2phylo.smk \
   --cores 4 \
   --use-conda \
   --configfile example_data/config_genes.yaml \
   --rerun-incomplete

echo Complete!

