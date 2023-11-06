#!/bin/bash

source activate skim2phylo

snakemake \
   --snakefile skim2phylo.smk \
   --cores 6 \
   --use-conda \
   --use-singularity \
   --configfile example_data/config_ribosomal.yaml \
   --rerun-incomplete

echo Complete!

