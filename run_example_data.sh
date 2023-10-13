#!/bin/bash

source activate skim2phylo

snakemake \
   --cores 4 \
   --use-conda \
   --use-singularity \
   --configfile example_data/config_example.yaml \
   --rerun-incomplete

echo Complete!

