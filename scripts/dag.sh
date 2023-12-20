#!/bin/bash

source activate skim2phylo

snakemake \
   --snakefile skim2phylo.smk \
   --configfile example_data/config_mitochondrion.yaml \
   --dag | dot -Tsvg > dag_skim2phylo.svg

snakemake \
   --snakefile gene2phylo.smk \
   --configfile example_data/config_genes.yaml \
   --dag | dot -Tsvg > dag_gene2phylo.svg

echo Complete!

