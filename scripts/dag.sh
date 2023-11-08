#!/bin/bash

source activate skim2phylo

snakemake \
   --configfile example_data/config_mitochondrion.yaml \
   --dag results_mitochondrion_example/summary/summary_contig.txt | dot -Tpdf > dag_part1.pdf

snakemake \
   --configfile config/config_example.yaml \
   --dag results_example/snakemake.ok | dot -Tpdf > dag_part2.pdf

echo Complete!

