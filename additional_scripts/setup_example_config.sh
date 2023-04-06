#!/bin/bash

# setup example config.yaml
cat config/config.yaml |
   sed -e 's:config/samples.csv:config/samples_example.csv:g' \
       -e 's:results:results_example:g' \
       -e 's:<target_type>:animal_mt:g' \
       -e 's:</path/to/blastdb>:/workspaces/groups/database/nt-2022-06-22_blastdb/nt:g' \
       -e 's:</path/to/new_taxdump/>:taxdump:g' \
       -e 's:</path/to/mitos/refseq>:mitos2_reference_data:g' \
       -e 's:<code>:5:g' \
       -e 's:<barrnap_kingdom>:NA:g' \
       -e 's:<threads>:6:g' > config/config_example.yaml

# setup example samples.csv
head -n 1 /home/oliw/software/genome_skimming_pipeline/config/samples.csv > config/samples_example.csv
for SRA in ERR5665185 ERR5665186 ERR5665187; do
   echo -e "A.mellifera_${SRA},example_data/${SRA}_1.fastq.gz,example_data/${SRA}_2.fastq.gz,apis_mellifera_db/mitochondrion/seed.fasta,apis_mellifera_db/mitochondrion/gene.fasta"
done >> config/samples_example.csv

echo Complete!

