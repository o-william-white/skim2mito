#!/bin/bash

# setup example config.yaml
cat config/config.yaml |
   sed -e 's:config/samples.csv:config/samples_example.csv:g' \
       -e 's:results:results_example:g' \
       -e 's:<target_type>:animal_mt:g' \
       -e 's:</path/to/blastdb>:/workspaces/groups/database/nt-2022-06-22_blastdb/nt:g' \
       -e 's:</path/to/new_taxdump/>:taxdump:g' \
       -e 's:</path/to/mitos/refseq>:mitos2_reference_data/refseq89m:g' \
       -e 's:<code>:2:g' \
       -e 's:<barrnap_kingdom>:NA:g' \
       -e 's:<threads>:6:g' > config/config_example.yaml

# setup example samples.csv
head -n 1 /home/oliw/software/genome_skimming_pipeline/config/samples.csv > config/samples_example.csv
for SRA in SRR5201684 SRR5201683 SRR5201682 SRR5201681 SRR5201680 SRR5201679 SRR5201678 SRR5201677 SRR5201676; do
   echo -e "M.glareolus_${SRA},example_data/${SRA}_1.fastq.gz,example_data/${SRA}_2.fastq.gz,myodes_db/mitochondrion/seed.fasta,myodes_db/mitochondrion/gene.fasta"
done >> config/samples_example.csv

echo Complete!

