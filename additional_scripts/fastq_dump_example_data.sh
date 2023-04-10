#!/bin/bash

source activate sra-tools

for SRA in SRR5201684 SRR5201683 SRR5201682 SRR5201681 SRR5201680 SRR5201679 SRR5201678 SRR5201677 SRR5201676; do
   echo fetching reads for $SRA
   fastq-dump --origfmt --split-files -X 250000 --gzip --outdir example_data $SRA
done

echo Complete!

