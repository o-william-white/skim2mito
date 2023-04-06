#!/bin/bash

source activate sra-tools

for SRA in ERR5665185 ERR5665186 ERR5665187; do 
   echo fetching reads for $SRA
   fastq-dump --origfmt --split-files -X 500000 --gzip --outdir example_data $SRA
done

echo Complete!

