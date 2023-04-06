#!/bin/bash

source activate sra-tools

for SRA in ERR5665185 ERR5665186 ERR5665187 ERR5665184 ERR5665183 ERR5665182 ERR5665181 ERR5665180; do 
   echo fetching reads for $SRA
   fastq-dump --origfmt --split-files -X 500000 --gzip --outdir example_data $SRA
done

echo Complete!

