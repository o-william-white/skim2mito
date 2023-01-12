#!/bin/bash

mkdir -p taxdump && cd taxdump

wget https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

gunzip -c new_taxdump.tar.gz | tar xf -

cd ..

echo Complete!

