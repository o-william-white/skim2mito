#!/bin/bash

mkdir -p mitos2_reference_data && cd mitos2_reference_data

wget https://zenodo.org/record/4284483/files/refseq39.tar.bz2
wget https://zenodo.org/record/4284483/files/refseq63f.tar.bz2
wget https://zenodo.org/record/4284483/files/refseq63m.tar.bz2
wget https://zenodo.org/record/4284483/files/refseq63o.tar.bz2
wget https://zenodo.org/record/4284483/files/refseq89f.tar.bz2
wget https://zenodo.org/record/4284483/files/refseq89m.tar.bz2
wget https://zenodo.org/record/4284483/files/refseq89o.tar.bz2

tar -xf refseq39.tar.bz2
tar -xf refseq63f.tar.bz2
tar -xf refseq63m.tar.bz2
tar -xf refseq63o.tar.bz2
tar -xf refseq89f.tar.bz2
tar -xf refseq89m.tar.bz2
tar -xf refseq89o.tar.bz2

cd ..

echo Complete!

