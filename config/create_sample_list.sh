#!/usr/bin/env bash

# script creates sample_list.csv which should be specified in config.yaml
# run it in the "config" folder

# accomodates:
# https://github.com/o-william-white/genome_skimming_pipeline

# assumption: 
# sequence data and reference genome db are located in "genome_skimming_pipeline"
# the first 14 characters of the sequence name can be used as ID (eg. RMNH_5117143_1, negativecontro)

# usage: ./create_sample_list.sh $1 $2 $3
# $1 = path to the folder that contains the R1 and R2 fastq.gz files (raw data)
# $2 = path to the seed file (reference genome db)
# $3 = path to the gene file (reference genome db)

# example: ./create_sample_list.sh \
# /data/dick.groenenberg/genome_skimming_pipeline/154286/raw_sequences \
# /data/dick.groenenberg/genome_skimming_pipeline/calliphora_db/mitochondrion/seed.fasta \
# /data/dick.groenenberg/genome_skimming_pipeline/calliphora_db/mitochondrion/gene.fasta

# outfile name
tempfile=tempfile.txt
tempfile2=tempfile2.txt
tempfile3=tempfile3.txt
outfile=sample_list.csv

# test if temp- and outfiles exist; exit if they do, create them if they don't exit
if [[ -e "$tempfile" || -e "$tempfile2" || -e "$outfile" ]]; then
	printf "script terminated:\n tempfile(s) or sample_list.csv already exist\n" && exit 1
else
	touch "$tempfile" "$tempfile2" "$outfile"
fi

# put forward and reverse in tempfile (assuming candidate files are selected by "RMNH" or "negative")
ls -1 "$1"/*.gz | awk -F"genome_skimming_pipeline/" '{print $2}' | egrep "RMNH|negative" | sed 'N;s/\(.*\)\n\(.*\)/\1,\2/' > "$tempfile"

# add seed and gene dbs (tempfile2)
seed=$(echo "$2" |  awk -F"genome_skimming_pipeline/" '{print $2}')
gene=$(echo "$2" |  awk -F"genome_skimming_pipeline/" '{print $2}')
while read line; do
    printf "$line,$seed,$gene\n" >> "$tempfile2"
done < $tempfile

# extract ID (tempfile3)
awk -F"raw_sequences/" '{print $2}' "$tempfile2"| cut -c 1-14 > "$tempfile3"

# merge tempfile3 and tempfile2
paste -d"," "$tempfile3" "$tempfile2" > "$outfile"

# add header
awk -i inplace 'BEGINFILE{print "ID,forward,reverse,seed,gene"}{print}' "$outfile"

# warn if ID (tempfile3) contains duplicates!
cat "$tempfile3" | sort -n | uniq -c | awk '{if ($1 > 1) print "\nDuplicate ID(s) detected!"; exit 1}'
cat "$tempfile3" | sort -n | uniq -c | awk '{if ($1 > 1) print $2}'

# remove tempfiles
rm "$tempfile" "$tempfile2" "$tempfile3"
