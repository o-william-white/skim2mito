#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_get_organelle_%j_%a.out
#SBATCH --error=job_get_organelle_%j_%a.err
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-26

export PATH=/home/oliw/software/miniconda3/bin/:$PATH
source activate getorganelle

REF=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 1)
TAX=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 2)
ORG=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 3)
DAT=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 4)
LAY=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 5)

mkdir -p get_organelle
if [ $LAY == "paired" ]; then
  
  get_organelle_from_reads.py \
    -1 fastp/${REF}/${REF}_R1.fastq.gz \
    -2 fastp/${REF}/${REF}_R2.fastq.gz \
    -t 1 \
    -o get_organelle/${REF} \
    -F ${DAT} \
    -s organelle_search_download/${REF}/annotated_regions/gene/gene.fasta \
    --genes organelle_search_download/${REF}/annotated_regions/gene/gene.fasta \
    --reduce-reads-for-coverage inf --max-reads inf \
    --overwrite
fi

if [ $LAY == "unpaired" ]; then
  get_organelle_from_reads.py \
    -u fastp/${REF}/${REF}.fastq.gz \
    -t 1 \
    -o get_organelle/${REF} \
    -F ${DAT} \
    -s organelle_search_download/${REF}/annotated_regions/gene/gene.fasta \
    --genes organelle_search_download/${REF}/annotated_regions/gene/gene.fasta \
    --reduce-reads-for-coverage inf --max-reads inf \
    --overwrite
fi

