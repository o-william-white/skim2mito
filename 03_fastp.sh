#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_fastp_%j_%a.out
#SBATCH --error=job_fastp_%j_%a.err
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-26

export PATH=/home/oliw/software/miniconda3/bin/:$PATH
source activate fastp

REF=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 1)
TAX=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 2)
ORG=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 3)
DAT=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 4)
LAY=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 5)

mkdir -p fastp
mkdir -p fastp/${REF}

if [ $LAY == "paired" ]; then
  fastp \
    -i data/${REF}_R1.fastq.gz \
    -I data/${REF}_R2.fastq.gz \
    -o fastp/${REF}/${REF}_R1.fastq.gz \
    -O fastp/${REF}/${REF}_R2.fastq.gz \
    -j fastp/${REF}/${REF}.json \
    -h fastp/${REF}/${REF}.html
fi 


if [ $LAY == "unpaired" ]; then
  fastp \
    -i data/${REF}.fastq.gz \
    -o fastp/${REF}/${REF}.fastq.gz \
    -j fastp/${REF}/${REF}.json \
    -h fastp/${REF}/${REF}.html
fi

echo Complete!

