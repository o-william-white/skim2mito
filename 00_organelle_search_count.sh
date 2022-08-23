#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_organelle_search_count_%j_%a.out
#SBATCH --error=job_organelle_search_count_%j_%a.err
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-26%1

export PATH=/home/oliw/software/miniconda3/bin/:$PATH
source activate genome_skim_pipeline

REF=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 1)
TAX=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 2)
ORG=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 3)
DAT=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 4)
LAY=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 5)

mkdir -p organelle_search_count

python additional_scripts/organelle_search_ncbi.py \
   --taxonomy "${TAX}" \
   --organelle ${ORG} \
   --email o.william.white@gmail.com \
   --action count > organelle_search_count/${REF}_count.txt

echo Complete!

