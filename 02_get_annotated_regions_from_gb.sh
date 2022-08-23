#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_get_annotated_regions_from_gb_%j_%a.out
#SBATCH --error=job_get_annotated_regions_from_gb_%j_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-26

export PATH=/home/oliw/software/miniconda3/bin/:$PATH
source activate genome_skim_pipeline

REF=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 1)
TAX=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 2)
ORG=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 3)
DAT=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 4)
LAY=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 5)

# see https://github.com/Kinggerm/GetOrganelle/wiki/FAQ#how-to-assemble-a-target-organelle-genome-using-my-own-reference for more details
python additional_scripts/get_annotated_regions_from_gb.py \
  $(echo organelle_search_download/${REF}/genbank/*.gb) \
  -o organelle_search_download/${REF}/annotated_regions \
  -t CDS \
  --mix

echo Complete!

