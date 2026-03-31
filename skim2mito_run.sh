#!/bin/bash
#SBATCH --job-name=S2M
#SBATCH --partition=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err



## Conda environment
source /mnt/apps/users/$USER/conda/etc/profile.d/conda.sh

conda activate skim2mito



## Snakemake run
# Unlock
snakemake --cores 8 --use-conda --unlock --profile workflow/profiles/slurm

# Run
snakemake --cores 8 --use-conda --profile workflow/profiles/slurm


