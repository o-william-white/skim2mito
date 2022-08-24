#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_annotation_%j_%a.out
#SBATCH --error=job_annotation_%j_%a.err
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-26%1

# chloe ape can only run one at a time

export PATH=/home/oliw/software/miniconda3/bin/:$PATH
source activate mitos

#export PATH=/home/oliw/software/miniconda3/bin/:$PATH
#source activate genome_skim_pipeline

REF=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 1)
TAX=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 2)
ORG=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 3)
DAT=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 4)
LAY=$(sed -n ${SLURM_ARRAY_TASK_ID}p data_metadata.txt | cut -f 5)

ASM=$(echo selected_assemblies/${REF}.fasta)

if [ -e ${ASM} ]; then

   mkdir -p annotation
   mkdir -p annotation/${REF}
   
   if [[ $ORG == "chloroplast" ]]; then
      echo Running chloe api
      ls -1 split_fasta/${REF}/*.fasta | while read FAS; do
         echo Annotating ${FAS}
	 DIR=$(basename ${FAS} | sed 's/.fasta//g')
	 mkdir -p annotation/${REF}/${DIR}
	 python additional_scripts/chloe_api.py \
            --input ${FAS} \
	    --output annotation/${REF}/${DIR}
      done
   fi
   
   if [[ $ORG == "mitochondrion" ]]; then
      echo Running mitos
      ls -1 split_fasta/${REF}/*.fasta | while read FAS; do
         echo Annotating ${FAS}
	 DIR=$(basename ${FAS} | sed 's/.fasta//g')
	 mkdir -p annotation/${REF}/${DIR}
         runmitos.py \
            --input ${FAS} \
            --code 5 \
            --outdir annotation/${REF}/${DIR} \
            --refseqver /home/oliw/software/miniconda3/envs/mitos/bin/refseq81m \
            --refdir . \
            --noplots
      done
   fi

else

   echo Assembly failed for ${REF} 

fi

echo Complete!

