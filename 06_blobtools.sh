#!/bin/bash
#SBATCH --partition=day
#SBATCH --output=job_blobtools_%j_%a.out
#SBATCH --error=job_blobtools_%j_%a.err
#SBATCH --time=1:00:00
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

ASM=$(echo selected_assemblies/${REF}.fasta)

if [ -e ${ASM} ]; then
    
      mkdir -p blobtools
      mkdir -p blobtools/${REF}

      echo Running blast
      blastn \
         -query $ASM \
         -db /workspaces/groups/database/nt-2021-09-07/nt \
         -out blobtools/${REF}/blastn.txt \
         -outfmt "6 qseqid staxids bitscore std" \
         -max_target_seqs 10 \
         -max_hsps 1 \
         -evalue 1e-25 \
         -num_threads 16
 
      echo Running minimap
      if [ $LAY == "paired" ]; then
         minimap2 \
            -ax sr \
            -t 16 $ASM \
            fastp/${REF}/${REF}_R1.fastq.gz fastp/${REF}/${REF}_R1.fastq.gz \
            | samtools sort -@16 -O BAM -o blobtools/${REF}/mapped.bam -
      fi
      if [ $LAY == "unpaired" ]; then
         minimap2 \
            -ax sr \
            -t 16 $ASM \
            fastp/${REF}/${REF}.fastq.gz \
            | samtools sort -@16 -O BAM -o blobtools/${REF}/mapped.bam -
      fi
       

      echo Running blobtools
      # singularity only works in home 
      mkdir -p /home/oliw/software/blobtools2/tmp_${REF} # rm 'tmp_' if repeated
      cp $ASM /home/oliw/software/blobtools2/tmp_${REF}/ref_${REF}.fasta
      cp blobtools/${REF}/blastn.txt /home/oliw/software/blobtools2/tmp_${REF}/
      cp blobtools/${REF}/mapped.bam /home/oliw/software/blobtools2/tmp_${REF}/

      cd /home/oliw/software/blobtools2/

      singularity exec blobtoolkit_latest.sif blobtools create --fasta tmp_${REF}/ref_${REF}.fasta tmp_${REF}

      singularity exec blobtoolkit_latest.sif blobtools add \
         --hits tmp_${REF}/blastn.txt \
         --taxrule bestsumorder \
         --taxdump /home/oliw/software/taxdump \
         --cov tmp_${REF}/mapped.bam \
         tmp_${REF}

      singularity exec blobtoolkit_latest.sif blobtools filter \
         --table tmp_${REF}/table.tsv \
         --table-fields gc,length,mapped_cov,mapped_read_cov,bestsumorder_superkingdom,bestsumorder_kingdom,bestsumorder_phylum,bestsumorder_class,bestsumorder_order,bestsumorder_family,bestsumorder_species \
         tmp_${REF}

      mv tmp_${REF} /workspaces/groups/clarkgroup/oliw/genome_skimming_pipeline/blobtools/${REF}/blobtools

   else

      echo Assembly failed for ${REF}

fi

echo Complete!

