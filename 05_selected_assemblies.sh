mkdir -p selected_assemblies 

# get assembled datasets, edit sequence names and copy to new dir
cut -f 1 data_metadata.txt | while read SAMPLE; do
   if [ $(ls -1 get_organelle/${SAMPLE}/*_sequence.fasta 2>/dev/null | wc -l) -ge 1 ]; then
      ASM=$(ls -1 get_organelle/${SAMPLE}/*_sequence.fasta | head -n 1)
      sed -e "s/>/>${SAMPLE}_/g" $ASM > selected_assemblies/${SAMPLE}.fasta
   else
      echo Assembly failed for ${SAMPLE}
   fi
done

