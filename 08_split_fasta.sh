cut -f 1 data_metadata.txt | while read REF; do
   
   mkdir -p split_fasta
   mkdir -p split_fasta/${REF}

   ASM=$(echo selected_assemblies/${REF}.fasta)

   if [ -e ${ASM} ]; then
      
      echo Splitting fasta for ${REF}
      python additional_scripts/split_fasta.py --input ${ASM} --output split_fasta/${REF} 

   else

      echo Assembly failed for ${REF}      

   fi

done

echo Complete!

