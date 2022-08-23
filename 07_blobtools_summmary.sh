
echo -e "sample\tindex\tidentifiers\tgc\tlength\tmapped_cov\tmapped_read_cov\tsuperkingdom\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies" > blobtools/blobtools_summary.txt

cut -f 1 data_metadata.txt | while read REF; do

   ASM=$(echo selected_assemblies/${REF}.fasta)

   if [ -e ${ASM} ]; then 
   
      tail -n +2 blobtools/${REF}/blobtools/table.tsv | while read LINE; do
         echo -e "$REF\t$LINE"
      done
   
   fi 

done >> blobtools/blobtools_summary.txt

echo Complete!

