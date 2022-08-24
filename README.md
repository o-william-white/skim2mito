
Set up input data
```
# get data
cp -r /workspaces/groups/clarkgroup/oliw/genome_skimming_simulations_2/simulation_benchmark/data .

# mv metadata to current dir
mv data/data_metadata.txt .
```

Set up conda env
```
# create env
conda create --name genome_skim_pipeline
conda activate genome_skim_pipeline

# python script to download references
conda install -c conda-forge biopython

# getorganelle
conda install -c bioconda getorganelle

# blobtools
conda install -c bioconda blast
conda install -c bioconda minimap2
conda install -c bioconda samtools
```

Blobtools installation
```
# blobtoolkit singularitiy and taxdump installed in home
mkdir /home/oliw/software/blobtools2/
cd /home/oliw/software/blobtools2/

# pull docker image
singularity pull docker://genomehubs/blobtoolkit

# download taxadump
cd /home/oliw/software/
mkdir taxdump
cd taxdump
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -
cd
```

Mitos and chloe installation
```
# mitos
# need to create a separate env for mitos
conda create -n mitos mitos=2.0.8 r-base=4.1

# download mitos referece database
cd /home/oliw/software/miniconda3/envs/mitos/bin
wget https://zenodo.org/record/2672835/files/refseq39.tar.bz2
wget https://zenodo.org/record/2672835/files/refseq63f.tar.bz2
wget https://zenodo.org/record/2672835/files/refseq63m.tar.bz2
wget https://zenodo.org/record/2672835/files/refseq63o.tar.bz2
wget https://zenodo.org/record/2672835/files/refseq81f.tar.bz2
wget https://zenodo.org/record/2672835/files/refseq81m.tar.bz2
wget https://zenodo.org/record/2672835/files/refseq81o.tar.bz2

tar xf refseq39.tar.bz2
tar xf refseq63f.tar.bz2
tar xf refseq63m.tar.bz2
tar xf refseq63o.tar.bz2
tar xf refseq81f.tar.bz2
tar xf refseq81m.tar.bz2
tar xf refseq81o.tar.bz2

# download additional libraries for the chloe api
conda activate mitos
conda install requests 
conda install pandas
conda deactivate
```

Main code
```
# count the number of organelle genomes available across taxa/lineages
sbatch 00_organelle_search_count.sh

# organelle assembles present for all, some have many (>100)
cut -f 1 data_metadata.txt | while read TAX; do  head -n 1 organelle_search_count/${TAX}_count.txt ; done

# download up to 5 references per taxon
sbatch 01_organelle_search_download.sh

# extract coding regions from genbank files
sbatch 02_get_annotated_regions_from_gb.sh

# check the number of annotated regions
grep -e "^>" -c organelle_search_download/*/annotated_regions/gene/gene.fasta

# fastp
sbatch sbatch 03_fastp.sh

# get organelle
sbatch 04_get_organelle.sh

# select first assembly option
bash 05_selected_assemblies.sh

# check for contamination with blobtools
sbatch 06_blobtools.sh
bash 07_blobtools_summmary.sh

# split assemblies containing multiple contigs
bash 08_split_fasta.sh

# annotate assemblies
09_annotation.sh

```


