# Genome skimming pipeline

Snakemake pipeline to assemble organelle or ribsomal sequences from genome skimming data using GetOrganelle.

See repo for go_fetch https://github.com/o-william-white/go_fetch for downloading reference databases formated for GetOrganelle.

```
# setup conda env using mamba
mamba env create -n genome_skimming_pipeline -f genome_skimming_pipeline/envs/conda_env.yaml

# wget new_taxdump
bash additional_scripts/wget_new_taxdump.sh

# get blast database
# TBC

# get mitos refseq
# TBC

```

