rule blobtools_filter:
    input:
        "results/blobtools/{sample}/bestsumorder_class_cindex.json",
        "results/blobtools/{sample}/bestsumorder_class.json",
        "results/blobtools/{sample}/bestsumorder_class_positions.json",
        "results/blobtools/{sample}/bestsumorder_class_score.json",
        "results/blobtools/{sample}/bestsumorder_family_cindex.json",
        "results/blobtools/{sample}/bestsumorder_family.json",
        "results/blobtools/{sample}/bestsumorder_family_positions.json",
        "results/blobtools/{sample}/bestsumorder_family_score.json",
        "results/blobtools/{sample}/bestsumorder_genus_cindex.json",
        "results/blobtools/{sample}/bestsumorder_genus.json",
        "results/blobtools/{sample}/bestsumorder_genus_positions.json",
        "results/blobtools/{sample}/bestsumorder_genus_score.json",
        "results/blobtools/{sample}/bestsumorder_kingdom_cindex.json",
        "results/blobtools/{sample}/bestsumorder_kingdom.json",
        "results/blobtools/{sample}/bestsumorder_kingdom_positions.json",
        "results/blobtools/{sample}/bestsumorder_kingdom_score.json",
        "results/blobtools/{sample}/bestsumorder_order_cindex.json",
        "results/blobtools/{sample}/bestsumorder_order.json",
        "results/blobtools/{sample}/bestsumorder_order_positions.json",
        "results/blobtools/{sample}/bestsumorder_order_score.json",
        "results/blobtools/{sample}/bestsumorder_phylum_cindex.json",
        "results/blobtools/{sample}/bestsumorder_phylum.json",
        "results/blobtools/{sample}/bestsumorder_phylum_positions.json",
        "results/blobtools/{sample}/bestsumorder_phylum_score.json",
        "results/blobtools/{sample}/bestsumorder_positions.json",
        "results/blobtools/{sample}/bestsumorder_species_cindex.json",
        "results/blobtools/{sample}/bestsumorder_species.json",
        "results/blobtools/{sample}/bestsumorder_species_positions.json",
        "results/blobtools/{sample}/bestsumorder_species_score.json",
        "results/blobtools/{sample}/bestsumorder_superkingdom_cindex.json",
        "results/blobtools/{sample}/bestsumorder_superkingdom.json",
        "results/blobtools/{sample}/bestsumorder_superkingdom_positions.json",
        "results/blobtools/{sample}/bestsumorder_superkingdom_score.json",
        "results/blobtools/{sample}/{sample}_cov.json",
        "results/blobtools/{sample}/gc.json",
        "results/blobtools/{sample}/identifiers.json",
        "results/blobtools/{sample}/length.json",
        "results/blobtools/{sample}/meta.json",
        "results/blobtools/{sample}/ncount.json",
    output:
        "results/blobtools/{sample}/table.tsv",
    log:
        "logs/blobtools_filter/{sample}.log",
    conda:
        "../envs/blobtools.yaml"
    shell:
        """
        blobtools filter \
            --table {output} \
            --table-fields gc,length,{wildcards.sample}_cov,bestsumorder_superkingdom,bestsumorder_kingdom,bestsumorder_phylum,bestsumorder_class,bestsumorder_order,bestsumorder_family,bestsumorder_species \
            results/blobtools/{wildcards.sample} &> {log}
        """
