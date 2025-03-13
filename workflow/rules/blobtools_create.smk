rule blobtools_create:
    input:
        "resources/taxdump/",
        "resources/taxdump/citations.dmp",
        "resources/taxdump/delnodes.dmp",
        "resources/taxdump/division.dmp",
        "resources/taxdump/excludedfromtype.dmp",
        "resources/taxdump/fullnamelineage.dmp",
        "resources/taxdump/gencode.dmp",
        "resources/taxdump/host.dmp",
        "resources/taxdump/images.dmp",
        "resources/taxdump/merged.dmp",
        "resources/taxdump/names.dmp",
        "resources/taxdump/nodes.dmp",
        "resources/taxdump/rankedlineage.dmp",
        "resources/taxdump/taxidlineage.dmp",
        "resources/taxdump/typematerial.dmp",
        "resources/taxdump/typeoftype.dmp",
        fas="results/assembled_sequence/{sample}.fasta",
        bla="results/blastn/{sample}.txt",
        bam="results/minimap/{sample}.bam",
    output:
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
    log:
        "logs/blobtools_create/{sample}.log",
    conda:
        "../envs/blobtools.yaml"
    shell:
        """
        blobtools create \
            --fasta {input.fas} \
            --hits {input.bla} \
            --taxrule bestsumorder \
            --taxdump resources/taxdump \
            --cov {input.bam} \
            results/blobtools/{wildcards.sample} &> {log}
        """
