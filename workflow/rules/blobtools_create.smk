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
        fas=f"{results_dir}/assembled_sequence/{{sample}}.fasta",
        bla=f"{results_dir}/blastn/{{sample}}.txt",
        bam=f"{results_dir}/minimap/{{sample}}.bam",
    output:
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_class_cindex.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_class.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_class_positions.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_class_score.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_family_cindex.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_family.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_family_positions.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_family_score.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_genus_cindex.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_genus.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_genus_positions.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_genus_score.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_kingdom_cindex.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_kingdom.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_kingdom_positions.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_kingdom_score.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_order_cindex.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_order.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_order_positions.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_order_score.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_phylum_cindex.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_phylum.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_phylum_positions.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_phylum_score.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_positions.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_species_cindex.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_species.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_species_positions.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_species_score.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_superkingdom_cindex.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_superkingdom.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_superkingdom_positions.json",
        f"{results_dir}/blobtools/{{sample}}/bestsumorder_superkingdom_score.json",
        f"{results_dir}/blobtools/{{sample}}/{{sample}}_cov.json",
        f"{results_dir}/blobtools/{{sample}}/gc.json",
        f"{results_dir}/blobtools/{{sample}}/identifiers.json",
        f"{results_dir}/blobtools/{{sample}}/length.json",
        f"{results_dir}/blobtools/{{sample}}/meta.json",
        f"{results_dir}/blobtools/{{sample}}/ncount.json",
    params:
        results_dir=results_dir,
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
            {params.results_dir}/blobtools/{wildcards.sample} &> {log}
        """
