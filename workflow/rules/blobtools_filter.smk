rule blobtools_filter:
    input:
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
    output:
        f"{results_dir}/blobtools/{{sample}}/table.tsv",
    params:
        results_dir=results_dir,
    log:
        "logs/blobtools_filter/{sample}.log",
    conda:
        "../envs/blobtools.yaml"
    shell:
        """
        blobtools filter \
            --table {output} \
            --table-fields gc,length,{wildcards.sample}_cov,bestsumorder_superkingdom,bestsumorder_kingdom,bestsumorder_phylum,bestsumorder_class,bestsumorder_order,bestsumorder_family,bestsumorder_species \
            {params.results_dir}/blobtools/{wildcards.sample} &> {log}
        """
