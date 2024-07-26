rule blobtools:
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
        "results/blobtools/{sample}/table.tsv",
    log:
        "logs/blobtools/{sample}.log",
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
        blobtools filter \
            --table {output} \
            --table-fields gc,length,{wildcards.sample}_cov,bestsumorder_superkingdom,bestsumorder_kingdom,bestsumorder_phylum,bestsumorder_class,bestsumorder_order,bestsumorder_family,bestsumorder_species \
            results/blobtools/{wildcards.sample} &>> {log}
        """
