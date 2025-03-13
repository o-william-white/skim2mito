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
<<<<<<< HEAD:workflow/rules/blobtools_create.smk
        ok="results/blobtools/{sample}/{sample}_create.ok",
=======
        "results/blobtools/{sample}/table.tsv",
>>>>>>> main:workflow/rules/blobtools.smk
    log:
        "logs/blobtools/{sample}_create.log",
    conda:
        "../envs/blobtools.yaml"
    shell:
        """
<<<<<<< HEAD:workflow/rules/blobtools_create.smk
        FAS=$(echo results/assembled_sequence/{wildcards.sample}.fasta)
        BLA=$(echo results/blastn/{wildcards.sample}.txt)
        MAP=$(echo results/minimap/{wildcards.sample}.bam)
        OUT=$(echo results/blobtools/{wildcards.sample}/table.tsv)
        if [ -e $FAS ]; then
            blobtools create \
                --fasta $FAS \
                --hits $BLA \
                --taxrule bestsumorder \
                --taxdump resources/taxdump \
                --cov $MAP \
                results/blobtools/{wildcards.sample} &> {log}
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
=======
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
>>>>>>> main:workflow/rules/blobtools.smk
        """
