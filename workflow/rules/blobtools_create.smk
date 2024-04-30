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
        "results/assembled_sequence/{sample}.ok",
        "results/blastn/{sample}.ok",
        "results/minimap/{sample}.ok",
    output:
        ok="results/blobtools/{sample}/{sample}_create.ok",
    log:
        "logs/blobtools/{sample}_create.log",
    conda:
        "../envs/blobtools.yaml"
    resources:
        mem_mb=7164
    shell:
        """
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
        """
