rule blobtools_filter:
    input:
        "results/blobtools/{sample}/{sample}_create.ok",
    output:
        ok="results/blobtools/{sample}/{sample}_filter.ok",
    log:
        "logs/blobtools/{sample}_filter.log",
    conda:
        "../envs/blobtools.yaml"
    shell:
        """
        FAS=$(echo results/assembled_sequence/{wildcards.sample}.fasta)
        BLA=$(echo results/blastn/{wildcards.sample}.txt)
        MAP=$(echo results/minimap/{wildcards.sample}.bam)
        OUT=$(echo results/blobtools/{wildcards.sample}/table.tsv)
        if [ -e $FAS ]; then
            blobtools filter \
                --table $OUT \
                --table-fields gc,length,{wildcards.sample}_cov,bestsumorder_superkingdom,bestsumorder_kingdom,bestsumorder_phylum,bestsumorder_class,bestsumorder_order,bestsumorder_family,bestsumorder_species \
                results/blobtools/{wildcards.sample} &> {log}
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """
