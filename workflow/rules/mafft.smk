rule mafft:
    input:
        "results/annotated_genes/{dataset}.fasta",
    output:
        "results/mafft/{dataset}.fasta",
    log:
        "logs/mafft/{dataset}.log",
    conda:
        "../envs/mafft.yaml"
    shell:
        """
        mafft \
            --maxiterate 1000 \
            --globalpair \
            --adjustdirectionaccurately \
            {input} 1> {output} 2> {log}
        """
