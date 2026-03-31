rule mafft:
    input:
        f"{results_dir}/annotated_genes/{{dataset}}.fasta",
    output:
        f"{results_dir}/mafft/{{dataset}}.fasta",
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
