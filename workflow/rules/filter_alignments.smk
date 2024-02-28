rule filter_alignments:
    input:
        "results/mafft/{dataset}.fasta"
    params:
        threshold = missing_threshold
    output:
        "results/mafft_filtered/{dataset}.fasta"
    log:
        "logs/mafft_filtered/{dataset}.log"
    conda:
        "../envs/conda_env.yaml"
    shell:
        """
        python workflow/scripts/alignments_filter.py --input {input} --output {output} --threshold {params.threshold} > {log}
        """
