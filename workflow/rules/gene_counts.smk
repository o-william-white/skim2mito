rule gene_counts:
    input:
        get_mafft_filtered_output,
    output:
        "results/summary/summary_gene_counts_mqc.txt",
    log:
        "logs/gene_counts/gene_counts.log",
    conda:
        "../envs/conda_env.yaml"
    shell:
        """
        python3 workflow/scripts/gene_counts.py --input results/mafft_filtered/ --output {output}  &> {log}
        """
