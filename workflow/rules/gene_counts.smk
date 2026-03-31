rule gene_counts:
    input:
        get_mafft_filtered_output,
    output:
        f"{results_dir}/summary/summary_gene_counts_mqc.txt",
    params:
        results_dir=results_dir,
    log:
        "logs/gene_counts/gene_counts.log",
    conda:
        "../envs/conda_env.yaml"
    shell:
        """
        python3 workflow/scripts/gene_counts.py --input {params.results_dir}/mafft_filtered/ --output {output}  &> {log}
        """
