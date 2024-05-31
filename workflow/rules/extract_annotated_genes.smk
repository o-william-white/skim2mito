checkpoint extract_annotated_genes:
    input:
        get_annotated_samples,
    output:
        directory("results/annotated_genes/"),
    log:
        "logs/annotated_genes/annotated_genes.log",
    conda:
        "../envs/conda_env.yaml"
    shell:
        """
        python workflow/scripts/mitos_alignments.py results/annotations/ results/annotated_genes &> {log}
        """
