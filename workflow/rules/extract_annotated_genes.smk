checkpoint extract_annotated_genes:
    input:
        get_annotated_samples,
    output:
        directory(f"{results_dir}/annotated_genes/"),
        summary=f"{results_dir}/annotated_genes/summary.txt"
    params:
        results_dir=results_dir,
    log:
        "logs/annotated_genes/annotated_genes.log",
    conda:
        "../envs/conda_env.yaml"
    shell:
        """
        python workflow/scripts/mitos_alignments.py {params.results_dir}/annotations/ {params.results_dir}/annotated_genes &> {log}
        """
