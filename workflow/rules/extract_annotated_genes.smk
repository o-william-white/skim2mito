checkpoint extract_annotated_genes:
    input:
        expand(
            "results/annotations/{sample}/{sample}.ok",
            sample=sample_data["ID"].tolist(),
        ),
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
