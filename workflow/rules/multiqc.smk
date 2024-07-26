rule multiqc:
    input:
        "results/summary/summary_samples_mqc.txt",
        "results/summary/summary_contigs_mqc.txt",
        "results/summary/summary_gene_counts_mqc.txt",
        get_plot_tree_output,
    output:
        "results/multiqc/multiqc_report.html"
    log:
        "logs/multiqc/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        # https://github.com/MultiQC/MultiQC/issues/2138
        sed -i -e 's/--in1/--i/g' results/fastp/*.json
        # multiqc
        multiqc \
            results/fastp \
            results/minimap/ \
            results/summary/ \
            --force \
            --config config/multiqc.yaml \
            --outdir results/multiqc &> {log}
        """
