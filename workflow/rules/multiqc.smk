rule multiqc:
    input:
        get_seqkit_output,
        get_blobtools_output,
        get_annotated_samples,
        get_assess_assembly_output,
        "results/summary/summary_samples.txt",
        "results/summary/summary_contigs.txt",
        "results/summary/summary_gene_counts.txt",
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
        # change file endings 
        # change file endings to match multiqc
        cp results/summary/summary_contigs.txt results/summary/summary_contigs_mqc.txt
        cp results/summary/summary_samples.txt results/summary/summary_samples_mqc.txt
        cp results/summary/summary_gene_counts.txt results/summary/summary_gene_counts_mqc.txt
        # multiqc
        multiqc \
            results/fastp \
            results/minimap/ \
            results/summary/ \
            --force \
            --config config/multiqc.yaml \
            --outdir results/multiqc &> {log}
        """
