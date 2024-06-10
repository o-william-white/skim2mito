rule summarise:
    input:
        get_seqkit_output,
        get_blobtools_output,
        get_assess_assembly_output,
    output:
        table_sample="results/summary/summary_samples.txt",
        table_contig="results/summary/summary_contigs.txt",
    log:
        "logs/summarise/summarise.log",
    conda:
        "../envs/r_env.yaml"
    shell:
        """
        Rscript workflow/scripts/summarise.R &> {log}
        """
