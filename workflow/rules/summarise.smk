rule summarise:
    input:
        get_seqkit_output,
        get_blobtools_output,
        get_assess_assembly_output,
    params:
        config = config["samples"]
    output:
        table_sample=f"{results_dir}/summary/summary_samples_mqc.txt",
        table_contig=f"{results_dir}/summary/summary_contigs_mqc.txt",
    log:
        "logs/summarise/summarise.log",
    conda:
        "../envs/r_env.yaml"
    shell:
        """
        Rscript workflow/scripts/summarise.R {params.config} &> {log}
        """
