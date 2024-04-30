rule summarise:
    input:
        expand("results/seqkit/{sample}.ok", sample=sample_data["ID"].tolist()),
        expand(
            "results/blobtools/{sample}/{sample}_filter.ok", sample=sample_data["ID"].tolist()
        ),
        expand(
            "results/annotations/{sample}/{sample}.ok",
            sample=sample_data["ID"].tolist(),
        ),
        expand("results/assess_assembly/{sample}.ok", sample=sample_data["ID"].tolist()),
    output:
        table_sample="results/summary/summary_sample.txt",
        table_contig="results/summary/summary_contig.txt",
    log:
        "logs/summarise/summarise.log",
    conda:
        "../envs/r_env.yaml"
    shell:
        """
        # cat seqkit output for each sample
        echo -e "sample format type num_seqs sum_len min_len avg_len max_len" > results/summary/tmp_summary_sample.txt
        cat results/seqkit/*.txt | grep file -v >> results/summary/tmp_summary_sample.txt
        column -t results/summary/tmp_summary_sample.txt > {output.table_sample}
        rm results/summary/tmp_summary_sample.txt
        # join blobtools with mitos annotations for each contig
        Rscript workflow/scripts/summarise.R results/ mitos {output.table_contig} &> {log}
        """
