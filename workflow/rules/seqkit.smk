rule seqkit:
    input:
        f"{results_dir}/assembled_sequence/{{sample}}.fasta",
    output:
        f"{results_dir}/seqkit/{{sample}}.txt",
    log:
        "logs/seqkit/{sample}.log",
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit stats -b {input} > {output} 2> {log}
        """
