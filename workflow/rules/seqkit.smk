rule seqkit:
    input:
        "results/assembled_sequence/{sample}.fasta",
    output:
        "results/seqkit/{wildcards.sample}.txt",
    log:
        "logs/seqkit/{sample}.log",
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit stats -b {input} > {output} 2> {log}
        """
