rule minimap:
    input:
        fas="results/assembled_sequence/{sample}.fasta",
        fwd="results/fastp/{sample}_R1.fq.gz",
        rev="results/fastp/{sample}_R2.fq.gz",
    output:
        bam = "results/minimap/{sample}.bam",
        stats = "results/minimap/{sample}_stats.txt"
    log:
        "logs/minimap/{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax sr {input.fas} {input.fwd} {input.rev} 2> {log} | samtools view -b -F 4 | samtools sort -O BAM -o {output.bam} - 2>> {log}
        samtools index {output.bam} 2>> {log}
        samtools index -c {output.bam} 2>> {log}
        samtools stats {output.bam} > {output.stats}
        """
