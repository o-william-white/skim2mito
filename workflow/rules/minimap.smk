rule minimap:
    input:
        fas=f"{results_dir}/assembled_sequence/{{sample}}.fasta",
        fwd=f"{results_dir}/fastp/{{sample}}_R1.fastq.gz",
        rev=f"{results_dir}/fastp/{{sample}}_R2.fastq.gz",
    output:
        bam = f"{results_dir}/minimap/{{sample}}.bam",
        stats = f"{results_dir}/minimap/{{sample}}_stats.txt"
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
