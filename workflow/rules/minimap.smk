rule minimap:
    input:
        "results/assembled_sequence/{sample}.ok",
        fwd="results/fastp/{sample}_R1.fastq",
        rev="results/fastp/{sample}_R2.fastq",
    output:
        ok="results/minimap/{sample}.ok",
    log:
        "logs/minimap/{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        FAS=$(echo results/assembled_sequence/{wildcards.sample}.fasta)
        OUT=$(echo results/minimap/{wildcards.sample}.bam)
        STA=$(echo results/minimap/{wildcards.sample}_stats.txt)
        if [ -e $FAS ]; then
            echo Running minimap for {wildcards.sample} > {log}
            minimap2 -ax sr $FAS {input.fwd} {input.rev} 2> {log} | samtools view -b -F 4 | samtools sort -O BAM -o $OUT - 2>> {log}
            samtools index $OUT 2>> {log}
            samtools index -c $OUT 2>> {log}
            samtools stats $OUT > $STA
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """
