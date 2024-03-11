rule seqkit:
    input:
        "results/assembled_sequence/{sample}.ok",
    output:
        ok="results/seqkit/{sample}.ok",
    log:
        "logs/seqkit/{sample}.log",
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        FAS=$(echo results/assembled_sequence/{wildcards.sample}.fasta)
        OUT=$(echo results/seqkit/{wildcards.sample}.txt)
        if [ -e $FAS ]; then
            echo Running seqkit for {wildcards.sample} > {log}
            seqkit stats -b $FAS > $OUT
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """
