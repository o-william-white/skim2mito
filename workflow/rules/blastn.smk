rule blastn:
    input:
        multiext(
            "resources/blastdb/refseq_mitochondrion/refseq_mitochondrion",
            ".ndb",
            ".nhr",
            ".nin",
            ".njs",
            ".nog",
            ".nos",
            ".not",
            ".nsq",
            ".ntf",
            ".nto",
        ),
        "results/assembled_sequence/{sample}.ok",
    output:
        ok="results/blastn/{sample}.ok",
    log:
        "logs/blastn/{sample}.log",
    conda:
        "../envs/blastn.yaml"
    shell:
        """
        FAS=$(echo results/assembled_sequence/{wildcards.sample}.fasta)
        OUT=$(echo results/blastn/{wildcards.sample}.txt)
        if [ -e $FAS ]; then
            echo Running blastn for {wildcards.sample} > {log}
            blastn \
               -query $FAS \
               -db resources/blastdb/refseq_mitochondrion/refseq_mitochondrion \
               -out $OUT \
               -outfmt '6 qseqid staxids bitscore std' \
               -max_target_seqs 10 \
               -max_hsps 1 \
               -evalue 1e-25 &> {log}
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """
