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
        fas=f"{results_dir}/assembled_sequence/{{sample}}.fasta",
    output:
        f"{results_dir}/blastn/{{sample}}.txt",
    log:
        "logs/blastn/{sample}.log",
    conda:
        "../envs/blastn.yaml"
    shell:
        """
        blastn \
            -query {input.fas} \
            -db resources/blastdb/refseq_mitochondrion/refseq_mitochondrion \
            -out {output} \
            -outfmt '6 qseqid staxids bitscore std' \
            -max_target_seqs 10 \
            -max_hsps 1 \
            -evalue 1e-25 &> {log}
        """
