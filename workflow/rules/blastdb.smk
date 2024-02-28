rule blastdb:
    output:
        multiext("resources/blastdb/refseq_mitochondrion/refseq_mitochondrion",
            ".ndb",
            ".nhr",
            ".nin",
            ".njs",
            ".nog",
            ".nos",
            ".not",
            ".nsq",
            ".ntf",
            ".nto")
    log:
        "logs/blastdb/blastdb.log"
    conda:
        "../envs/conda_env.yaml"
    shell:
        """
        wget --wait 10 --random-wait -P resources/blastdb/ https://zenodo.org/records/8424777/files/refseq_mitochondrion.tar.gz &> {log}
        tar xvzf resources/blastdb/refseq_mitochondrion.tar.gz --directory resources/blastdb/ &>> {log}
        rm resources/blastdb/refseq_mitochondrion.tar.gz &>> {log}
        """
