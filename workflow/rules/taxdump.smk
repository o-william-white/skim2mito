rule taxdump:
    # if we give the blastdb database as input, it stops wget trying to download the blastdb at the same time as taxdump which causes an error
    input:
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
    output:
        directory("resources/taxdump"),
        "resources/taxdump/citations.dmp",
        "resources/taxdump/delnodes.dmp",
        "resources/taxdump/division.dmp",
        "resources/taxdump/excludedfromtype.dmp",
        "resources/taxdump/fullnamelineage.dmp",
        "resources/taxdump/gencode.dmp",
        "resources/taxdump/host.dmp",
        "resources/taxdump/images.dmp",
        "resources/taxdump/merged.dmp",
        "resources/taxdump/names.dmp",
        "resources/taxdump/nodes.dmp",
        "resources/taxdump/rankedlineage.dmp",
        "resources/taxdump/taxidlineage.dmp",
        "resources/taxdump/typematerial.dmp",
        "resources/taxdump/typeoftype.dmp"
    log:
        "logs/taxdump/taxdump.log"
    conda:
        "../envs/conda_env.yaml"
    shell:
        """
        wget --wait 10 --random-wait -P resources/taxdump/ https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz &> {log}
        tar xvzf resources/taxdump/new_taxdump.tar.gz --directory resources/taxdump/ &>> {log}
        rm resources/taxdump/new_taxdump.tar.gz &>> {log}
        """
