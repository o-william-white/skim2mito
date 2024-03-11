rule mitos_db:
    # if we give the taxdump database as input, it stops wget trying to download the taxdump mitos_db at the same time which causes an error
    input:
        "resources/taxdump/",
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
        "resources/taxdump/typeoftype.dmp",
    params:
        refseq=mitos_refseq,
    output:
        directory("resources/mitos_db/{mitos_refseq}"),
    log:
        "logs/mitos_db/mitos_db_{mitos_refseq}.log",
    conda:
        "../envs/conda_env.yaml"
    shell:
        """
        wget --wait 10 --random-wait -P resources/mitos_db https://zenodo.org/record/4284483/files/{params.refseq}.tar.bz2  &> {log}
        tar xf resources/mitos_db/{params.refseq}.tar.bz2 --directory resources/mitos_db &>> {log}
        rm resources/mitos_db/{params.refseq}.tar.bz2 >> {log}
        """
