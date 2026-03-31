rule go_fetch:
    params:
        email=config["user_email"],
        api=config["user_api"],
        gf_min=config["gf_min"],
        gf_max=config["gf_max"],
        results_dir=results_dir,
    output:
        f"{results_dir}/go_fetch/{{taxids}}/gene.fasta",
        f"{results_dir}/go_fetch/{{taxids}}/seed.fasta",
    log:
        "logs/go_fetch/{taxids}.log",
    conda:
        "../envs/go_fetch.yaml"
    threads:
        workflow.cores * 0.2
    shell:
        """
        python3 workflow/scripts/go_fetch.py \
            --taxonomy {wildcards.taxids} \
            --target mitochondrion \
            --db genbank \
            --min {params.gf_min}  \
            --max {params.gf_max} \
            --output {params.results_dir}/go_fetch/{wildcards.taxids} \
            --getorganelle \
            --email {params.email} \
            --api {params.api} \
            --overwrite &> {log}
        """
