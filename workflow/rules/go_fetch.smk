rule go_batch:
    params: 
        email = user_email
    output:
        "results/go_fetch/{taxids}/gene.fasta",
        "results/go_fetch/{taxids}/seed.fasta"
    log:
        "logs/go_fetch/{taxids}.log"
    conda:
        "../envs/go_fetch.yaml"
    shell:
        """
        python3 workflow/scripts/go_fetch.py \
            --taxonomy {wildcards.taxids} \
            --target mitochondrion \
            --db genbank \
            --min 5  \
            --max 10 \
            --output results/go_fetch/{wildcards.taxids} \
            --getorganelle \
            --email o.william.white@gmail.com \
            --overwrite &> {log}
        """
