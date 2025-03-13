rule fastqc_raw_rev:
    input:
        get_reverse,
    output:
        html="results/fastqc_raw/{sample}_R2.html",
        zip="results/fastqc_raw/{sample}_R2_fastqc.zip",
    log:
        "logs/fastqc_raw/{sample}_R2.log",
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v4.3.0/bio/fastqc"

