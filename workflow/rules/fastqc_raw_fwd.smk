rule fastqc_raw_fwd:
    input:
        get_forward,
    output:
        html="results/fastqc_raw/{sample}_R1.html",
        zip="results/fastqc_raw/{sample}_R1_fastqc.zip",
    log:
        "logs/fastqc_raw/{sample}_R1.log",
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v4.3.0/bio/fastqc"
