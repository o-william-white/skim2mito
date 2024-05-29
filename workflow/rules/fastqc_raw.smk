rule fastqc_raw_fwd:
    input:
        get_forward,
    output:
        html="results/fastqc_raw/{sample}_R1.html",
        zip="results/fastqc_raw/{sample}_R1_fastqc.zip",
    log:
        "logs/fastqc_raw/{sample}_R1.log",
    threads: 1
    wrapper:
        "v3.3.6/bio/fastqc"


rule fastqc_raw_rev:
    input:
        get_reverse,
    output:
        html="results/fastqc_raw/{sample}_R2.html",
        zip="results/fastqc_raw/{sample}_R2_fastqc.zip",
    log:
        "logs/fastqc_raw/{sample}_R2.log",
    threads: 1
    wrapper:
        "v3.3.6/bio/fastqc"
