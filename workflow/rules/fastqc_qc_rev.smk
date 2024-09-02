rule fastqc_qc_rev:
    input:
        "results/fastp/{sample}_R2.fastq.gz",
    output:
        html="results/fastqc_qc/{sample}_R2.html",
        zip="results/fastqc_qc/{sample}_R2_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/fastqc_qc/{sample}_R2.log",
    threads: 1
    resources:
        mem_mb = 1024
    wrapper:
        "v4.3.0/bio/fastqc"
