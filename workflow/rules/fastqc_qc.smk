rule fastqc_qc_fwd:
    input:
        "results/fastp/{sample}_R1.fastq",
    output:
        html="results/fastqc_qc/{sample}_R1.html",
        zip="results/fastqc_qc/{sample}_R2_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/fastqc_qc/{sample}_R2.log",
    threads: 1
    resources:
        mem_mb=1024,
    wrapper:
        "v3.3.6/bio/fastqc"


rule fastqc_qc_rev:
    input:
        "results/fastp/{sample}_R2.fastq",
    output:
        html="results/fastqc_qc/{sample}_R2.html",
        zip="results/fastqc_qc/{sample}_R2_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/fastqc_qc/{sample}.log",
    threads: 1
    resources:
        mem_mb=1024,
    wrapper:
        "v3.3.6/bio/fastqc"
