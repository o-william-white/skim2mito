rule fastqc_qc_rev:
    input:
        f"{results_dir}/fastp/{{sample}}_R2.fastq.gz",
    output:
        html=f"{results_dir}/fastqc_qc/{{sample}}_R2.html",
        zip=f"{results_dir}/fastqc_qc/{{sample}}_R2_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/fastqc_qc/{sample}_R2.log",
    wrapper:
        "v7.5.0/bio/fastqc"
