rule fastqc_qc_fwd:
    input:
        f"{results_dir}/fastp/{{sample}}_R1.fastq.gz",
    output:
        html=f"{results_dir}/fastqc_qc/{{sample}}_R1.html",
        zip=f"{results_dir}/fastqc_qc/{{sample}}_R1_fastqc.zip",
    params:
        extra="--quiet",
    log:
        "logs/fastqc_qc/{sample}_R1.log",
    wrapper:
        "v7.5.0/bio/fastqc"
