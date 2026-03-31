rule fastqc_raw_rev:
    input:
        get_reverse,
    output:
        html=f"{results_dir}/fastqc_raw/{{sample}}_R2.html",
        zip=f"{results_dir}/fastqc_raw/{{sample}}_R2_fastqc.zip",
    log:
        "logs/fastqc_raw/{sample}_R2.log",
    wrapper:
        "v7.5.0/bio/fastqc"
