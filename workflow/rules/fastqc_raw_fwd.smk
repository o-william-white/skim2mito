rule fastqc_raw_fwd:
    input:
        get_forward,
    output:
        html=f"{results_dir}/fastqc_raw/{{sample}}_R1.html",
        zip=f"{results_dir}/fastqc_raw/{{sample}}_R1_fastqc.zip",
    log:
        "logs/fastqc_raw/{sample}_R1.log",
    wrapper:
        "v7.5.0/bio/fastqc"
