if fastp_dedup == "True":
    extra_params = "--dedup  --trim_poly_g"
else:
    extra_params = "--trim_poly_g"


rule fastp:
    input:
        sample=get_fastq,
    output:
        trimmed=[f"{results_dir}/fastp/{{sample}}_R1.fastq.gz", f"{results_dir}/fastp/{{sample}}_R2.fastq.gz"],
        unpaired1=f"{results_dir}/fastp/{{sample}}_u1.fastq.gz",
        unpaired2=f"{results_dir}/fastp/{{sample}}_u2.fastq.gz",
        failed=f"{results_dir}/fastp/{{sample}}.failed.fastq.gz",
        html=f"{results_dir}/fastp/{{sample}}_fastp.html",
        json=f"{results_dir}/fastp/{{sample}}_fastp.json",
    log:
        "logs/fastp/{sample}.log",
    params:
        adapters=expand(
            "--adapter_sequence {fwd} --adapter_sequence_r2 {rev}",
            fwd=forward_adapter,
            rev=reverse_adapter,
        ),
        extra=extra_params,
    threads: 2
    wrapper:
        "v7.5.0/bio/fastp"
