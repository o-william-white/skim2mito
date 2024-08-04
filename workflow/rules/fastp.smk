if fastp_dedup == "True":
    extra_params = "--dedup  --trim_poly_g"
else:
    extra_params = "--trim_poly_g"


rule fastp:
    input:
        sample=get_fastq,
    output:
        trimmed=["results/fastp/{sample}_R1.fq.gz", "results/fastp/{sample}_R2.fq.gz"],
        unpaired1="results/fastp/{sample}_u1.fq.gz",
        unpaired2="results/fastp/{sample}_u2.f.gz",
        failed="results/fastp/{sample}.failed.fq.gz",
        html="results/fastp/{sample}_fastp.html",
        json="results/fastp/{sample}_fastp.json",
    log:
        "logs/fastp/{sample}.log",
    params:
        adapters=expand(
            "--adapter_sequence {fwd} --adapter_sequence_r2 {rev}",
            fwd=forward_adapter,
            rev=reverse_adapter,
        ),
        extra=extra_params,
    threads: 1
    wrapper:
        "v3.3.6/bio/fastp"
