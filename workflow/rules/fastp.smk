if fastp_dedup == "True":
    extra_params = "--dedup  --trim_poly_g"
else:
    extra_params = "--trim_poly_g"

rule fastp_pe:
    input:
        sample = get_fastq
    output:
        trimmed = ["results/fastp/{sample}_R1.fastq", "results/fastp/{sample}_R2.fastq"],
        unpaired1 = "results/fastp/{sample}_u1.fastq",
        unpaired2 = "results/fastp/{sample}_u2.fastq",
        failed = "results/fastp/{sample}.failed.fastq",
        html = "results/fastp/{sample}.html",
        json = "results/fastp/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        adapters = expand("--adapter_sequence {fwd} --adapter_sequence_r2 {rev}", fwd=forward_adapter, rev=reverse_adapter),
        extra = extra_params
    threads: 2
    wrapper:
        "v3.3.6/bio/fastp"
