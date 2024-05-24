checkpoint assembled_sequence:
    input:
        expand("results/getorganelle/{sample}/", sample=sample_data["ID"]),
    output:
        directory("results/assembled_sequence/"),
    log:
        "logs/assembled_sequence.log",
    conda:
        "../envs/conda_env.yaml"
    script:
        "../scripts/rename_assembled.py"
