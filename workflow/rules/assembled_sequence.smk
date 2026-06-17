checkpoint assembled_sequence:
    input:
        expand(f"{results_dir}/getorganelle/{{sample}}/", sample=sample_data["ID"]),
    output:
        directory(f"{results_dir}/assembled_sequence/"),
    log:
        "logs/assembled_sequence.log",
    script:
        "../scripts/rename_assembled.py"
