rule assembled_sequence:
    input:
        "results/getorganelle/{sample}/getorganelle.ok"
    output:
        ok = "results/assembled_sequence/{sample}.ok"
    log:
        "logs/assembled_sequence/{sample}.log"
    conda:
        "../envs/conda_env.yaml"
    shell:
        """
        # find selected path(s) fasta
        FAS=$(find results/getorganelle/{wildcards.sample}/ -name *path_sequence.fasta)
        # z option: true if length if string is zero.
        if [[ -z $FAS ]]; then
            echo No assembly produced for {wildcards.sample} > {log}
        # more than one selected path
        elif [ "$(echo $FAS | tr ' ' '\\n' | wc -l)" -gt 1 ]; then
            FAS1=$(echo $FAS | tr ' ' '\\n' | head -n 1)
            echo More than one assembly produced for {wildcards.sample} > {log}
            echo Selecting the first assembly $FAS1 > {log}
            python workflow/scripts/rename_assembled.py \
                --input $FAS1 \
                --sample {wildcards.sample} \
                --output results/assembled_sequence
        elif [ "$(echo $FAS | tr ' ' '\\n' | wc -l)" -eq 1 ]; then
            echo One assembly produced for {wildcards.sample} > {log}
            python workflow/scripts/rename_assembled.py \
                --input $FAS \
                --sample {wildcards.sample} \
                --output results/assembled_sequence
        fi
        touch {output.ok}
        """
