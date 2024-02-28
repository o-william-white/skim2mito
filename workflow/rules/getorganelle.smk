if go_reference == "go_fetch":
    rule getorganelle:
        input:
            fwd = "results/fastp/{sample}_R1.fastq",
            rev = "results/fastp/{sample}_R2.fastq"
        params:
            taxid = get_taxid,
        output:
            ok = "results/getorganelle/{sample}/getorganelle.ok"
        log:
            "logs/getorganelle/{sample}.log"
        conda:
            "../envs/getorganelle.yaml"
        shell:
            """
            get_organelle_from_reads.py \
                -1 {input.fwd} \
                -2 {input.rev} \
                -o results/getorganelle/{wildcards.sample} \
                -F animal_mt \
                -s results/go_fetch/{params.taxid}/seed.fasta \
                --genes results/go_fetch/{params.taxid}/gene.fasta \
                --reduce-reads-for-coverage inf \
                --max-reads inf \
                -R 20 \
                --overwrite &> {log}
            touch {output.ok}
            """
else: 
    if go_reference == "custom":
        rule getorganelle:
            input:
                fwd = "results/fastp/{sample}_R1.fastq",
                rev = "results/fastp/{sample}_R2.fastq"
            params:
                seed = get_seed,
                gene = get_gene,
            output:
                ok = "results/getorganelle/{sample}/getorganelle.ok"
            log:
                "logs/getorganelle/{sample}.log"
            conda:
                "../envs/getorganelle.yaml"
            shell:
                """
                get_organelle_from_reads.py \
                    -1 {input.fwd} \
                    -2 {input.rev} \
                    -o results/getorganelle/{wildcards.sample} \
                    -F animal_mt \
                    -s {params.seed} \
                    --genes {params.gene} \
                    --reduce-reads-for-coverage inf \
                    --max-reads inf \
                    -R 20 \
                    --overwrite &> {log}
                touch {output.ok}
                """

