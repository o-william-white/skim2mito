if go_reference == "go_fetch":

    rule getorganelle:
        input:
            expand("results/go_fetch/{taxids}/gene.fasta", taxids=list(set(sample_data["taxid"]))),
            expand("results/go_fetch/{taxids}/seed.fasta", taxids=list(set(sample_data["taxid"]))),
            fwd="results/fastp/{sample}_R1.fastq.gz",
            rev="results/fastp/{sample}_R2.fastq.gz",
        params:
            taxid=get_taxid,
        output:
            directory("results/getorganelle/{sample}/"),
        log:
            "logs/getorganelle/{sample}.log",
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
            """

else:
    if go_reference == "custom":

        rule getorganelle:
            input:
                fwd="results/fastp/{sample}_R1.fastq.gz",
                rev="results/fastp/{sample}_R2.fastq.gz",
            params:
                seed=get_seed,
                gene=get_gene,
            output:
                directory("results/getorganelle/{sample}/"),
            log:
                "logs/getorganelle/{sample}.log",
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
                """
