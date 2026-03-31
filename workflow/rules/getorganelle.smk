if go_reference == "go_fetch":

    rule getorganelle:
        input:
            expand(f"{results_dir}/go_fetch/{{taxids}}/gene.fasta", taxids=list(set(sample_data["taxid"]))),
            expand(f"{results_dir}/go_fetch/{{taxids}}/seed.fasta", taxids=list(set(sample_data["taxid"]))),
            fwd=f"{results_dir}/fastp/{{sample}}_R1.fastq.gz",
            rev=f"{results_dir}/fastp/{{sample}}_R2.fastq.gz",
        params:
            taxid=get_taxid,
            results_dir=results_dir,
            kmer_vals=config["kmer_vals"],
        output:
            directory(f"{results_dir}/getorganelle/{{sample}}/"),
        log:
            "logs/getorganelle/{sample}.log",
        conda:
            "../envs/getorganelle.yaml"
        shell:
            """
            get_organelle_from_reads.py \
                -1 {input.fwd} \
                -2 {input.rev} \
                -o {params.results_dir}/getorganelle/{wildcards.sample} \
                -F animal_mt \
                -s {params.results_dir}/go_fetch/{params.taxid}/seed.fasta \
                --genes {params.results_dir}/go_fetch/{params.taxid}/gene.fasta \
                --reduce-reads-for-coverage inf \
                --max-reads inf \
                -k {params.kmer_vals} \
                -R 20 \
                --overwrite &> {log}
            """

else:
    if go_reference == "custom":

        rule getorganelle:
            input:
                fwd=f"{results_dir}/fastp/{{sample}}_R1.fastq.gz",
                rev=f"{results_dir}/fastp/{{sample}}_R2.fastq.gz",
            params:
                seed=get_seed,
                gene=get_gene,
                results_dir=results_dir,
                kmer_vals=kmer_vals,
            output:
                directory(f"{results_dir}/getorganelle/{{sample}}/"),
            log:
                "logs/getorganelle/{sample}.log",
            conda:
                "../envs/getorganelle.yaml"
            shell:
                """
                get_organelle_from_reads.py \
                    -1 {input.fwd} \
                    -2 {input.rev} \
                    -o {params.results_dir}/getorganelle/{wildcards.sample} \
                    -F animal_mt \
                    -s {params.seed} \
                    --genes {params.gene} \
                    --reduce-reads-for-coverage inf \
                    --max-reads inf \
                    -k {params.kmer_vals} \
                    -R 20 \
                    --overwrite &> {log}
                """
