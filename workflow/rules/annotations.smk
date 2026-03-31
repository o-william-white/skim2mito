rule annotations:
    input:
        expand("resources/mitos_db/{refseq}", refseq=mitos_refseq),
        fasta=f"{results_dir}/assembled_sequence/{{sample}}.fasta",
    params:
        refseq=mitos_refseq,
        code=mitos_code,
        results_dir=results_dir,
    output:
        directory(f"{results_dir}/annotations/{{sample}}/"),
    log:
        "logs/annotations/{sample}.log",
    conda:
        "../envs/annotations.yaml"
    shell:
        """
        mkdir -p {output}
        # if the "grep -c" count is 0, it returns an exit code of 1 which causes snakemake to error silently
        # || true makes grep always return an exit code of 0  
        circular_count=$(grep -c circular {input.fasta} || true)
        sequence_count=$(grep -c "^>" {input.fasta} || true)
        
        if [ "$circular_count" -eq 1 ] && [ "$sequence_count" -eq 1 ] ; then
            echo Treating mitochondrial sequence as circular &> {log}
            runmitos.py \
                --input {input.fasta} \
                --code {params.code} \
                --outdir {params.results_dir}/annotations/{wildcards.sample}/ \
                --refseqver resources/mitos_db/{params.refseq} \
                --refdir . \
                --noplots &>> {log}
        else
            echo Treating mitochondrial sequence as linear &> {log}
            runmitos.py \
                --input {input.fasta} \
                --code {params.code} \
                --outdir {params.results_dir}/annotations/{wildcards.sample}/ \
                --refseqver resources/mitos_db/{params.refseq} \
                --refdir . \
                --linear \
                --noplots &>> {log}
        fi
        """
