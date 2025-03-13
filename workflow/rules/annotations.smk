rule annotations:
    input:
        expand("resources/mitos_db/{refseq}", refseq=mitos_refseq),
        fasta="results/assembled_sequence/{sample}.fasta",
    params:
        refseq=mitos_refseq,
        code=mitos_code,
    output:
        directory("results/annotations/{sample}/"),
    log:
        "logs/annotations/{sample}.log",
    conda:
        "../envs/annotations.yaml"
    shell:
        """
        mkdir -p {output}
        if [ $(grep circular -c {input.fasta}) -eq 1 ] ; then
            echo Treating mitochondrial sequence as circular &> {log}
            runmitos.py \
                --input {input.fasta} \
                --code {params.code} \
                --outdir results/annotations/{wildcards.sample}/ \
                --refseqver resources/mitos_db/{params.refseq} \
                --refdir . &>> {log}
        else
            echo Treating mitochndrial sequence as linear &> {log}
            runmitos.py \
                --input {input.fasta} \
                --code {params.code} \
                --outdir results/annotations/{wildcards.sample}/ \
                --refseqver resources/mitos_db/{params.refseq} \
                --refdir . \
                --linear &>> {log}
        fi
        """
