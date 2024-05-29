rule annotations:
    input:
        expand("resources/mitos_db/{refseq}", refseq=mitos_refseq),
        fas="results/assembled_sequence/{sample}.fasta",
    params:
        refseq=mitos_refseq,
        code=mitos_code,
    output:
        directory("results/annotations/{sample}/"),
        ok="results/annotations/{sample}/{sample}.ok",
    log:
        "logs/annotations/{sample}.log",
    conda:
        "../envs/annotations.yaml"
    shell:
        """
        if [ $(grep circular -c {input.fasta}) -eq 1 ] ; then
            echo Treating mitochondrial seqeunce as circular &> {log}
            runmitos.py \
                --input {input.fas} \
                --code {params.code} \
                --outdir results/annotations/{wildcards.sample}/ \
                --refseqver resources/mitos_db/{params.refseq} \
                --refdir . &>> {log}
        else
            echo Treating mitochndrial seqeunce as linear &> {log}
            runmitos.py \
                --input {input.fas} \
                --code {params.code} \
                --outdir results/annotations/{wildcards.sample}/ \
                --refseqver resources/mitos_db/{params.refseq} \
                --refdir . \
                --linear &>> {log}
        fi
        touch {output.ok}
        """
