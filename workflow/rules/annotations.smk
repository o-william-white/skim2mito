rule annotations:
    input:
        expand("resources/mitos_db/{refseq}", refseq=mitos_refseq),
        "results/assembled_sequence/{sample}.ok",
    params:
        refseq=mitos_refseq,
        code=mitos_code,
    output:
        ok="results/annotations/{sample}/{sample}.ok",
    log:
        "logs/annotations/{sample}.log",
    conda:
        "../envs/annotations.yaml"
    shell:
        """
        FAS=$(echo results/assembled_sequence/{wildcards.sample}.fasta)
        if [ -e $FAS ]; then
            if [ $(grep circular -c $FAS) -eq 1 ] ; then
                echo Treating mitochondrial seqeunce as circular &> {log}
                runmitos.py \
                    --input $FAS \
                    --code {params.code} \
                    --outdir results/annotations/{wildcards.sample}/ \
                    --refseqver resources/mitos_db/{params.refseq} \
                    --refdir . &>> {log}
            else
                echo Treating mitochndrial seqeunce as linear &> {log}
                runmitos.py \
                    --input $FAS \
                    --code {params.code} \
                    --outdir results/annotations/{wildcards.sample}/ \
                    --refseqver resources/mitos_db/{params.refseq} \
                    --refdir . \
                    --linear &>> {log}
            fi
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """
