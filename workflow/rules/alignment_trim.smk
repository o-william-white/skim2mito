rule alignment_trim:
    input:
        fasta = "results/mafft_filtered/{dataset}.fasta"
    params:
        trim = alignment_trim,
        tmp = "results/alignment_trim/{dataset}_tmp.fasta"
    output:
        fasta = "results/alignment_trim/{dataset}.fasta"
    log:
        "logs/alignment_trim/{dataset}.log"
    conda:
        "../envs/alignment_trim.yaml"
    shell:
        """
        if [ $(grep -c "^>" {input.fasta}) -lt "5" ]; then
            cp {input.fasta} {output.fasta}
        else
            if [[ {params.trim} == "gblocks" ]]; then
                # gblocks add reuslts to same dir as input
                cp {input.fasta} {params.tmp}
                # gblocks always gives error code of 1. Ignore.
                Gblocks {params.tmp} -t=d -b4=5 -b5=h &> {log} || true
                # sed to remove gaps
                sed 's/ //g' {params.tmp}-gb > {output.fasta}
                # rm tmp
                rm {params.tmp}
                rm {params.tmp}-gb
            else
                if [[ {params.trim} == "clipkit" ]]; then
                    clipkit {input.fasta} -o {output.fasta} &> {log}
                fi
            fi
        fi
        """
