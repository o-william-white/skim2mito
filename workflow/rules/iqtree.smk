rule iqtree:
    input:
        fasta="results/alignment_trim/{dataset}.fasta",
    output:
        tree="results/iqtree/{dataset}.treefile",
        fasta_renamed="results/iqtree/{dataset}.fasta",
    log:
        "logs/iqtree/{dataset}.log",
    conda:
        "../envs/iqtree.yaml"
    shell:
        """
        # remove special characters from sample names
        sed -e 's/;/_/g' -e 's/+//g' \
            {input.fasta} > {output.fasta_renamed}

        # iqtree will not bootstrap if less than 5 samples in alignment
        if [ $(grep -c "^>" {input.fasta}) -lt "5" ] || [ $(grep -e "^>" -v {input.fasta} | sort | uniq | wc -l)  -lt 5 ] ; then
            touch {output.tree}
        else
            iqtree -s {output.fasta_renamed} -B 1000 --prefix results/iqtree/{wildcards.dataset} -redo &> {log}
        fi
        """
