import os

# set configfile
configfile: "config/config.yaml"

# configfile parameters
input_dir = config["input_dir"]
output_dir = config["output_dir"]
realign = config["realign"]
missing_threshold = config["missing_threshold"]
alignment_trim = config["alignment_trim"]
outgroup = config["outgroup"]
plot_height = config["plot_height"]
plot_width = config["plot_width"]

# paramter checks
if realign != True and realign != False: 
    sys.exit("Error: realign must be True or False")
if realign == True:
    if not isinstance(missing_threshold, float) or missing_threshold < 0.0 or missing_threshold > 1.0:
        sys.exit("Error: if realign == True, missing_threshold must be a float between 0.0 and 1.0")
    if alignment_trim not in ["gblocks", "clipkit"]:
        sys.exit("Error: if realign == True, alignment_trim must be 'gblocks' or 'clipkit'")

# get a list of all possible input genes
GENES = []
for file in os.listdir(input_dir):
    if file.endswith(".fasta"):
       g = file.replace(".fasta", "")
       GENES.append(g)

# one rule to rule them all :)
rule all:
    input:
        expand("{out}/iqtree_plots/{gene}.png", out=output_dir, gene=GENES),
        f"{output_dir}/iqtree_partitioned_plot/iqtree_partitioned_plot.png",
        f"{output_dir}/astral_plot/astral_plot.png"

# to do
# test it works with and without realign 
# add param checks to smk

###

if realign == True: 

    rule mafft:
        input:
            input_dir+"/{genes}.fasta"
        output:
            "{output_dir}/mafft/{genes}.fasta"
        log:
            "{output_dir}/logs/mafft/{genes}.log"
        conda:
            "envs/mafft.yaml"
        shell:
            """
            mafft \
                --maxiterate 1000 \
                --globalpair \
                --adjustdirectionaccurately \
                {input} 1> {output} 2> {log}
            """

    rule filter_alignments:
        input:
            "{output_dir}/mafft/{genes}.fasta"
        params:
            threshold = missing_threshold
        output:
            "{output_dir}/mafft_filtered/{genes}.fasta"
        log:
            "{output_dir}/logs/mafft_filtered/{genes}.log"
        shell:
            """
            python scripts/alignments_filter.py --input {input} --output {output} --threshold {params.threshold} > {log}
            """

    rule alignment_trim:
        input:
            "{output_dir}/mafft_filtered/{genes}.fasta"
        params:
            tmp = "{output_dir}/alignment_trim/{genes}_tmp.fasta"
        output:
            out = "{output_dir}/alignment_trim/{genes}.fasta"
        log:
            "{output_dir}/logs/alignment_trim/{genes}.log"
        conda:
            "envs/alignment_trim.yaml"
        shell:
            """
            if [ $(grep -c "^>" {input}) -lt "5" ]; then
                cp {input} {output.out}
            else
                # if [ $(grep -c "^>" {input[0]}) -lt "0" ]; then
                if [[ {alignment_trim} == "gblocks" ]]; then
                    # gblocks add reuslts to same dir as input
                    cp {input} {params.tmp}
                    # gblocks always gives error code of 1. Ignore.
                    Gblocks {params.tmp} -t=d -b4=5 -b5=h &> {log} || true
                    # sed to remove gaps
                    sed 's/ //g' {params.tmp}-gb > {output.out}
                    # rm tmp
                    rm {params.tmp}
                    rm {params.tmp}-gb
                else
                    if [[ {alignment_trim} == "clipkit" ]]; then
                        clipkit {input} -o {output.out} &> {log}
                    fi
                fi
            fi
            """

    rule iqtree:
        input:
            fasta = "{output_dir}/alignment_trim/{genes}.fasta"
        output:
            tree = "{output_dir}/iqtree/{genes}.treefile"
        log:
            "{output_dir}/logs/iqtree/{genes}.log"
        conda:
            "envs/iqtree.yaml"
        shell:
            """
            iqtree -s {input.fasta} -B 1000 --prefix {output_dir}/iqtree/{wildcards.genes} &> {log}
            """

else:

    rule iqtree:
        input:
            fasta = input_dir+"/{genes}.fasta"
        output:
            tree = "{output_dir}/iqtree/{genes}.treefile"
        log:
            "{output_dir}/logs/iqtree/{genes}.log"
        conda:
            "envs/iqtree.yaml"
        shell:
            """
            iqtree -s {input.fasta} -B 1000 --prefix {output_dir}/iqtree/{wildcards.genes} &> {log}
            """

rule root_iqtree:
    input:
        tree = "{output_dir}/iqtree/{genes}.treefile"
    output:
        tree = "{output_dir}/iqtree/{genes}.treefile.rooted.newick"
    log:
        "{output_dir}/logs/root_iqtree/{genes}.txt"
    conda:
        "envs/ete3.yaml"
    shell:
        """
        if [ {outgroup} == "NA" ]; then
            echo "Outgroup not specified. Leaving as unrooted" > {log}
            cp {input.tree} {output.tree}
        else
            python scripts/root_newick.py \
                --input {input.tree} \
                --output {output.tree} \
                --outgroup {outgroup} &> {log}
        fi
        """

rule plot_iqtree:
    input:
        tree = "{output_dir}/iqtree/{genes}.treefile.rooted.newick"
    output:
        png = "{output_dir}/iqtree_plots/{genes}.png"
    log:
        "{output_dir}/logs/iqtree_plots/{genes}.txt"
    conda:
        "envs/r_env.yaml"
    shell:
        """
        Rscript scripts/plot_tree.R {input.tree} {output.png} {plot_height} {plot_width} &> {log}
        """

rule concatenate_alignments: 
    output:
        "{output_dir}/concatenate_alignments/output.fasta",
        "{output_dir}/concatenate_alignments/output.txt"
    log: 
        "{output_dir}/logs/concatentate_alignments/log.txt"
    shell:
        """
        python scripts/concatenate_alignments.py \
            --input {input_dir} \
            --type DNA \
            --output {output_dir}/concatenate_alignments/ \
            --prefix output \
            --overwrite
        """

rule iqtree_partitioned:
    input:
        fasta = "{output_dir}/concatenate_alignments/output.fasta",
        partitions = "{output_dir}/concatenate_alignments/output.txt" 
    output:
        "{output_dir}/iqtree_partitioned/output.treefile"
    log:
        "{output_dir}/logs/iqtree_partitioned/log.txt"
    conda:
        "envs/iqtree.yaml"
    shell:
        """
        iqtree \
            -s {input.fasta} \
            -p {input.partitions} \
            -B 1000 \
            --prefix {output_dir}/iqtree_partitioned/output &> {log}
        """

rule root_iqtree_partitioned:
    input: 
        tree = "{output_dir}/iqtree_partitioned/output.treefile"
    output:
        tree = "{output_dir}/iqtree_partitioned/output.treefile.rooted.newick"
    log:
        "{output_dir}/logs/root_iqtree_partitioned/log.txt"
    conda:
        "envs/ete3.yaml"
    shell:
        """
        if [ {outgroup} == "NA" ]; then
            echo "Outgroup not specified. Leaving as unrooted" > {log}
            cp {input.tree} {output.tree}
        else
            python scripts/root_newick.py \
                 --input {input.tree} \
                 --output {output.tree} \
                 --outgroup {outgroup} &> {log}
        fi
        """

rule plot_iqtree_partitioned:
    input:
        tree = "{output_dir}/iqtree_partitioned/output.treefile.rooted.newick"
    output:
        png = "{output_dir}/iqtree_partitioned_plot/iqtree_partitioned_plot.png"
    log:
        "{output_dir}/logs/iqtree_partitioned_plot/log.txt"
    conda:
        "envs/r_env.yaml"
    shell:
        """
        Rscript scripts/plot_tree.R {input.tree} {output.png} {plot_height} {plot_width} &> {log}
        """

rule astral_cat_input:
    input:
        tree = expand("{out}/iqtree/{gene}.treefile", out=output_dir, gene=GENES)
    output:
        input_trees = "{output_dir}/astral/input_trees.tree"
    log:
        "{output_dir}/logs/astral/log.txt"
    conda:
        "envs/astral.yaml"
    shell:
        """
        cat {input.tree} > {output.input_trees} 2> {log}
        """

rule astral:
    input:
        tree = "{output_dir}/astral/input_trees.tree"
    output:
        tree = "{output_dir}/astral/output.tree"
    log:
        "{output_dir}/logs/astral/log.txt"
    conda:
        "envs/astral.yaml"
    shell:
        """
        (java -jar $CONDA_PREFIX/share/astral-tree-5.7.8-0/astral.5.7.8.jar -i {input.tree} -o {output.tree}) &> {log}
        #astral --input {input.tree} --output {output.tree}
        """

rule root_astral:
    input:
        tree = "{output_dir}/astral/output.tree"
    output:
        tree = "{output_dir}/astral/output.tree.rooted.newick"
    log:
        "{output_dir}/logs/root_astral/log.txt"
    conda:
        "envs/ete3.yaml"
    shell:
        """
        if [ {outgroup} == "NA" ]; then
            echo "Outgroup not specified. Leaving as unrooted" > {log}
            cp {input.tree} {output.tree}
        else
            python scripts/root_newick.py \
                --input {input.tree} \
                --output {output.tree} \
                --outgroup {outgroup} &> {log}
        fi
        """

rule plot_astral:
    input:
        tree = "{output_dir}/astral/output.tree.rooted.newick"
    output:
        png = "{output_dir}/astral_plot/astral_plot.png"
    log:
        "{output_dir}/logs/astral_plot/log.txt"
    conda:
        "envs/r_env.yaml"
    shell:
        """
        Rscript scripts/plot_tree.R {input.tree} {output.png} {plot_height} {plot_width} &> {log}
        """

