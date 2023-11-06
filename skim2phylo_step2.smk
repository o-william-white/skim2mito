import os

# input params
input_dir = "results_genes_example"
output_dir = "results_genes_example_output"
outgroup = "Eurema_blanda"

# get a list of all possible input genes
GENES = []
for file in os.listdir(input_dir):
    if file.endswith(".fasta"):
       g = file.replace(".fasta", "")
       GENES.append(g)

print(GENES)

# one rule to rule them all :)
rule all:
    input:
        expand("{out}/iqtree/{gene}.treefile", out=output_dir, gene=GENES),
        f"{output_dir}/astral/output.tree", 
        f"{output_dir}/concatenate_alignments/output.fasta",
        f"{output_dir}/concatenate_alignments/output.txt",
        f"{output_dir}/iqtree_partitioned/output.treefile",
        f"{output_dir}/iqtree_partitioned/output.treefile.rooted.newick",
        expand("{out}/iqtree/{gene}.treefile.rooted.newick", out=output_dir, gene=GENES)

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
        if [ $(grep {outgroup} -c {input.tree}) == 0 ]
        then  
            echo "Outgroup not present in tree. Leaving as unrooted" > {log}
            cp {input.tree} {output.tree}
        else
            python scripts/root_newick.py \
                --input {input.tree} \
                --output {output.tree} \
                --outgroup {outgroup} &> {log}
        fi
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
        python scripts/root_newick.py \
            --input {input.tree} \
            --output {output.tree} \
            --outgroup {outgroup} &> {log}
        """


