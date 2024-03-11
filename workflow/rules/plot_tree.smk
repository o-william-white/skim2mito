rule plot_tree:
    input:
        tree="results/iqtree/{dataset}.treefile.rooted.newick",
    params:
        height=plot_height,
        width=plot_width,
    output:
        png="results/plot_tree/{dataset}.png",
    log:
        "logs/plot_tree/{dataset}.log",
    conda:
        "../envs/r_env.yaml"
    shell:
        """
        # check if file empty
        if [ -s {input} ]; then
           # file not empty
            Rscript workflow/scripts/plot_tree.R {input.tree} {output.png} {params.height} {params.width} &> {log}
        else
           # file empty
           touch {output}
        fi
        """
