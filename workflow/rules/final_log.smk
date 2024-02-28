# create final log when complete
rule final_log:
    input:
        get_plot_tree_output
    output:
        "results/snakemake.ok"
    conda:
        "../envs/conda_env.yaml"
    shell:
        """
        touch {output}
        """
