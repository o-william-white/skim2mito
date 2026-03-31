rule assess_assembly:
    input:
        get_sucessfully_annotated_samples,
        fasta=f"{results_dir}/assembled_sequence/{{sample}}.fasta",
        bam=f"{results_dir}/minimap/{{sample}}.bam",
    params:
        results_dir=results_dir,
    output:
        directory(f"{results_dir}/assess_assembly/{{sample}}"),
    log:
        "logs/assess_assembly/{sample}.log",
    conda:
        "../envs/assess_assembly.yaml"
    shell:
        """
        if [ $(grep -e "^>" -c {input.fasta}) -eq 1 ] ; then
            echo Single sequence found in fasta > {log}
            python workflow/scripts/assess_assembly.py \
                --fasta {input.fasta} \
                --bam {input.bam} \
                --bed {params.results_dir}/annotations/{wildcards.sample}/result.bed \
                --sample {wildcards.sample} \
                --output {params.results_dir}/assess_assembly/{wildcards.sample} &> {log}
        else
            echo More than one sequence found in fasta > {log}
            # mitos creates subdirectories for each contig
            # find bed files and cat
            mkdir -p {params.results_dir}/assess_assembly/{wildcards.sample}
            find {params.results_dir}/annotations/{wildcards.sample}/ -type f -name result.bed | while read line; do 
                cat $line
            done > {params.results_dir}/assess_assembly/{wildcards.sample}/{wildcards.sample}.bed
            python workflow/scripts/assess_assembly.py \
                --fasta {input.fasta} \
                --bam {input.bam} \
                --bed {params.results_dir}/assess_assembly/{wildcards.sample}/{wildcards.sample}.bed \
                --sample {wildcards.sample} \
                --output {params.results_dir}/assess_assembly/{wildcards.sample} &> {log}
        fi
        """
