rule assess_assembly:
    input:
        get_sucessfully_annotated_samples,
        fasta="results/assembled_sequence/{sample}.fasta",
        bam="results/minimap/{sample}.bam",
    output:
        directory("results/assess_assembly/{sample}"),
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
                --bed results/annotations/{wildcards.sample}/result.bed \
                --sample {wildcards.sample} \
                --output results/assess_assembly/{wildcards.sample} &> {log}
        else
            echo More than one sequence found in fasta > {log}
            # mitos creates subdirectories for each contig
            # find bed files and cat
            mkdir -p results/assess_assembly/{wildcards.sample}
            find results/annotations/{wildcards.sample}/ -type f -name result.bed | while read line; do 
                cat $line
            done > results/assess_assembly/{wildcards.sample}/{wildcards.sample}.bed
            python workflow/scripts/assess_assembly.py \
                --fasta {input.fasta} \
                --bam {input.bam} \
                --bed results/assess_assembly/{wildcards.sample}/{wildcards.sample}.bed \
                --sample {wildcards.sample} \
                --output results/assess_assembly/{wildcards.sample} &> {log}
        fi
        """
