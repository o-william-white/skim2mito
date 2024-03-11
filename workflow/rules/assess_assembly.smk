rule assess_assembly:
    input:
        "results/assembled_sequence/{sample}.ok",
        "results/annotations/{sample}/{sample}.ok",
        "results/minimap/{sample}.ok",
    output:
        ok="results/assess_assembly/{sample}.ok",
    log:
        "logs/assess_assembly/{sample}.log",
    conda:
        "../envs/assess_assembly.yaml"
    shell:
        """
        FAS=$(echo results/assembled_sequence/{wildcards.sample}.fasta)
        if [ -e $FAS ]; then
            if [ $(grep -e "^>" -c $FAS) -eq 1 ] ; then
                echo Single sequence found in fasta > {log}
                python workflow/scripts/assess_assembly.py \
                    --fasta results/assembled_sequence/{wildcards.sample}.fasta \
                    --bam results/minimap/{wildcards.sample}.bam \
                    --bed results/annotations/{wildcards.sample}/result.bed \
                    --sample {wildcards.sample} \
                    --output results/assess_assembly/
            else
                echo More than one sequence found in fasta > {log}
                # mitos creates subdirectories for each contig
                # find bed files and cat
                find results/annotations/{wildcards.sample}/ -type f -name result.bed | while read line; do  cat $line; done > results/assess_assembly/{wildcards.sample}.bed

                python workflow/scripts/assess_assembly.py \
                    --fasta results/assembled_sequence/{wildcards.sample}.fasta \
                    --bam results/minimap/{wildcards.sample}.bam \
                    --bed results/assess_assembly/{wildcards.sample}.bed \
                    --sample {wildcards.sample} \
                    --output results/assess_assembly/
            fi
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """
