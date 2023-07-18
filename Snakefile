import pandas as pd

# set configfile
configfile: "config/config.yaml"

# configfile parameters
target_type = config["target_type"]
output_dir = config["output_dir"]
blast_db = config["blast_db"]
taxdump = config["taxdump"]
mitos_refseq = config["mitos_refseq"]
mitos_code = config["mitos_code"]
barrnap_kingdom = config["barrnap_kingdom"]
alignment_trim = config["alignment_trim"]
threads = config["threads"]

# read sample data
sample_data = pd.read_csv(config["samples"]).set_index("ID", drop=False)

# functions to get forward and reverse reads from sample data
def get_forward(wildcards):
    return sample_data.loc[wildcards.sample, "forward"]

def get_reverse(wildcards):
    return sample_data.loc[wildcards.sample, "reverse"]

def get_seed(wildcards):
    return sample_data.loc[wildcards.sample, "seed"]

def get_gene(wildcards):
    return sample_data.loc[wildcards.sample, "gene"]

# one rule to rule them all :)
rule all:
    input:
        output_dir+"/summary/summary_sample.txt",
        output_dir+"/summary/summary_contig.txt",
        output_dir+"/snakemake.ok",
        expand(output_dir+"/fastqc/{sample}_R1.html", 
            sample=sample_data.index.tolist())

# fastqc
# quality check fastq files
rule fastqc:
    input:
        fwd = get_forward,
        rev = get_reverse
    output:
        fwd = output_dir+"/fastqc/{sample}_R1.html",
        rev = output_dir+"/fastqc/{sample}_R2.html",
    params:
        outdir=lambda wildcards, output: os.path.abspath(os.path.dirname(output[0])) + "/",
        fwd_outfile = lambda wildcards, input: os.path.basename(input[0]).replace('.fastq.gz', '_fastqc.html').replace('.fq.gz','_fastqc.html'),
        rev_outfile = lambda wildcards, input: os.path.basename(input[1]).replace('.fastq.gz', '_fastqc.html').replace('.fq.gz','_fastqc.html')
    log:
        output_dir+"/logs/fastqc/{sample}.log"
    conda:
        "envs/fastqc.yaml"
    threads: 1
    shell:
        """
        fastqc -o "{params.outdir}" {input.fwd} -t {threads} &> {log} &&
        mv {params.outdir}/{params.fwd_outfile} {output.fwd} &&
        fastqc -o "{params.outdir}" {input.rev} -t {threads} &> {log} &&
        mv {params.outdir}/{params.rev_outfile} {output.rev}
        """

# rule fastqc:
#     input:
#         get_forward
#     output:
#         html=output_dir+"/fastqc/{sample}.html",
#         zip=output_dir+"/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
#     params:
#         extra = "--quiet"
#     log:
#         "logs/fastqc/{sample}.log"
#     threads: 1
#     resources:
#         mem_mb = 1024
#     wrapper:
#         "v2.1.1/bio/fastqc"

# convert fastq to fasta
# add option for forward and reverse primers if known
rule fastp:
    input:
        fwd = get_forward,
        rev = get_reverse
    output:
        fwd = output_dir+"/fastp/{sample}_R1.fq.gz",
        rev = output_dir+"/fastp/{sample}_R2.fq.gz",
        html = output_dir+"/fastp/{sample}.html",
        json = output_dir+"/fastp/{sample}.json"
    log:
        output_dir+"/logs/fastp/{sample}.log"
    conda:
        "envs/fastp.yaml"
    threads: threads
    shell:
        """
        fastp --in1 {input.fwd} --in2 {input.rev} \
            --out1 {output.fwd} --out2 {output.rev} \
            --html {output.html} --json {output.json} \
            --disable_quality_filtering \
            --thread {threads} &> {log}
        """

rule getorganelle:
    input:
        fwd = output_dir+"/fastp/{sample}_R1.fq.gz",
        rev = output_dir+"/fastp/{sample}_R2.fq.gz"
    params:
        seed = get_seed,
        gene = get_gene
    output:
        ok = output_dir+"/getorganelle/{sample}/getorganelle.ok"
    log:
        output_dir+"/logs/getorganelle/{sample}.log"
    conda:
        "envs/getorganelle.yaml"
    threads: threads
    shell:
        """
        if [[ {target_type} == "animal_mt" || {target_type} == "embplant_pt" ]]; then 
            get_organelle_from_reads.py \
                -1 {input.fwd} -2 {input.rev} \
                -o {output_dir}/getorganelle/{wildcards.sample} \
                -F {target_type} \
                -s {params.seed} \
                --genes {params.gene} \
                --reduce-reads-for-coverage inf --max-reads inf \
                -R 20 \
                --overwrite -t {threads} &> {log}
        else 
            if [ {target_type} == "anonym" ]; then
                get_organelle_from_reads.py \
                    -1 {input.fwd} -2 {input.rev} \
                    -o {output_dir}/getorganelle/{wildcards.sample} \
                    -F {target_type} \
                    -s {params.seed} \
                    --genes {params.gene} \
                    --reduce-reads-for-coverage inf --max-reads inf \
                    -R 10 \
                    --max-extending-len 100 \
                    -P 0 \
                    --overwrite -t {threads} &> {log}
            fi
        fi
        touch {output.ok}
        """

rule assembled_sequence:
    input:
        output_dir+"/getorganelle/{sample}/getorganelle.ok"
    output:
        ok = output_dir+"/assembled_sequence/{sample}.ok"
    log:
        output_dir+"/logs/assembled_sequence/{sample}.log"
    shell:
        """
        # find selected path(s) fasta
        FAS=$(find {output_dir}/getorganelle/{wildcards.sample}/ -name *path_sequence.fasta)
        # z option: true if length if string is zero.
        if [[ -z $FAS ]]; then
            echo No assembly produced for {wildcards.sample} > {log}
        # more than one selected path
        elif [ "$(echo $FAS | tr ' ' '\\n' | wc -l)" -gt 1 ]; then
            FAS1=$(echo $FAS | tr ' ' '\\n' | head -n 1)
            echo More than one assembly produced for {wildcards.sample} > {log}
            echo Selecting the first assembly $FAS1 > {log}
            python scripts/rename_assembled.py \
                --input $FAS1 \
                --sample {wildcards.sample} \
                --output {output_dir}/assembled_sequence
        elif [ "$(echo $FAS | tr ' ' '\\n' | wc -l)" -eq 1 ]; then
            echo One assembly produced for {wildcards.sample} > {log}
            python scripts/rename_assembled.py \
                --input $FAS \
                --sample {wildcards.sample} \
                --output {output_dir}/assembled_sequence
        fi
        touch {output.ok}
        """

rule seqkit:
    input:
        output_dir+"/assembled_sequence/{sample}.ok"
    output:
        ok = output_dir+"/seqkit/{sample}.ok"
    log:
        output_dir+"/logs/seqkit/{sample}.log"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        FAS=$(echo {output_dir}/assembled_sequence/{wildcards.sample}.fasta)
        OUT=$(echo {output_dir}/seqkit/{wildcards.sample}.txt)
        if [ -e $FAS ]; then
            echo Running seqkit for {wildcards.sample} > {log}
            seqkit stats -b $FAS > $OUT
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """

rule blastn:
    input: 
        output_dir+"/assembled_sequence/{sample}.ok"
    output:
        ok = output_dir+"/blastn/{sample}.ok"
    log:
        output_dir+"/logs/blastn/{sample}.log"
    conda:
        "envs/blastn.yaml"
    shell:
        """
        FAS=$(echo {output_dir}/assembled_sequence/{wildcards.sample}.fasta)
        OUT=$(echo {output_dir}/blastn/{wildcards.sample}.txt)
        if [ -e $FAS ]; then
            echo Running blastn for {wildcards.sample} > {log}
            blastn \
                -query $FAS \
                -db {blast_db} \
                -out $OUT \
                -outfmt '6 qseqid staxids bitscore std' \
                -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 &> {log} 
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """

rule minimap:
    input:
        output_dir+"/assembled_sequence/{sample}.ok",
        fwd = output_dir+"/fastp/{sample}_R1.fq.gz",
        rev = output_dir+"/fastp/{sample}_R2.fq.gz"
    output:
        ok = output_dir+"/minimap/{sample}.ok"
    log:
        output_dir+"/logs/minimap/{sample}.log"
    conda:
        "envs/minimap2.yaml"
    shell:
        """
        FAS=$(echo {output_dir}/assembled_sequence/{wildcards.sample}.fasta)
        OUT=$(echo {output_dir}/minimap/{wildcards.sample}.bam)
        if [ -e $FAS ]; then
            echo Running minimap for {wildcards.sample} > {log}
            minimap2 -ax sr $FAS {input.fwd} {input.rev} 2> {log} | samtools sort -O BAM -o $OUT - 2>> {log}
            samtools index $OUT 2>> {log}
            samtools index -c $OUT 2>> {log}
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """

# need to have taxdump in same dir 
rule blobtools:
    input:
        output_dir+"/assembled_sequence/{sample}.ok",
        output_dir+"/blastn/{sample}.ok",
        output_dir+"/minimap/{sample}.ok"
    output:
        ok = output_dir+"/blobtools/{sample}/{sample}.ok"
    log:
        output_dir+"/logs/blobtools/{sample}.log"
    container:
        # "docker://genomehubs/blobtoolkit"
        "docker://genomehubs/blobtoolkit"
    shell:
        """
        FAS=$(echo {output_dir}/assembled_sequence/{wildcards.sample}.fasta)
        BLA=$(echo {output_dir}/blastn/{wildcards.sample}.txt)
        MAP=$(echo {output_dir}/minimap/{wildcards.sample}.bam)
        OUT=$(echo {output_dir}/blobtools/{wildcards.sample}/table.tsv)
        if [ -e $FAS ]; then            
            blobtools create \
                --fasta $FAS \
                --hits $BLA \
                --taxrule bestsumorder \
                --taxdump {taxdump} \
                --cov $MAP \
                {output_dir}/blobtools/{wildcards.sample} &> {log}
            blobtools filter \
                --table $OUT \
                --table-fields gc,length,{wildcards.sample}_cov,bestsumorder_superkingdom,bestsumorder_kingdom,bestsumorder_phylum,bestsumorder_class,bestsumorder_order,bestsumorder_family,bestsumorder_species \
                {output_dir}/blobtools/{wildcards.sample} &>> {log}
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """

rule annotations:
    input:
        output_dir+"/assembled_sequence/{sample}.ok"
    output:
        ok = output_dir+"/annotations/{sample}/{sample}.ok"
    log:
        output_dir+"/logs/annotations/{sample}.log"
    conda:
        "envs/annotations.yaml"
    shell:
        """
        FAS=$(echo {output_dir}/assembled_sequence/{wildcards.sample}.fasta)
        if [ -e $FAS ]; then
            if [[ {target_type} == "animal_mt" ]]; then
                if [ $(grep circular -c $FAS) -eq 1 ] ; then 
                    echo Treating mitochondrial seqeunce as circular &> {log}
                    runmitos.py \
                        --input $FAS \
                        --code {mitos_code} \
                        --outdir {output_dir}/annotations/{wildcards.sample}/ \
                        --refseqver {mitos_refseq} \
                        --refdir . \
                        --noplots &>> {log} 
                else
                    echo Treating mitochndrial seqeunce as linear &> {log}
                    runmitos.py \
                        --input $FAS \
                        --code {mitos_code} \
                        --outdir {output_dir}/annotations/{wildcards.sample}/ \
                        --refseqver {mitos_refseq} \
                        --refdir . \
                        --noplots \
                        --linear &>> {log}
                fi
            else
                if [[ {target_type} == "anonym" ]]; then
                    barrnap \
                        --kingdom {barrnap_kingdom} \
                        --outseq {output_dir}/annotations/{wildcards.sample}/result.fas $FAS 1> {output_dir}/annotations/{wildcards.sample}/result.gff 2> {log}
                fi
            fi

        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """

rule assess_assembly:
    input:
        output_dir+"/assembled_sequence/{sample}.ok",
        output_dir+"/annotations/{sample}/{sample}.ok",
        output_dir+"/minimap/{sample}.ok"
    output:
        ok = output_dir+"/assess_assembly/{sample}.ok"
    log:
        output_dir+"/logs/assess_assembly/{sample}.log"
    conda:
        "envs/assess_assembly.yaml"
    shell:
        """
        FAS=$(echo {output_dir}/assembled_sequence/{wildcards.sample}.fasta)
        if [ -e $FAS ]; then
            if [[ {target_type} == "animal_mt" ]]; then
                if [ $(grep -e "^>" -c $FAS) -eq 1 ] ; then
                    echo Single sequence found in fasta > {log}
                    python scripts/assess_assembly.py \
                        --fasta {output_dir}/assembled_sequence/{wildcards.sample}.fasta \
                        --bam {output_dir}/minimap/{wildcards.sample}.bam \
                        --bed {output_dir}/annotations/{wildcards.sample}/result.bed \
                        --sample {wildcards.sample} \
                        --output {output_dir}/assess_assembly/
                else
                    echo More than one sequence found in fasta > {log}
                    # mitos creates subdirectories for each contig
                    # find bed files and cat
                    find {output_dir}/annotations/{wildcards.sample}/ -type f -name result.bed | while read line; do  cat $line; done > {output_dir}/assess_assembly/{wildcards.sample}.bed

                    python scripts/assess_assembly.py \
                        --fasta {output_dir}/assembled_sequence/{wildcards.sample}.fasta \
                        --bam {output_dir}/minimap/{wildcards.sample}.bam \
                        --bed {output_dir}/assess_assembly/{wildcards.sample}.bed \
                        --sample {wildcards.sample} \
                        --output {output_dir}/assess_assembly/
                fi
            else
                if [[ {target_type} == "anonym" ]]; then
                    echo TBC > {log}
                fi
            fi
        fi
        touch {output.ok}
        """

rule summarise:
    input: 
        expand(output_dir+"/seqkit/{sample}.ok", sample=sample_data["ID"].tolist()),
        expand(output_dir+"/blobtools/{sample}/{sample}.ok", sample=sample_data["ID"].tolist()),
        expand(output_dir+"/annotations/{sample}/{sample}.ok", sample=sample_data["ID"].tolist()),
        expand(output_dir+"/assess_assembly/{sample}.ok", sample=sample_data["ID"].tolist())
    output:
        table_sample = output_dir+"/summary/summary_sample.txt",
        table_contig = output_dir+"/summary/summary_contig.txt"
    log:
        output_dir+"/logs/summarise/summarise.log"
    conda:
        "envs/r_env.yaml"
    shell:
        """
        # cat seqkit output for each sample
        echo -e "sample format type num_seqs sum_len min_len avg_len max_len" > {output_dir}/summary/tmp_summary_sample.txt
        cat {output_dir}/seqkit/*.txt | grep file -v >> {output_dir}/summary/tmp_summary_sample.txt
        column -t {output_dir}/summary/tmp_summary_sample.txt > {output.table_sample}
        rm {output_dir}/summary/tmp_summary_sample.txt
        
        # join blobtools with mitos annotations for each contig
        if [[ {target_type} == "animal_mt" ]]; then
            Rscript scripts/summarise.R {output_dir}/ mitos {output.table_contig} &> {log}
        else
            if [[ {target_type} == "anonym" ]]; then
                Rscript scripts/summarise.R {output_dir}/ barrnap {output.table_contig} &> {log}
            fi
        fi
        """

checkpoint extract_protein_coding_genes:
    input: 
        expand(output_dir+"/annotations/{sample}/{sample}.ok", sample=sample_data["ID"].tolist())
    output:
        directory(output_dir+"/protein_coding_genes/")
    log:
        output_dir+"/logs/protein_coding_genes/protein_coding_genes.log"
    shell:
        """
        if [[ {target_type} == "animal_mt" ]]; then
            python scripts/mitos_alignments.py {output_dir}/annotations/ {output_dir}/protein_coding_genes &> {log}
        else
            if [[ {target_type} == "anonym" ]]; then
                python scripts/barrnap_alignments.py {output_dir}/annotations/ {output_dir}/protein_coding_genes &> {log}
            fi
        fi
        """

rule mafft:
    input:
        output_dir+"/protein_coding_genes/{dataset}.fasta"
    output:
        output_dir+"/mafft/{dataset}.fasta"
    log:
        output_dir+"/logs/mafft/{dataset}.log"
    conda:
        "envs/mafft.yaml"
    shell:
        """
        mafft \
            --maxiterate 1000 \
            --globalpair \
            --adjustdirection \
            {input} 1> {output} 2> {log}
        """

rule filter_alignments:
    input:
        output_dir+"/mafft/{dataset}.fasta"
    output:
        output_dir+"/mafft_filtered/{dataset}.fasta"
    log:
        output_dir+"/logs/mafft_filtered/{dataset}.log"
    shell:
        """
        python scripts/alignments_filter.py  --input {input} --output {output} --threshold 0.5 > {log}
        """

rule alignment_trim:
    input:
        output_dir+"/mafft_filtered/{dataset}.fasta"
    output:
        tmp = output_dir+"/alignment_trim/{dataset}_tmp.fasta",
        out = output_dir+"/alignment_trim/{dataset}.fasta"
    log:
        output_dir+"/logs/alignment_trim/{dataset}.log"
    conda:
        "envs/alignment_trim.yaml"
    shell:
        """
        if [ $(grep -c "^>" {input}) -lt "3" ]; then
            cp -R {input} {output.tmp}
            cp -R {input} {output.out}
        else
            # if [ $(grep -c "^>" {input[0]}) -lt "0" ]; then
            if [[ {alignment_trim} == "gblocks" ]]; then
                # gblocks add reuslts to same dir as input
                cp {input} {output.tmp}
                # gblocks always gives error code of 1. Ignore.
                Gblocks {output.tmp} -t=d &> {log} || true
                # sed to remove gaps
                sed 's/ //g' {output.tmp}-gb > {output.out}        
                # rm tmp
                # rm {output.tmp}-gb
            else
                if [[ {alignment_trim} == "clipkit" ]]; then
                    clipkit {input} -o {output.out} &> {log}
                fi
            fi
        fi
        """

rule iqtree:
    input:
        output_dir+"/alignment_trim/{dataset}.fasta"        
    output:
        output_dir+"/iqtree/{dataset}.treefile",
        renamed = output_dir+"/iqtree/{dataset}.fasta"
    log:
        output_dir+"/logs/iqtree/{dataset}.log"
    conda:
        "envs/iqtree.yaml"
    shell:
        """
        # only run iqtree if there are at least samples in alignment
        #if [ $(grep -e "^>" -c {input}) -ge 5 ] ; then  
        # remove special characters from sample names
        sed -e 's/;/_/g' -e 's/+//g' \
            {input} > {output.renamed}

        # iqtree
        # iqtree will not bootstrap if less than 5 samples in alignment
        # add if else statement here 
        #iqtree -s {output.renamed} --prefix {output_dir}/iqtree/{wildcards.dataset} &> {log}
        if [ $(grep -c "^>" {input}) -lt "3" ]; then
            touch {output[0]}
        else
            iqtree -s {output.renamed} -B 1000 --prefix {output_dir}/iqtree/{wildcards.dataset} &> {log}
        fi
        """

rule plot_tree:
    input:
        output_dir+"/iqtree/{dataset}.treefile"
    output:
        output_dir+"/plot_tree/{dataset}.png"
    log:
        output_dir+"/logs/plot_tree/{dataset}.log"
    conda:
        "envs/r_env.yaml"
    shell:
        """
        if [ $(grep -cvP '\S' {input[0]}) -eq "0" ]; then
            touch {output}
        else
            Rscript scripts/plot_tree.R {input} {output} &> {log}
        fi
        """

def get_plot_tree_output(wildcards):
    checkpoint_output = checkpoints.extract_protein_coding_genes.get(**wildcards).output[0]
    return expand(output_dir+"/plot_tree/{i}.png", i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i)

# create final log when complete 
rule final_log:
    input:
        get_plot_tree_output
    output:
        output_dir+"/snakemake.ok"
    shell:
        """
        touch {output}
        """


