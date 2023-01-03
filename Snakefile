import pandas as pd

# set configfile
configfile: "config/config.yaml"

# configfile parameters
output_dir = config["output_dir"]
blast_db = config["blast_db"]
taxdump = config["taxdump"]
mitos_refseq = config["mitos_refseq"]
mitos_code = config["mitos_code"]
threads = config["threads"]

# read sample data
sample_data = pd.read_csv(config["samples"]).set_index("ID", drop=False)

# functions to get forward and reverse reads from sample data
def get_forward(wildcards):
    return sample_data.loc[wildcards.sample, "forward"]

def get_reverse(wildcards):
    return sample_data.loc[wildcards.sample, "reverse"]

def get_organelle_type(wildcards):
    return sample_data.loc[wildcards.sample, "organelle_type"]

def get_seed(wildcards):
    return sample_data.loc[wildcards.sample, "seed"]

def get_gene(wildcards):
    return sample_data.loc[wildcards.sample, "gene"]

# one rule to rule them all :)
rule all:
    input:
        output_dir+"/summary/table.txt"

# convert fastq to fasta
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
        org = get_organelle_type,
        seed = get_seed,
        gene = get_gene
    output:
        #output_dir+"/getorganelle/{sample}/get_org.log.txt"
        ok = output_dir+"/getorganelle/{sample}/getorganelle.ok"
    log:
        output_dir+"/logs/getorganelle/{sample}.log"
    conda:
        "envs/getorganelle.yaml"
    threads: threads
    shell:
        """
        get_organelle_from_reads.py \
            -1 {input.fwd} -2 {input.rev} \
            -o {output_dir}/getorganelle/{wildcards.sample} \
            -F {params.org} \
            -s {params.seed} \
            --genes {params.gene} \
            --reduce-reads-for-coverage inf --max-reads inf \
            -R 20 \
            --overwrite -t {threads} &> {log}
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
            sed 's/>/>{wildcards.sample};/g' $FAS1 | awk NF > {output_dir}/assembled_sequence/{wildcards.sample}.fasta
        elif [ "$(echo $FAS | tr ' ' '\\n' | wc -l)" -eq 1 ]; then
            echo One assembly produced for {wildcards.sample} > {log}
            sed 's/>/>{wildcards.sample};/g' $FAS  | awk NF > {output_dir}/assembled_sequence/{wildcards.sample}.fasta
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
            minimap2 -ax sr $FAS {input.fwd} {input.rev} 2> {log} | samtools sort -O BAM -o $OUT - 2> {log}
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
                --table-fields gc,length,{wildcards.sample}_cov,{wildcards.sample}_read_cov,bestsumorder_superkingdom,bestsumorder_kingdom,bestsumorder_phylum,bestsumorder_class,bestsumorder_order,bestsumorder_family,bestsumorder_species \
                {output_dir}/blobtools/{wildcards.sample} &> {log}
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """

# treats all assemblies as circular
rule mitos:
    input:
        output_dir+"/assembled_sequence/{sample}.ok"
    output:
        ok = output_dir+"/mitos/{sample}/{sample}.ok"
    log:
        output_dir+"/logs/mitos/{sample}.log"
    conda:
        "envs/mitos.yaml"
    shell:
        """
        FAS=$(echo {output_dir}/assembled_sequence/{wildcards.sample}.fasta)
        if [ -e $FAS ]; then
            runmitos.py \
                --input $FAS \
                --code {mitos_code} \
                --outdir {output_dir}/mitos/{wildcards.sample}/ \
                --refseqver {mitos_refseq} \
                --refdir . \
                --noplots &> {log}
        else
            echo No assembled sequence for {wildcards.sample} > {log}
        fi
        touch {output.ok}
        """

rule extract_protein_coding_genes:
    input:
        output_dir+"/mitos/{sample}/{sample}.ok"
    output:
        ok = output_dir+"/protein_coding_genes/{sample}.ok"
    log:
        output_dir+"/logs/protein_coding_genes/{sample}.log"
    conda:
        "envs/seqkit.yaml"
    shell:
        """
        # find selected path(s) fasta
        FAS=$(find {output_dir}/mitos/{wildcards.sample}/ -name *result.fas)
        BED=$(find {output_dir}/mitos/{wildcards.sample}/ -name *result.bed)
        # z option: true if length if string is zero.
        if [[ -z $FAS ]]; then
           echo No annotated genes produced for {wildcards.sample} > {log}
        else
            # extract annotated protein coding genes in fasta format
            cat $FAS | sed -e 's/ //g' > results/protein_coding_genes/tmp_{wildcards.sample}.fasta
            grep -e "^>" results/protein_coding_genes/tmp_{wildcards.sample}.fasta | grep -e ";OH" -e ";OL" -e ";trn" -v | sed -e 's/>//g' > results/protein_coding_genes/tmp_{wildcards.sample}.txt
            if [[ $(cat results/protein_coding_genes/tmp_{wildcards.sample}.txt | wc -l) -eq 0 ]]; then
                echo No annotated protein coding genes produced for {wildcards.sample} > {log}
            else
                echo {wildcards.sample} genes
                seqkit grep -f results/protein_coding_genes/tmp_{wildcards.sample}.txt \
                    results/protein_coding_genes/tmp_{wildcards.sample}.fasta \
                    -o results/protein_coding_genes/{wildcards.sample}.fasta &> {log}
            fi
            rm results/protein_coding_genes/tmp_{wildcards.sample}.fasta
            rm results/protein_coding_genes/tmp_{wildcards.sample}.txt
            # extract annotated protein coding genes in bed format
            cat $BED | grep -P "\t[O|t]" -v > results/protein_coding_genes/{wildcards.sample}.bed
        fi
        touch {output.ok}
        """

rule summarise:
    input: 
        expand(output_dir+"/blobtools/{sample}/{sample}.ok", sample=sample_data["ID"].tolist()),
        expand(output_dir+"/mitos/{sample}/{sample}.ok", sample=sample_data["ID"].tolist())    
    output:
        table = output_dir+"/summary/table.txt"
    log:
        output_dir+"/logs/summarise/summarise.log"
    conda:
        "envs/r_env.yaml"
    shell:
        """ 
        Rscript scripts/summarise.R results/ {output.table} &> {log}
        """


