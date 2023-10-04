import argparse
import pandas as pd
import pysam
import re
import statistics
import matplotlib.pyplot as plt

# argparse
parser = argparse.ArgumentParser()
parser.add_argument("--fasta",     help = "Input fasta",                  required=True)
parser.add_argument("--bam",       help = "Input BAM",                    required=True)
parser.add_argument("--bed",       help = "Input bed",                    required=True)
parser.add_argument("--sample",    help = "Sample name to use as prefix", required=True)
parser.add_argument("--output",    help = "Output directory",             required=True)
args = parser.parse_args()

# set additional params
bed_interval = 100

### functions

# create a bed file as a pandas dataframe
def create_bed_file(ref_name, ref_length, interval_size):
    bed_data = []
    for i in range(0, ref_length, interval_size):
        chrom_start = i
        chrom_end = min(i + interval_size, ref_length)
        bed_data.append([ref_name, chrom_start, chrom_end])
    bed_df = pd.DataFrame(bed_data, columns=["chrom", "start", "end"])
    return bed_df

# calculate mean depth within a region of interest in a bam file
def mean_depth(bam, chrom, start, end):
    depth = [0] * (end - start)
    for pileupcolumn in bam.pileup(chrom, start, end):
        pos = pileupcolumn.pos - start
        if pos >= 0 and pos < len(depth):
            depth[pos] = pileupcolumn.nsegments
    return statistics.mean(depth)

# count the proportion of mismatched bases within a region of interest in a bam file
def proportion_mismatches(bam, fas, chrom, start, end):
    matches = 0
    mismatches = 0
    for pileupcolumn in bam.pileup(chrom, start, end):
        for pileupread in pileupcolumn.pileups:
            if pileupread.query_position is None:
                continue
            query_base = pileupread.alignment.query_sequence[pileupread.query_position]
            reference_base = fas.fetch(reference=pileupcolumn.reference_name,
                    start=pileupcolumn.reference_pos,
                    end=pileupcolumn.reference_pos + 1).upper()
            if query_base == reference_base:
                matches += 1
            else:
                mismatches += 1
    total_bases = matches + mismatches
    if total_bases == 0:
        return 0.0
    else:
        return float(mismatches) / float(total_bases)

# calculate gc content within a region of interest in a fasta file
def gc_content(fas, chrom, start, end):
    gc_count = 0
    total_count = 0
    seq = fas.fetch(chrom, start, end)
    for base in seq:
        if base.upper() == "G" or base.upper() == "C":
            gc_count += 1
        total_count += 1
    if total_count == 0:
        return 0.0
    else:
        return float(gc_count) / float(total_count)

### main

# open input files
fas = pysam.FastaFile(args.fasta)
bam = pysam.AlignmentFile(filename = args.bam)
# bed_annotations = open(args.bed, "r")

# iterate through bed_annotations file and append gene locations to list
#bed_annotations_list = []
#for line in bed_annotations:
#    line = line.rstrip("\n").split("\t")
#    chrom, start, end, annotation = line[0], line[1], line[2], line[3]
#    # protein coding genes only
#    if not annotation.startswith("trn"):
#        bed_annotations_list.append([chrom, start, end, annotation])

# iterate through nreferences
for n in range(0, bam.nreferences):

    # get reference name and lengh
    reference_name = bam.references[n]
    reference_length = bam.lengths[n]

    # create bed file as pandas df
    bed_df = create_bed_file(ref_name = bam.references[n], ref_length = bam.lengths[n], interval_size = bed_interval)
    
    # iterate through bed_annotations file and append gene locations to list
   
    with open(args.bed, "r") as bed_annotations:
        bed_annotations_list = []
        for line in bed_annotations:
            line = line.rstrip("\n").split("\t")
            chrom, start, end, annotation = line[0], line[1], line[2], line[3]
            # protein coding genes only
            if chrom == reference_name and not re.search("^trn|^OH" , annotation):
                bed_annotations_list.append([chrom, start, end, annotation])
                #print(f"{chrom}, {start}, {end}, {annotation}")
    #print(reference_name)
    #print(reference_length)
    #print(bed_df)
    #print(bed_annotations)

    # create emtpy lists for mean depth, gc content and proportion of mismatches
    md = []
    pm = []
    gc = []

    # iterate through rows in bed file
    for index, row in bed_df.iterrows():
        # get values from each row
        chrom, start, end = row['chrom'], row['start'], row['end']
        # append values to list
        md.append(mean_depth(bam, chrom, start, end))
        pm.append(proportion_mismatches(bam, fas, chrom, start, end))
        gc.append(gc_content(fas, chrom, start, end))

    # add lists to dataframe
    bed_df["mean_depth"] = md
    bed_df["gc_content"] = gc
    bed_df["proportion_mismatches"] = pm

    # write .txt file
    bed_df.to_csv(f"{args.output}/{args.sample}_{n}.txt", sep="\t")

    # plot
    fig, axes = plt.subplots(nrows = 4, ncols = 1, figsize=(8, 8))

    # plot each gene as a horizontal line
    contig_annotations = 0
    for annotation in bed_annotations_list:
        contig = annotation[0]
        start = annotation[1]
        end = annotation[2]
        annotation_name = annotation[3]
        if contig == reference_name: 
            contig_annotations =+ 1
            axes[0].hlines(y=annotation_name, xmin=int(start), xmax=int(end), linewidth=2, color='gray')
            annotation_position = int(start) + ((int(end)-int(start))/2)
            axes[0].annotate(annotation_name, xy=(annotation_position, annotation_name), fontsize=8, ha='center', va='bottom')

    axes[0].set_ylim(-2, len(bed_annotations_list))
    axes[0].yaxis.set_ticklabels([])
    axes[0].yaxis.set_ticks([])
    axes[0].invert_yaxis()
    axes[0].set_ylabel('Annotations')

    axes[1].plot(bed_df['start'], bed_df['mean_depth'], color='tab:blue')
    axes[1].set_ylabel('Mean depth')

    axes[2].plot(bed_df['start'], bed_df['gc_content'], color='tab:orange')
    axes[2].set_ylabel('GC content')

    axes[3].plot(bed_df['start'], bed_df['proportion_mismatches'], color='tab:green')
    axes[3].set_ylabel('Proportion mismatches')

    # adjust spacing between subplots
    fig.tight_layout()

    # Save the plot to a file
    plt.savefig(f"{args.output}/{args.sample}_{n}.png")

# close files
bam.close()
fas.close()
bed_annotations.close()
