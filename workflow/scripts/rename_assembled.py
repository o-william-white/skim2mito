import os
import logging

# simple python script to rename assembled sequence
# a two column table is written with new and old names

logger = logging.getLogger('logging')
fh = logging.FileHandler(str(snakemake.log[0]))
fh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)


def read_fasta(filename):
    name, seq = None, ""
    with open(filename) as fasta:
        for line in fasta:
            if line.startswith('>') and name is None:
                name = line.rstrip('\n').replace('>', '')
            else:
                if line.startswith('>') and name is not None:
                    yield name, seq
                    name = line.rstrip('\n').replace('>', '')
                    seq = ''
                else:
                    seq = seq + line.rstrip('\n')
        yield name, seq


# function to write a renamed fasta and summary of name changes
def rename_fasta(fasta, sample_name, output_prefix):
    # open output files
    with (open(f"{output_prefix}.fasta", "w") as new_fasta,
          open(os.path.join(snakemake.output[0], "rename.txt"), "a+") as name_file):
        # iterate through fasta
        for contig_number, (name, sequence) in enumerate(fasta):
            if "circular" in name:
                # define short name if sequence circular
                short_name = f"{sample_name}_circular"
            else:
                # define short name if contig
                short_name = f"{sample_name}_contig{contig_number}"
            # write to output files
            new_fasta.write(f">{short_name}\n{sequence}\n")
            name_file.write(f"{short_name}\t{name}\n")


try:
    # make output dir if not already present
    if not os.path.exists(snakemake.output[0]):
        os.mkdir(snakemake.output[0])

    for getorganelle_output in snakemake.input:
        sample = os.path.basename(getorganelle_output)
        path_sequences = [os.path.join(getorganelle_output, f)
                          for f in os.listdir(getorganelle_output) if "path_sequence" in f]
        if len(path_sequences) != 0:
            # read fasta
            fas = read_fasta(path_sequences[0])  # only read the first file
            # rename and write fasta and summary file
            rename_fasta(fas, sample, os.path.join(snakemake.output[0], sample))

except Exception as e:
    print(e)
    logger.error(e, exc_info=True)
    exit(1)
